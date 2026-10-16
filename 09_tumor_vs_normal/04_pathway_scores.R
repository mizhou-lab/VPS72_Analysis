# Per-sample Hallmark and composite proliferation scores, then the VPS72
# median-split comparison within each cancer and tissue state.
#
# RECONSTRUCTED STEP. The script that produced these tables was not present in
# the working folder; the notebook downstream read its outputs without ever
# writing them. What follows was rebuilt from the schema those outputs had to
# satisfy, so it reproduces the structure of the published tables but not
# necessarily the exact scoring method used. Check it against your own record
# before relying on the values. The two choices to verify are marked below.
source("00_config.R")

master <- fread(require_input(inputs$master))
tcga <- master[source=="TCGA" & cancer %chin% eligible & state %chin% c("Normal","Tumor")]
if (!nrow(tcga)) stop("No eligible TCGA samples in the master table")

sets <- lapply(names(hallmark_files),
               function(nm) read_geneset(file.path(inputs$genesets,hallmark_files[[nm]])))
names(sets) <- names(hallmark_files)
members <- unique(unlist(sets,use.names=FALSE))

# Xena rows are Ensembl IDs; the gencode probeMap maps them to symbols when it
# sits beside the matrix.
expr <- fread(require_input(inputs$xena_expr))
setnames(expr,1,"gene_id")

sym_file <- file.path(dirname(inputs$xena_expr),"probeMap_gencode.v23.annotation.gene.probemap")
if (file.exists(sym_file)) {
  pm <- fread(sym_file)
  setnames(pm,1:2,c("gene_id","symbol"))
  pm[,symbol:=toupper(symbol)]
  expr <- merge(expr,pm[,.(gene_id,symbol)],by="gene_id")
} else {
  expr[,symbol:=toupper(sub("\\..*$","",gene_id))]
  cat("no probeMap at",sym_file,"- treating row IDs as symbols\n")
}

expr <- expr[symbol %chin% members]
if (nrow(expr)<50L) stop("Only ",nrow(expr)," member genes matched the expression matrix; ",
                         "check that the probeMap is present and that symbols map")

long <- melt(expr,id.vars=c("gene_id","symbol"),
             measure.vars=setdiff(names(expr),c("gene_id","symbol")),
             variable.name="sample",value.name="expr")
long[,sample:=tcga15(sample)]
long[,expr:=suppressWarnings(as.numeric(expr))]
long <- long[is.finite(expr) & sample %chin% tcga$sample]
rm(expr); gc()

# CHOICE 1 to verify: genes are z-scored within cancer, across tumour and
# normal together, matching the VPS72 standardization in 03. A set score is
# then the mean z of its member genes. Other scoring methods (GSVA, singscore)
# would give different absolute values and a different composite scale.
long <- merge(long,tcga[,.(sample,cancer)],by="sample")
long[,gene_z:=as.numeric(scale(expr)),by=.(symbol,cancer)]

score_one <- function(nm) {
  sub <- long[symbol %chin% sets[[nm]]]
  if (!nrow(sub)) stop("No member gene recovered for ",nm)
  s <- sub[,.(score=mean(gene_z,na.rm=TRUE)),by=sample]
  setnames(s,"score",nm)
  s
}
scores <- Reduce(function(a,b) merge(a,b,by="sample",all=TRUE),
                 lapply(names(sets),score_one))

# CHOICE 2 to verify: the composite is the unweighted mean of the five
# proliferation sets named in prolif_sets, not of all eight Hallmark sets.
prolif_ok <- intersect(prolif_sets,names(scores))
if (length(prolif_ok)<2L) stop("Fewer than two proliferation sets scored")
scores[,proliferation_score:=rowMeans(.SD,na.rm=TRUE),.SDcols=prolif_ok]

dat <- merge(tcga[,.(sample,patient,cancer,state,VPS72,VPS72_z,paired)],scores,by="sample")
if (!nrow(dat)) stop("No sample carried both VPS72 and a pathway score")

fp <- file.path(OUT_DIR,"pathway_scores_per_sample.tsv.gz")
fwrite(dat,fp,sep="\t")
cat("wrote",fp,"(",nrow(dat),"rows )\n")

# Median split on VPS72 within each cancer and state, so the comparison is
# always within one lineage and one tissue state.
score_cols <- c("proliferation_score",names(sets))
split_one <- function(pw,ca,st) {
  z <- dat[cancer==ca & state==st,.(VPS72,score=get(pw))]
  z <- z[is.finite(VPS72) & is.finite(score)]
  if (nrow(z)<min_cor_n) return(NULL)
  z[,grp:=fifelse(VPS72>=median(VPS72),"High","Low")]
  if (uniqueN(z$grp)<2L) return(NULL)
  hi <- z[grp=="High"]$score; lo <- z[grp=="Low"]$score
  if (!length(hi) || !length(lo)) return(NULL)
  wt <- suppressWarnings(wilcox.test(hi,lo))
  data.table(pathway=pw,Program=pretty_name(pw),cancer=ca,state=st,
             n_high=length(hi),n_low=length(lo),
             median_high=median(hi),median_low=median(lo),
             delta_high_low=mean(hi)-mean(lo),p=wt$p.value)
}

grid <- CJ(pw=score_cols,ca=eligible,st=c("Normal","Tumor"),sorted=FALSE)
res <- rbindlist(Map(split_one,grid$pw,grid$ca,grid$st))
if (!length(res) || !nrow(res)) stop("No cancer/state cell reached ",min_cor_n," samples")

# BH within pathway and state: each pathway in each state is its own family of
# tests across cancer types.
res[,FDR:=p.adjust(p,method="BH"),by=.(pathway,state)]
res[,positive:=FDR<fdr_cut & delta_high_low>0]
setorder(res,pathway,state,p)

write_tsv(res,"pathway_high_low_results.tsv")
cat("Pathways:",uniqueN(res$pathway),"; cancer/state cells:",nrow(res),
    "; positive at FDR <",fdr_cut,":",sum(res$positive,na.rm=TRUE),"\n")
print(res[pathway=="proliferation_score",.(state,positive=sum(positive,na.rm=TRUE),
                                           tested=.N),by=state])
