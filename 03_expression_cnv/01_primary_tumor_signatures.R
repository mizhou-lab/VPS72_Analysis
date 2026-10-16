# Primary TCGA tumors: aligned expression/GISTIC signature scores and associations.
# Outputs numeric data and statistical results only.
source("00_config.R")

get_signatures <- function() {
  out <- lapply(signature_files,function(x) {
    fp <- require_input(file.path(input_dir,x))
    dt <- fread(fp,sep="\t",header=FALSE)
    if (ncol(dt)<2L) stop("Malformed gene-set file: ",fp)
    v <- dt[V1=="GENE_SYMBOLS",V2]
    if (length(v)!=1L) stop("Expected one GENE_SYMBOLS row: ",fp)
    unique(toupper(trimws(unlist(strsplit(v,",")))))
  })
  if (any(lengths(out)<3L)) stop("At least one signature contains fewer than 3 genes")
  out
}

sigs <- get_signatures()
model_genes <- unique(c("VPS72",srcap_genes,tip60_genes))
required <- unique(c(unlist(sigs,use.names=FALSE),model_genes))
expr <- load_gene_matrix(inputs$tcga_expression,required,"Expression")
cnv <- load_gene_matrix(inputs$tcga_gistic,required,"GISTIC")
expr <- expr[,keep_primary(colnames(expr)),drop=FALSE]
cnv  <- cnv[,keep_primary(colnames(cnv)),drop=FALSE]
expr <- collapse_samples(expr); cnv <- collapse_samples(cnv)
if (!"VPS72" %in% rownames(expr) || !"VPS72" %in% rownames(cnv))
  stop("VPS72 missing from expression or GISTIC")

clinical <- fread(require_input(inputs$tcga_clinical))
find_col <- function(pat) { x <- grep(pat,names(clinical),ignore.case=TRUE,value=TRUE); if (!length(x)) stop("Clinical column missing: ",pat); x[1L] }
id_col <- find_col("^Sample[ _]?ID$")
study_col <- find_col("^Study[ _]?ID$")
clinical_map <- unique(data.table(patient=tcga12(clinical[[id_col]]),
                                  cohort=toupper(trimws(as.character(clinical[[study_col]])))))
clinical_map <- clinical_map[!is.na(patient) & !is.na(cohort) & nzchar(cohort)]
if (clinical_map[,uniqueN(cohort),by=patient][V1>1L,.N]>0L)
  stop("Clinical metadata assigns multiple cancer cohorts to one patient")
clinical_map <- unique(clinical_map,by="patient")
common <- intersect(colnames(expr),colnames(cnv))
common <- common[substr(common,1L,12L) %in% clinical_map$patient]
if (!length(common)) stop("No common primary-tumor samples")
expr <- expr[,common,drop=FALSE]; cnv <- cnv[,common,drop=FALSE]
samples <- data.table(sample=common,patient=substr(common,1L,12L))
samples[,cohort:=clinical_map$cohort[match(patient,clinical_map$patient)]]
stopifnot(!anyNA(samples$cohort),!anyDuplicated(samples$sample))
write_tsv(samples,"primary_tumor_manifest.tsv")

zrows <- function(m) {
  s <- matrixStats::rowSds(m,na.rm=TRUE)
  s[!is.finite(s) | s==0] <- NA_real_
  sweep(sweep(m,1L,rowMeans(m,na.rm=TRUE),"-"),1L,s,"/")
}
calc_signature <- function(mat, genes, method=c("zmean","mean")) {
  method <- match.arg(method)
  g <- intersect(genes,rownames(mat))
  if (length(g)<3L) return(rep(NA_real_,ncol(mat)))
  m <- mat[g,,drop=FALSE]
  if (method=="zmean") m <- zrows(m)
  res <- colMeans(m,na.rm=TRUE)
  res[!is.finite(res)] <- NA_real_
  res
}

expr_scores <- data.table(sample=common)
cnv_scores <- data.table(sample=common)
coverage <- rbindlist(lapply(names(sigs),function(nm) data.table(
  signature=nm,genes_in_set=length(sigs[[nm]]),
  expression_found=length(intersect(sigs[[nm]],rownames(expr))),
  gistic_found=length(intersect(sigs[[nm]],rownames(cnv)))
)))
for (nm in names(sigs)) {
  expr_scores[,(nm):=calc_signature(expr,sigs[[nm]],"zmean")]
  cnv_scores[,(nm):=calc_signature(cnv,sigs[[nm]],"mean")]
}
expr_scores[,VPS72:=as.numeric(expr["VPS72",])]
expr_scores <- merge(expr_scores,samples[,.(sample,cohort)],by="sample",sort=FALSE)
cnv_scores <- merge(cnv_scores,samples[,.(sample,cohort)],by="sample",sort=FALSE)
expr_scores[,group:= {
  med <- median(VPS72,na.rm=TRUE)
  ifelse(is.finite(VPS72),ifelse(VPS72>=med,"High","Low"),NA_character_)
},by=cohort]
expr_scores[,composite:=rowMeans(.SD,na.rm=TRUE),.SDcols=names(sigs)]
expr_scores[!is.finite(composite),composite:=NA_real_]
cnv_scores <- merge(cnv_scores,expr_scores[,.(sample,group,VPS72)],by="sample",sort=FALSE)
write_tsv(coverage,"signature_gene_coverage.tsv")
write_tsv(expr_scores,"primary_expression_signature_scores.tsv")
write_tsv(cnv_scores,"primary_gistic_signature_scores.tsv")

# Every cohort x signature tested with two-sided Wilcoxon; BH across that family.
stat_set <- c(names(sigs),"composite")
long_e <- melt(expr_scores,id.vars=c("sample","cohort","group","VPS72"),
               measure.vars=stat_set,variable.name="signature",value.name="score")
res <- long_e[,.(n_high=sum(is.finite(score) & group=="High",na.rm=TRUE),
                n_low=sum(is.finite(score) & group=="Low",na.rm=TRUE),
                median_high=if (any(is.finite(score) & group=="High",na.rm=TRUE)) median(score[group=="High"],na.rm=TRUE) else NA_real_,
                median_low=if (any(is.finite(score) & group=="Low",na.rm=TRUE)) median(score[group=="Low"],na.rm=TRUE) else NA_real_,
                p_value=safe_wilcox(score,group)),by=.(cohort,signature)]
res[,delta_median:=median_high-median_low]
res[,FDR:=p.adjust(p_value,method="BH")]
write_tsv(res,"per_cancer_signature_high_low.tsv")

# Continuous association avoids discarding intermediate VPS72 values.
continuous <- long_e[,as.list(safe_spearman(VPS72,score,min_n=min_group_n)),by=.(cohort,signature)]
continuous[,FDR:=p.adjust(p,method="BH")]
write_tsv(continuous,"per_cancer_signature_spearman.tsv")

# Source for residual models: small, aligned genuine sample-level matrices.
saveRDS(expr[intersect(model_genes,rownames(expr)),,drop=FALSE],
        file.path(output_dir,"model_expression_samples.rds"))
saveRDS(cnv[intersect("VPS72",rownames(cnv)),,drop=FALSE],
        file.path(output_dir,"model_gistic_samples.rds"))
cat("Primary samples:",length(common),"; cancer cohorts:",uniqueN(samples$cohort),"\n")
cat("Analyzed signatures:",paste(names(sigs),collapse=", "),"\n")
