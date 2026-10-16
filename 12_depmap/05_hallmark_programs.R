# VPS72 expression and copy number against all 50 MSigDB Hallmark programmes,
# pan-cancer and per lineage, plus the quartile dose-response tables.
source("00_config.R")
suppressPackageStartupMessages({library(msigdbr)})

model <- read_model()

hall <- as.data.table(msigdbr(species="Homo sapiens",category="H"))
if (!nrow(hall)) stop("msigdbr returned no Hallmark sets")
sets <- split(toupper(hall$gene_symbol),hall$gs_name)
cat("Hallmark sets:",length(sets),"\n")

# Matrix columns are "SYMBOL (ENTREZ)". Symbols that collapse to the same name
# would silently take whichever column came first, so duplicates are dropped.
load_matrix <- function(p) {
  m <- fread(require_input(p))
  setnames(m,1,"ModelID")
  syms <- strip_entrez(names(m)[-1L])
  bad <- is.na(syms) | !nzchar(syms) | duplicated(syms)
  if (any(bad)) cat("  dropping",sum(bad),"unnamed or duplicate gene columns\n")
  keep <- c("ModelID",names(m)[-1L][!bad])
  m <- m[,..keep]
  setnames(m,c("ModelID",syms[!bad]))
  m
}

expr <- load_matrix(inputs$expr_matrix)
if (!"VPS72" %chin% names(expr)) stop("VPS72 is not a column of ",inputs$expr_matrix)
cat("expression matrix:",nrow(expr),"lines,",ncol(expr)-1L,"genes\n")

# A programme score is the mean of its member genes' z-scores across cell
# lines, so each gene contributes equally regardless of expression level.
score_set <- function(nm) {
  present <- intersect(sets[[nm]],names(expr))
  if (length(present)<min_set_n) return(NULL)
  z <- scale(as.matrix(expr[,..present]))
  data.table(ModelID=expr$ModelID,program=nm,score=rowMeans(z,na.rm=TRUE),
             n_genes=length(present))
}
scored <- rbindlist(lapply(names(sets),score_set))
if (!nrow(scored)) stop("No Hallmark set reached ",min_set_n," member genes")
cat("programmes scored:",uniqueN(scored$program),"\n")

master <- dcast(scored,ModelID~program,value.var="score")
master <- merge(master,data.table(ModelID=expr$ModelID,vps72_expr=expr$VPS72),by="ModelID")
master <- merge(master,model[,.(ModelID,lineage,is_cancer)],by="ModelID")
master <- master[!is.na(lineage) & nzchar(lineage) & lineage!="Unknown"]
if (!nrow(master)) stop("No scored cell line carried a lineage label")

# Copy number is optional; the WGS matrix is large and not every release has it.
cnv_ok <- FALSE
if (file.exists(inputs$cnv_matrix)) {
  cnv <- load_matrix(inputs$cnv_matrix)
  if ("VPS72" %chin% names(cnv)) {
    master <- merge(master,data.table(ModelID=cnv$ModelID,vps72_cn=cnv$VPS72),
                    by="ModelID",all.x=TRUE)
    cnv_ok <- sum(is.finite(master$vps72_cn))>=min_lineage_n
    cat("VPS72 copy number on",sum(is.finite(master$vps72_cn)),"lines\n")
  }
  rm(cnv); gc()
} else cat("no copy-number matrix at",inputs$cnv_matrix,"- copy-number arm skipped\n")

programs <- sort(unique(scored$program))
rm(expr,scored); gc()

## ---- Pan-cancer ------------------------------------------------------------
pan_one <- function(prog,xvar) {
  s <- spearman(master[[xvar]],master[[prog]])
  if (is.null(s)) NULL else { s[,`:=`(program=prog,predictor=xvar)]; s[] }
}
pan <- rbindlist(lapply(programs,pan_one,xvar="vps72_expr"))
if (cnv_ok) pan <- rbind(pan,rbindlist(lapply(programs,pan_one,xvar="vps72_cn")))
if (!nrow(pan)) stop("No programme could be correlated pan-cancer")
pan[,FDR:=p.adjust(p,method="BH"),by=predictor]
setorder(pan,predictor,p)
write_tsv(pan,"hallmark_pancancer.tsv")
print(pan[predictor=="vps72_expr",.(n_programs=.N,n_sig=sum(FDR<fdr_cut)),by=predictor])

## ---- Per lineage -----------------------------------------------------------
# The lineage label is attached inside the loop. Assigning it afterwards with
# rep(lineages, each = n_programs) is wrong whenever a lineage drops a
# programme for want of data, which silently shifts every label after it.
lineages <- master[,.N,by=lineage][N>=min_lineage_n]$lineage
if (!length(lineages)) stop("No lineage reaches ",min_lineage_n," cell lines")

lin <- rbindlist(lapply(lineages,function(l) {
  d <- master[lineage==l]
  rbindlist(lapply(programs,function(prog) {
    s <- spearman(d[[prog]],d$vps72_expr,min_n=min_lineage_n)
    if (is.null(s)) NULL else { s[,`:=`(lineage=l,program=prog)]; s[] }
  }))
}))
if (!nrow(lin)) stop("No lineage-programme pair had enough paired values")
lin[,FDR:=p.adjust(p,method="BH"),by=lineage]
setorder(lin,lineage,p)
write_tsv(lin,"hallmark_by_lineage.tsv")
cat("Lineages:",uniqueN(lin$lineage),"; lineage-programme pairs:",nrow(lin),
    "; at FDR <",fdr_cut,":",sum(lin$FDR<fdr_cut,na.rm=TRUE),"\n")

# Wide rho matrix, the table the lineage heatmap was drawn from.
write_tsv(dcast(lin,program~lineage,value.var="rho"),"hallmark_rho_matrix.tsv")

## ---- Quartile dose-response ------------------------------------------------
# Medians per quartile with a Kruskal-Wallis test across the four groups. The
# continuous Spearman above is the primary estimate; this is the shape of it.
quartile_table <- function(xvar,label) {
  d <- master[is.finite(get(xvar))]
  if (nrow(d)<40L) return(NULL)
  br <- quantile(d[[xvar]],probs=seq(0,1,0.25),na.rm=TRUE)
  if (anyDuplicated(br)) { cat("quartile breaks not unique for",label,"- skipped\n"); return(NULL) }
  d[,quartile:=cut(get(xvar),breaks=br,labels=c("Q1","Q2","Q3","Q4"),include.lowest=TRUE)]
  rbindlist(lapply(programs,function(prog) {
    z <- d[is.finite(get(prog)),.(quartile,score=get(prog))]
    if (uniqueN(z$quartile)<4L) return(NULL)
    kw <- suppressWarnings(kruskal.test(score~quartile,data=z))
    med <- dcast(z[,.(median=median(score)),by=quartile],.~quartile,value.var="median")
    med[,`.`:=NULL]
    cbind(data.table(predictor=label,program=prog,n=nrow(z)),med,
          data.table(kruskal_p=kw$p.value))
  }))
}

q_expr <- quartile_table("vps72_expr","vps72_expr")
q_all <- if (cnv_ok) rbind(q_expr,quartile_table("vps72_cn","vps72_cn"),fill=TRUE) else q_expr
if (!is.null(q_all) && nrow(q_all)) {
  q_all[,kruskal_FDR:=p.adjust(kruskal_p,method="BH"),by=predictor]
  setorder(q_all,predictor,kruskal_p)
  write_tsv(q_all,"hallmark_quartile_dose_response.tsv")
}

## ---- Expression against copy number ----------------------------------------
# Whether a programme tracking VPS72 expression also tracks VPS72 dosage. If it
# does, the association is more easily read as dosage-driven.
if (cnv_ok) {
  cmp <- merge(pan[predictor=="vps72_expr",.(program,rho_expr=rho,FDR_expr=FDR)],
               pan[predictor=="vps72_cn",.(program,rho_cn=rho,FDR_cn=FDR)],by="program")
  cmp[,`:=`(sig_expr=FDR_expr<fdr_cut,sig_cn=FDR_cn<fdr_cut)]
  cmp[,agreement:=fcase(sig_expr & sig_cn & sign(rho_expr)==sign(rho_cn),"both_same_direction",
                        sig_expr & sig_cn,"both_opposite_direction",
                        sig_expr,"expression_only",
                        sig_cn,"copy_number_only",
                        default="neither")]
  setorder(cmp,-abs(rho_expr))
  write_tsv(cmp,"hallmark_expression_vs_copynumber.tsv")
  r <- spearman(cmp$rho_expr,cmp$rho_cn)
  if (!is.null(r)) cat("Expression rho vs copy-number rho across programmes: rho =",
                       round(r$rho,3),"over",r$n,"programmes\n")
  print(cmp[,.N,by=agreement])
}

write_tsv(master[,.(ModelID,lineage,is_cancer,vps72_expr,
                    vps72_cn=if (cnv_ok) vps72_cn else NA_real_)],
          "hallmark_program_scores.tsv")
