# VPS72 expression in cancer versus non-cancerous cell lines, and the
# copy-number to expression relationship.
source("00_config.R")

model <- read_model()
expr <- read_gene_export(inputs$expr_gene,"expr_log1p_tpm")
cnv  <- read_gene_export(inputs$cnv_gene,"cn_log2")

dat <- merge(model,expr,by="ModelID")
dat <- dat[is.finite(expr_log1p_tpm)]
if (!nrow(dat)) stop("No cell line carried a VPS72 expression value")

# The non-cancerous arm of DepMap is small, so the group sizes are reported
# with the test rather than left to the figure.
counts <- dat[,.(n=.N,median=median(expr_log1p_tpm),mean=mean(expr_log1p_tpm),
                 sd=sd(expr_log1p_tpm)),by=is_cancer]
if (uniqueN(dat$is_cancer)<2L) {
  cat("only one group present; cancer vs normal test skipped\n")
} else {
  wt <- wilcox.test(expr_log1p_tpm~is_cancer,data=dat)
  res <- data.table(comparison="Cancer_vs_Normal",
                    n_cancer=counts[is_cancer=="Cancer"]$n,
                    n_normal=counts[is_cancer=="Normal"]$n,
                    median_cancer=counts[is_cancer=="Cancer"]$median,
                    median_normal=counts[is_cancer=="Normal"]$median,
                    W=as.numeric(wt$statistic),p=wt$p.value)
  write_tsv(res,"expression_cancer_vs_normal.tsv")
  print(res)
}
write_tsv(counts,"expression_group_summary.tsv")

cn <- merge(dat,cnv,by="ModelID")
cn <- cn[is.finite(cn_log2)]
if (!nrow(cn)) stop("No cell line carried both expression and copy number")
write_tsv(cn[,.(ModelID,name,lineage,is_cancer,expr_log1p_tpm,cn_log2)],
          "expression_cnv_per_line.tsv")

r <- spearman(cn$cn_log2,cn$expr_log1p_tpm)
if (is.null(r)) stop("Too few lines to correlate copy number with expression")
r[,comparison:="CN_vs_expression"]
write_tsv(r,"cnv_expression_correlation.tsv")
cat("Copy number vs expression: rho =",round(r$rho,3),"over",r$n,"lines, p =",
    signif(r$p,3),"\n")

# Per lineage, since a pan-cancer dosage relationship can be driven by a few
# lineages with recurrent 1q gain.
by_lineage <- rbindlist(lapply(
  cn[,.N,by=lineage][N>=min_lineage_n]$lineage,
  function(l) {
    d <- cn[lineage==l]
    s <- spearman(d$cn_log2,d$expr_log1p_tpm,min_n=min_lineage_n)
    if (is.null(s)) NULL else { s[,lineage:=l]; s[] }
  }))
if (nrow(by_lineage)) {
  by_lineage[,FDR:=p.adjust(p,method="BH")]
  setorder(by_lineage,p)
  write_tsv(by_lineage,"cnv_expression_by_lineage.tsv")
  cat("Lineages tested:",nrow(by_lineage),"; at FDR <",fdr_cut,":",
      sum(by_lineage$FDR<fdr_cut,na.rm=TRUE),"\n")
}
