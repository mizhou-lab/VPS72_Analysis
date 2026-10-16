# Hallmark GSEA for each perturbation, and the pairwise NES concordance that
# asks whether VPS72 moves the same programmes as SRCAP and as H2A.Z.1.
source("00_config.R")
suppressPackageStartupMessages({library(fgsea);library(msigdbr)})

deg <- fread(require_input(inputs$deg))
hall <- as.data.table(msigdbr(species="Homo sapiens",category="H"))
pathways <- split(toupper(hall$gene_symbol),hall$gs_name)
if (!length(pathways)) stop("No Hallmark sets returned by msigdbr")

# Genes are ranked by log2 fold change alone, not by a signed significance
# statistic. A low-count gene with a large fold change therefore ranks as highly
# as a well-measured one; this is the ranking the published analysis used.
run_one <- function(ct) {
  d <- deg[contrast==ct & !is.na(gene) & is.finite(log2fc)]
  d <- d[,.SD[1L],by=gene]
  if (nrow(d)<500L) stop("Only ",nrow(d)," ranked genes for ",ct)
  r <- setNames(d$log2fc,d$gene)
  r <- sort(r,decreasing=TRUE)
  res <- as.data.table(fgseaMultilevel(pathways=pathways,stats=r))
  res[,`:=`(contrast=ct,leadingEdge=vapply(leadingEdge,paste,character(1),collapse="/"))]
  res[]
}

fg <- rbindlist(lapply(c("VPS72","SRCAP","H2AZ1"),run_one),fill=TRUE)
setorder(fg,contrast,padj)
write_tsv(fg,"fgsea_hallmark.tsv")
print(fg[,.(n_pathways=.N,n_sig=sum(padj<padj_cut,na.rm=TRUE)),by=contrast])

# Concordance: a pathway absent from one contrast is carried as NES 0 rather
# than dropped, so the comparison covers the same pathway set on both axes.
concordance <- function(a,b) {
  A <- fg[contrast==a,.(pathway,NES_a=NES,padj_a=padj)]
  B <- fg[contrast==b,.(pathway,NES_b=NES,padj_b=padj)]
  m <- merge(A,B,by="pathway",all=TRUE)
  m[is.na(NES_a),NES_a:=0]; m[is.na(NES_b),NES_b:=0]
  m[,`:=`(sig_a=!is.na(padj_a) & padj_a<padj_cut,
          sig_b=!is.na(padj_b) & padj_b<padj_cut)]
  m[,sig_any:=sig_a|sig_b]
  m[,agreement:=fcase(
    sig_any & NES_a>0 & NES_b>0,"up_concordant",
    sig_any & NES_a<0 & NES_b<0,"down_concordant",
    sig_any & sign(NES_a)!=sign(NES_b),"discordant",
    default="not_significant")]
  m[,`:=`(contrast_a=a,contrast_b=b)]
  setcolorder(m,c("contrast_a","contrast_b","pathway"))
  m[]
}

conc <- rbindlist(list(concordance("VPS72","SRCAP"),concordance("VPS72","H2AZ1")))
write_tsv(conc,"fgsea_nes_concordance.tsv")

# One correlation per pair over the pathways significant in either contrast:
# the headline number behind the concordance panels.
summ <- conc[sig_any==TRUE,.(n_pathways=.N,
                             spearman_rho=cor(NES_a,NES_b,method="spearman"),
                             pearson_r=cor(NES_a,NES_b),
                             up_concordant=sum(agreement=="up_concordant"),
                             down_concordant=sum(agreement=="down_concordant"),
                             discordant=sum(agreement=="discordant")),
             by=.(contrast_a,contrast_b)]
write_tsv(summ,"fgsea_concordance_summary.tsv")
print(summ)
