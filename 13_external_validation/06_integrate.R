# Bring the two assays together: for each signature gene, its expression
# response to YL1 knockdown next to its promoter H4ac change, and the
# correlation between them.
source("00_config.R")

deg <- fread(require_input(inputs$deg))
h4 <- fread(require_input(file.path(OUT_DIR,"promoter_h4ac_ranking.tsv")))

wide <- dcast(deg[contrast %chin% c("VPS72","SRCAP","H2AZ1")],
              gene~contrast,value.var=c("log2fc","padj","direction"))
m <- merge(h4[,.(gene,signatures,control_mean,kd_mean,h4ac_delta=delta,
                 h4ac_log2_ratio=log2_ratio,h4ac_rank=rank)],
           wide,by="gene")
if (!nrow(m)) stop("No gene carried both an H4ac measurement and a DEG result")
setorder(m,h4ac_rank)
write_tsv(m,"expression_vs_h4ac.tsv")

# If YL1 loss reduces promoter acetylation and that drives transcription down,
# the two deltas should move together. Reported over all measured genes and
# again over those called significant in the YL1 contrast.
cor_block <- function(d,label) {
  d <- d[is.finite(h4ac_delta) & is.finite(log2fc_VPS72)]
  if (nrow(d)<10L) return(NULL)
  ct <- suppressWarnings(cor.test(d$h4ac_delta,d$log2fc_VPS72,method="spearman",exact=FALSE))
  data.table(gene_set=label,n=nrow(d),spearman_rho=as.numeric(ct$estimate),p=ct$p.value)
}
summ <- rbindlist(list(cor_block(m,"all_signature_genes"),
                       cor_block(m[direction_VPS72!="NS"],"YL1_responsive_only")))
if (!nrow(summ)) stop("Too few genes with both measurements to correlate")
write_tsv(summ,"expression_h4ac_correlation.tsv")

# Per signature, how the YL1-responsive genes behave on both axes.
by_sig <- m[direction_VPS72!="NS",.(n=.N,
             median_lfc=median(log2fc_VPS72,na.rm=TRUE),
             median_h4ac_delta=median(h4ac_delta,na.rm=TRUE),
             n_lost_h4ac=sum(h4ac_delta<0,na.rm=TRUE)),
            by=.(signature=signatures)]
setorder(by_sig,-n)
write_tsv(by_sig,"signature_level_summary.tsv")

print(summ)
cat("Genes with both measurements:",nrow(m),
    "; YL1-responsive among them:",sum(m$direction_VPS72!="NS"),"\n")
