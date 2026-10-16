# Composite disruption score ranking the complex genes.
# The three component means are z-scored across the panel and summed, so the
# score is a position within this gene set, not an absolute quantity.
source("00_config.R")

read_long <- function(p,value_name) {
  d <- fread(require_input(p))
  setnames(d,1,"gene")
  melt(d,id.vars="gene",variable.name="cancer",value.name=value_name)
}

merged <- Reduce(function(a,b) merge(a,b,by=c("gene","cancer"),all=TRUE),
                 list(read_long(inputs$expr_mat,"expr"),
                      read_long(inputs$cnv_mat,"cnv"),
                      read_long(inputs$amp_mat,"amp")))
if (!nrow(merged)) stop("No rows after merging the three component matrices")

# n_cancers travels with the score because a gene measured in fewer cancer
# types has a less stable mean, which the score itself does not express.
score <- merged[,.(expr_mean=mean(expr,na.rm=TRUE),
                   cnv_mean=mean(cnv,na.rm=TRUE),
                   amp_mean=mean(amp,na.rm=TRUE),
                   expr_frac_ge_z1=mean(expr>=1,na.rm=TRUE),
                   cnv_frac_ge_0.1=mean(cnv>=0.10,na.rm=TRUE),
                   amp_frac_ge_5pct=mean(amp>=5,na.rm=TRUE),
                   n_cancers=sum(!is.na(expr))),by=gene]
if (nrow(score)<3L) stop("Too few genes to z-score the components")

score[,`:=`(expr_z=as.numeric(scale(expr_mean)),
            cnv_z=as.numeric(scale(cnv_mean)),
            amp_z=as.numeric(scale(amp_mean)))]
score[,disruption:=expr_z+cnv_z+amp_z]
score[,rank:=rank(-disruption,ties.method="min")]
# Fraction of panel genes scoring at or above each gene: a within-panel
# percentile over 15 genes, not a test against a null model.
score[,empirical_p:=vapply(disruption,function(x) mean(disruption>=x,na.rm=TRUE),numeric(1))]
setorder(score,rank)

write_tsv(score,"complex_disruption_score.tsv")
print(score[,.(gene,disruption=round(disruption,3),rank,empirical_p)])
