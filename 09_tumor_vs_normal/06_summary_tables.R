# Source data behind the manuscript table: the median-split and continuous
# analyses side by side, per cancer type and pooled across cancer types.
source("00_config.R")

split_res <- fread(require_input(inputs$highlow))
cors <- fread(require_input(inputs$spearman))
meta_state <- fread(require_input(inputs$meta_state))
ints <- fread(require_input(inputs$interact))
diff_meta <- fread(require_input(file.path(OUT_DIR,"tumor_vs_normal_meta_regression.tsv")))
int_meta <- fread(require_input(file.path(OUT_DIR,"interaction_meta.tsv")))

pw <- "proliferation_score"
if (!pw %chin% cors$pathway) stop("Composite proliferation absent from the correlation table")

# Per cancer type: the split result, both correlations, and the interaction.
wide_split <- dcast(split_res[pathway==pw],cancer~state,
                    value.var=c("delta_high_low","p","FDR","positive"))
wide_cors <- dcast(cors[pathway==pw],cancer~state,
                   value.var=c("n","rho","ci_low","ci_high","p"))
per_cancer <- Reduce(function(a,b) merge(a,b,by="cancer",all=TRUE),
                     list(wide_split,wide_cors,
                          ints[pathway==pw,.(cancer,interaction_beta,interaction_CI_low=CI_low,
                                             interaction_CI_high=CI_high,
                                             interaction_p=p,interaction_FDR=FDR)]))
per_cancer[,ord:=match(cancer,eligible)]
setorder(per_cancer,ord); per_cancer[,ord:=NULL]
write_tsv(per_cancer,"table_source_per_cancer.tsv")

# Pooled across cancer types, for every programme scored. The two analyses
# answer the same question differently: the split counts in how many cancer
# types VPS72-high tumours carry a higher programme score, the continuous arm
# pools the correlation itself and tests the tumour-normal difference formally.
recurrence <- dcast(split_res[,.(tested=.N,positive=sum(positive,na.rm=TRUE),
                                 median_delta=median(delta_high_low,na.rm=TRUE)),
                              by=.(pathway,Program,state)],
                    pathway+Program~state,value.var=c("tested","positive","median_delta"))
if (!all(c("positive_Tumor","positive_Normal") %chin% names(recurrence)))
  stop("Median-split results are missing one tissue state")
recurrence[,`:=`(pct_Tumor=100*positive_Tumor/tested_Tumor,
                 pct_Normal=100*positive_Normal/tested_Normal)]
recurrence[,recurrence_difference_pp:=pct_Tumor-pct_Normal]

pooled <- merge(
  dcast(meta_state,pathway+Program~state,
        value.var=c("n_cancers","pooled_rho","CI_low","CI_high","p","FDR","I2")),
  diff_meta[,.(pathway,fisher_z_difference,diff_CI_low=CI_low,diff_CI_high=CI_high,
               diff_p=p,diff_FDR=FDR)],by="pathway",all.x=TRUE)
pooled <- merge(pooled,
                int_meta[,.(pathway,pooled_interaction_beta,int_CI_low=CI_low,
                            int_CI_high=CI_high,int_p=p,int_FDR=FDR)],
                by="pathway",all.x=TRUE)
pooled <- merge(recurrence,pooled,by=c("pathway","Program"),all=TRUE)

# Composite first, then the individual programmes in the order they are listed
# in 00_config.R.
pooled[,ord:=match(pathway,names(program_labels))]
setorder(pooled,ord); pooled[,ord:=NULL]
write_tsv(pooled,"table_source_pooled.tsv")

# The formatted columns the manuscript table reports, kept alongside the raw
# estimates rather than in place of them.
formatted <- pooled[,.(Program,
  `VPS72-high tumours`=sprintf("%d/%d (%.0f%%)",positive_Tumor,tested_Tumor,pct_Tumor),
  `Tumour median delta`=round(median_delta_Tumor,2),
  `VPS72-high normal`=sprintf("%d/%d (%.0f%%)",positive_Normal,tested_Normal,pct_Normal),
  `Normal median delta`=round(median_delta_Normal,2),
  `Recurrence difference (pp)`=round(recurrence_difference_pp),
  `Pooled rho, tumour`=sprintf("%.2f (%.2f to %.2f)",pooled_rho_Tumor,CI_low_Tumor,CI_high_Tumor),
  `Pooled rho, normal`=sprintf("%.2f (%.2f to %.2f)",pooled_rho_Normal,CI_low_Normal,CI_high_Normal),
  `Tumour vs normal, Fisher z`=sprintf("%.3f (p = %.3g)",fisher_z_difference,diff_p),
  `Pooled interaction beta`=sprintf("%.3f (p = %.3g)",pooled_interaction_beta,int_p))]
write_tsv(formatted,"table_formatted.tsv")

cat("Programmes summarised:",nrow(pooled),"; cancer types in the per-cancer table:",
    nrow(per_cancer),"\n")
print(formatted[,.(Program,`VPS72-high tumours`,`VPS72-high normal`,`Recurrence difference (pp)`)])
