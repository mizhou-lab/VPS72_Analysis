# A1 gain frequency, A2 CN-expression slopes, A3 gain effects, A6 dosage sensitivity.
# Deliberately excludes the INVALID in-sample OLS-residual and all-sample-z-score tests
# from the original 01_1q21_gain_dosage_distinctness.R (formerly A4 and A5).
source("00_config.R")
check_packages("data.table")
library(data.table)
set.seed(42)
check_required(unlist(FILES_1Q21[c("gistic","log2cn","z_dip","z_all")]))

raw_list <- list(
  read_cbio_long(FILES_1Q21$gistic, "gistic"),
  read_cbio_long(FILES_1Q21$log2cn, "log2cn"),
  read_cbio_long(FILES_1Q21$z_dip, "z_dip"),
  read_cbio_long(FILES_1Q21$z_all, "z_all")
)
keys <- c("study","sample","gene")
dat <- Reduce(function(a,b) merge(a,b, by=keys, all=FALSE), raw_list)
if (file.exists(FILES_1Q21$z_norm)) {
  zn <- read_cbio_long(FILES_1Q21$z_norm, "z_norm")
  dat <- merge(dat, zn, by=keys, all.x=TRUE)
} else {
  message("Normal-reference z-score input absent. Script 02 will run held-out models but skip A5b.")
  dat[, z_norm := NA_real_]
}
if (!nrow(dat)) stop("No common samples across 1q21 molecular inputs")
dat[, cancer := toupper(sub("_tcga.*$", "", study))]
dat[, patient := patient_id(sample)]
# One tumour sample per patient/study, selected reproducibly by sample barcode.
ids <- unique(dat[, .(study, sample, patient)])[order(study, patient, sample)]
ids <- ids[!duplicated(ids[, .(study, patient)]), .(study, sample)]
dat <- dat[ids, on=.(study, sample), nomatch=0]
dat[, cn_group := fcase(is.na(gistic), NA_character_, gistic <= -1, "loss",
                        gistic == 0, "neutral", gistic == 1, "gain", gistic == 2, "amp",
                        default = NA_character_)]
dat[, gain := as.integer(gistic >= 1)]
dat[, amp := as.integer(gistic == 2)]
setorder(dat, cancer, study, sample, gene)
saveRDS(dat, out_file("1q21_merged_long.rds"))
fwrite(dat[, .(study,sample,cancer,gene,gistic,log2cn,z_dip,z_all,z_norm)],
       out_file("00_primary_tumor_molecular_long.tsv.gz"), sep="\t")
qc <- data.table(metric=c("primary_tumor_samples","cancer_types","genes","normal_reference_available"),
                 value=c(uniqueN(dat$sample),uniqueN(dat$cancer),uniqueN(dat$gene),
                         file.exists(FILES_1Q21$z_norm)))
fwrite(qc, out_file("00_QC.tsv"),sep="\t")

# A1: GISTIC gain >= +1, high-level amplification = +2.
a1_pan <- dat[!is.na(gistic), .(n=.N,gain_pct=100*mean(gain),amp_pct=100*mean(amp)),by=gene][order(-gain_pct)]
a1_ct  <- dat[!is.na(gistic), .(n=.N,gain_pct=100*mean(gain),amp_pct=100*mean(amp)),by=.(cancer,gene)]
fwrite(a1_pan,out_file("A1_gain_freq_pancancer.tsv"),sep="\t")
fwrite(a1_ct,out_file("A1_gain_freq_by_cancer.tsv"),sep="\t")

# A2: within-cancer linear slope (z relative to diploids ~ continuous log2 CN).
a2 <- dat[is.finite(log2cn) & is.finite(z_dip), {
  if (.N >= 30L && sd(log2cn) > 0 && sd(z_dip) > 0) {
    fit <- lm(z_dip ~ log2cn)
    z <- summary(fit)$coefficients["log2cn",]
    .(n=.N,slope=unname(z["Estimate"]),se=unname(z["Std. Error"]),
      p=unname(z["Pr(>|t|)"]),r2=summary(fit)$r.squared)
  } else .(n=.N,slope=NA_real_,se=NA_real_,p=NA_real_,r2=NA_real_)
},by=.(cancer,gene)]
a2[, fdr := p.adjust(p,"BH"),by=cancer]
a2_pan <- a2[is.finite(slope),.(n_cancers=.N,median_slope=median(slope),
                                    mean_slope=mean(slope),median_r2=median(r2),
                                    n_sig=sum(fdr<0.05,na.rm=TRUE)),by=gene][order(-median_slope)]
a2_pan[,rank:=seq_len(.N)]
fwrite(a2,out_file("A2_dosage_slope_by_cancer.tsv"),sep="\t")
fwrite(a2_pan,out_file("A2_dosage_slope_pancancer.tsv"),sep="\t")
a2_wide <- dcast(a2, cancer~gene, value.var="slope")
pairwise <- rbindlist(lapply(NEIGHBOURS,function(g) {
  d <- a2_wide[[GOI]]-a2_wide[[g]]; d <- d[is.finite(d)]
  if (length(d)<5L) return(NULL)
  data.table(neighbour=g,n_cancers=length(d),median_diff=median(d),
             p=tryCatch(suppressWarnings(wilcox.test(d,mu=0,exact=FALSE)$p.value),
                        error=function(e) NA_real_))
}),fill=TRUE)
if (nrow(pairwise)) pairwise[,fdr:=p.adjust(p,"BH")]
fwrite(pairwise,out_file("A2_pairwise_vs_neighbours.tsv"),sep="\t")

# A3: gain/amplification vs copy-neutral within cancer; Cliff's delta and Wilcoxon.
a3 <- dat[!is.na(cn_group) & is.finite(z_dip), {
  hi<-z_dip[cn_group %in% c("gain","amp")]; lo<-z_dip[cn_group=="neutral"]
  if (length(hi)>=MIN_GROUP && length(lo)>=MIN_GROUP) {
    .(n_hi=length(hi),n_lo=length(lo),delta_median=median(hi)-median(lo),
      cliffs_d=cliffs_delta(hi,lo),
      p=suppressWarnings(wilcox.test(hi,lo,exact=FALSE)$p.value))
  } else .(n_hi=length(hi),n_lo=length(lo),delta_median=NA_real_,cliffs_d=NA_real_,p=NA_real_)
},by=.(cancer,gene)]
a3[,fdr:=p.adjust(p,"BH"),by=cancer]
a3_pan<-a3[is.finite(cliffs_d),.(n_cancers=.N,median_cliffs=median(cliffs_d),
                                  median_delta=median(delta_median),
                                  n_sig_positive=sum(fdr<0.05 & cliffs_d>0,na.rm=TRUE)),
           by=gene][order(-median_cliffs)]
a3_pan[,rank:=seq_len(.N)]
fwrite(a3,out_file("A3_gain_vs_diploid_by_cancer.tsv"),sep="\t")
fwrite(a3_pan,out_file("A3_gain_vs_diploid_pancancer.tsv"),sep="\t")

# A6: descriptive sensitivity index (NOT a significance test).
a6<-dat[cn_group %in% c("gain","amp") & is.finite(z_dip),
        .(n_gained=.N,pct_overexpr=100*mean(z_dip>1)),by=gene]
a6<-merge(a6,a2_pan[,.(gene,median_slope)],by="gene")
a6<-merge(a6,a3_pan[,.(gene,median_cliffs)],by="gene")
a6[,`:=`(r_overexpr=frank(-pct_overexpr,ties.method="average"),
            r_slope=frank(-median_slope,ties.method="average"),
            r_cliffs=frank(-median_cliffs,ties.method="average"))]
a6[,mean_rank:=(r_overexpr+r_slope+r_cliffs)/3]
setorder(a6,mean_rank,gene)
fwrite(a6,out_file("A6_dosage_sensitivity_index.tsv"),sep="\t")

# Correct cancer-type bootstrap: preserve resampled cancer multiplicity.
cts<-sort(unique(a2$cancer))
boot <- rbindlist(lapply(seq_len(1000L),function(i) {
  resampled<-sample(cts,length(cts),replace=TRUE)
  d<-rbindlist(lapply(resampled,function(ct) a2[cancer==ct & is.finite(slope),.(gene,slope)]))
  summ<-d[,.(median_slope=median(slope)),by=gene]
  summ[,rank:=frank(-median_slope,ties.method="min")]
  data.table(bootstrap=i,rank_VPS72=summ[gene==GOI,rank][1])
}))
fwrite(boot,out_file("A6_bootstrap_VPS72_slope_rank.tsv"),sep="\t")
writeLines(capture.output(sessionInfo()),out_file("sessionInfo.txt"))
message("01 complete: primary-tumor gain, dosage slopes, gain effects and descriptive ranking; outputs: ",OUT_DIR)
