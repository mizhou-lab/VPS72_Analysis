# VPS72 expression against genomic-instability measures: CIN, aneuploidy,
# ploidy, LOH fraction, MSI score, and whole-genome doubling.
source("00_config.R")

model <- read_model()
expr <- read_gene_export(inputs$expr_gene,"expr_log1p_tpm")

# The global signature table is keyed by profile, not model, so it joins
# through OmicsProfiles.
glob <- clean_names(fread(require_input(inputs$global_sig)))
setnames(glob,1,"ProfileID")
prof <- clean_names(fread(require_input(inputs$profiles)))
setnames(prof,c(find_col(prof,c("ProfileID","profile_id")),
                find_col(prof,c("ModelID","DepMap_ID"))),c("ProfileID","ModelID"))

dat <- merge(glob,unique(prof[,.(ProfileID,ModelID)]),by="ProfileID")
dat <- merge(dat,expr,by="ModelID")
dat <- merge(dat,model[,.(ModelID,lineage,is_cancer)],by="ModelID")
dat <- dat[is.finite(expr_log1p_tpm)]
if (!nrow(dat)) stop("No cell line carried both a global signature profile and VPS72 expression")
cat("lines with instability metrics:",nrow(dat),"\n")

present <- intersect(instability_metrics,names(dat))
if (!length(present)) stop("None of the expected instability metrics is in ",inputs$global_sig,
                           "\nAvailable: ",paste(names(dat),collapse=", "))
absent <- setdiff(instability_metrics,present)
if (length(absent)) cat("not in this release:",paste(absent,collapse=", "),"\n")

pan <- rbindlist(lapply(present,function(m) {
  s <- spearman(suppressWarnings(as.numeric(dat[[m]])),dat$expr_log1p_tpm)
  if (is.null(s)) NULL else { s[,metric:=m]; s[] }
}))
if (!nrow(pan)) stop("No instability metric had enough paired values")
pan[,FDR:=p.adjust(p,method="BH")]
setorder(pan,p)
write_tsv(pan,"instability_pancancer.tsv")
print(pan)

# Per lineage, with BH applied within each metric across lineages: each metric
# is its own family of tests.
by_lineage <- rbindlist(lapply(present,function(m) {
  v <- suppressWarnings(as.numeric(dat[[m]]))
  d <- data.table(lineage=dat$lineage,metric_value=v,expr=dat$expr_log1p_tpm)
  d <- d[is.finite(metric_value) & is.finite(expr) & !is.na(lineage) & nzchar(lineage)]
  rbindlist(lapply(d[,.N,by=lineage][N>=min_lineage_n]$lineage,function(l) {
    s <- spearman(d[lineage==l]$metric_value,d[lineage==l]$expr,min_n=min_lineage_n)
    if (is.null(s)) NULL else { s[,`:=`(metric=m,lineage=l)]; s[] }
  }))
}))
if (nrow(by_lineage)) {
  by_lineage[,FDR:=p.adjust(p,method="BH"),by=metric]
  setorder(by_lineage,metric,p)
  write_tsv(by_lineage,"instability_by_lineage.tsv")
  cat("Metric-lineage pairs tested:",nrow(by_lineage),"; at FDR <",fdr_cut,":",
      sum(by_lineage$FDR<fdr_cut,na.rm=TRUE),"\n")
}

# Whole-genome doubling is a flag, not a score, so it is a group comparison.
wgd_col <- names(dat)[grepl("^WGD$|WholeGenomeDoubling",names(dat),ignore.case=TRUE)][1L]
if (is.na(wgd_col)) {
  cat("no WGD column in this release; comparison skipped\n")
} else {
  d <- dat[,.(wgd=suppressWarnings(as.integer(get(wgd_col))),expr=expr_log1p_tpm)]
  d <- d[wgd %in% c(0L,1L) & is.finite(expr)]
  if (uniqueN(d$wgd)<2L || nrow(d)<20L) {
    cat("too few lines in one WGD group; comparison skipped\n")
  } else {
    wt <- wilcox.test(expr~wgd,data=d)
    out <- data.table(comparison="WGD_1_vs_0",
                      n_wgd0=sum(d$wgd==0L),n_wgd1=sum(d$wgd==1L),
                      median_wgd0=median(d[wgd==0L]$expr),
                      median_wgd1=median(d[wgd==1L]$expr),
                      W=as.numeric(wt$statistic),p=wt$p.value)
    write_tsv(out,"instability_wgd.tsv")
    print(out)
  }
}
