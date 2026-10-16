# Optional: PCAWG / ICGC pan-cancer VPS72 survival (cBioPortal clinical + FPKM-UQ).
# Unlike TCGA Xena, this is an independent cohort; never silently pool datasets.
source("00_config.R"); need(c("data.table","survival"))
cl <- data.table::fread(required(PCAWG_CLIN_FILE)); ex <- data.table::fread(required(PCAWG_EXPR_FILE))
patient_cl <- pick_col(names(cl),c("Patient ID","PATIENT_ID","patient_id"))
patient_ex <- pick_col(names(ex),c("Patient ID","PATIENT_ID","patient_id"))
time_col <- pick_col(names(cl),c("Overall Survival (Months)","OS Months","overall_survival_months"))
event_col <- pick_col(names(cl),c("Overall Survival Status","OS Status","overall_survival_status"))
expr_cols <- grep("VPS72",names(ex),ignore.case=TRUE,value=TRUE)
expr_cols <- setdiff(expr_cols,c("VPS72_group","VPS72_status"))
if(length(expr_cols)!=1) stop("Identify exactly one VPS72 FPKM-UQ column; observed: ",paste(expr_cols,collapse=", "))
a <- data.table::data.table(patient=trimws(as.character(ex[[patient_ex]])),FPKM_UQ=as.numeric(ex[[expr_cols]]))
a <- a[patient!=""&!is.na(patient)&is.finite(FPKM_UQ)&FPKM_UQ>=0]
a <- a[,.(FPKM_UQ=median(FPKM_UQ)),by=patient]
b <- data.table::data.table(patient=trimws(as.character(cl[[patient_cl]])),time_months=as.numeric(cl[[time_col]]),event=parse_event(cl[[event_col]]))
b <- b[patient!=""&!is.na(patient)&is.finite(time_months)&time_months>0&!is.na(event)]
if(b[,.(n_time=uniqueN(round(time_months,5)),n_event=uniqueN(event)),by=patient][n_time>1|n_event>1,.N]) stop("PCAWG duplicate patients have conflicting survival")
b <- unique(b,by="patient")
d <- merge(a,b,by="patient"); if(nrow(d)<30) stop("Insufficient matched PCAWG samples")
d[,expression:=log2(FPKM_UQ+1)]; d[,vps72_z:=as.numeric(scale(expression))]
d[,group:=factor(ifelse(expression>=median(expression),"High","Low"),levels=c("Low","High"))]
d[,group_high:=as.integer(group=="High")]
write_table(d,"06_PCAWG_survival_samples.tsv")
r1 <- run_cox(d,"vps72_z",label="continuous_per_SD")
r2 <- run_cox(d,"group_high",label="median_split_high_vs_low")
res <- data.table::rbindlist(Filter(Negate(is.null),list(r1,r2)),fill=TRUE)
if(nrow(res)) write_table(res,"06_PCAWG_survival_Cox_models.tsv")
# Original PCAWG script also calculated unadjusted RMST. Preserve as an optional,
# explicitly dependency-checked analysis instead of auto-installing packages.
if(requireNamespace("survRM2",quietly=TRUE) && data.table::uniqueN(d$group_high)==2) {
  max_by_group <- d[,.(last=max(time_months)),by=group_high]
  tau <- floor(min(max_by_group$last))-1
  if(is.finite(tau) && tau>0 && all(d[,.(N=.N),by=group_high]$N>=10)) {
    m <- tryCatch(survRM2::rmst2(d$time_months,d$event,d$group_high,tau=tau),error=function(e) NULL)
    if(!is.null(m)) {
      z <- data.table::as.data.table(m$unadjusted.result,keep.rownames="contrast")
      z[,tau_months:=tau]
      write_table(z,"06_PCAWG_unadjusted_RMST.tsv")
    }
  }
} else message("PCAWG RMST skipped: install survRM2 to reproduce the optional RMST analysis")
