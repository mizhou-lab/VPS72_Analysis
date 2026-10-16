# Primary-tumor TCGA Xena: Kaplan-Meier source data, pooled/per-cancer Cox, quartile comparisons.
source("00_config.R"); need(c("data.table","survival"))
x <- data.table::fread(required(XENA_EXPR_FILE)); c <- data.table::fread(required(XENA_CLIN_FILE))
for(k in c("sample","sample_type","VPS72","cancer type abbreviation")) if(!k %in% names(x)) stop("Expression missing: ",k)
for(k in c("sample","OS","OS.time")) if(!k %in% names(c)) stop("Clinical missing: ",k)
x <- x[sample_type=="Primary Tumor",.(sample=as.character(sample),cancer=as.character(`cancer type abbreviation`),expression=as.numeric(VPS72))]
x[,`:=`(sample15=tcga15(sample),patient=tcga12(sample))]
x <- x[!is.na(patient) & !is.na(sample15) & substr(sample15,14,15)=="01" & is.finite(expression) & expression>=0]
x <- x[order(patient,sample15,sample)][!duplicated(patient)]
# Xena clinical 'sample' may be an aliquot (or have an extra trailing character), so link
# on unique patient ID only after checking for incompatible survival records.
c[,patient:=tcga12(sample)]; c <- c[!is.na(patient)]
c[,`:=`(time_months=as.numeric(OS.time)/30.4375,event=parse_event(OS))] # Xena OS.time is in DAYS.
c <- c[is.finite(time_months) & time_months>0 & !is.na(event)]
conflict <- c[,.(times=data.table::uniqueN(round(time_months,5)),events=data.table::uniqueN(event)),by=patient][times>1|events>1]
if(nrow(conflict)) stop("Clinical records disagree within patient; resolve: ",paste(head(conflict$patient),collapse=", "))
c <- unique(c[,.(patient,time_months,event)],by="patient")
d <- merge(x,c,by="patient",all=FALSE); d <- d[!cancer %in% LIQUID & !is.na(cancer) & cancer!=""]
if(!nrow(d)) stop("No matching primary-tumor clinical/expression records")
d[,expr_log2:=log2(expression+1)]; d[,expr_z:=if(.N>2&&sd(expr_log2)>0) as.numeric(scale(expr_log2)) else NA_real_,by=cancer]
d[,group_cancer:=ifelse(expr_log2>=median(expr_log2),"High","Low"),by=cancer]
write_table(d,"01_TCGA_survival_eligible_samples.tsv")
# The pooled unstratified KM split is descriptive, as in the uploaded figure scripts.
# Continuous Cox with cancer strata estimates association accounting for baseline risks.
res <- list()
add <- function(a,ct,co) if(!is.null(a)) res[[length(res)+1L]] <<- cbind(cancer=ct,cohort=co,a)
for(ct in sort(unique(d$cancer))) {
  s <- d[cancer==ct]; s[,med_grp:=as.integer(expr_log2>=median(expr_log2))]
  add(run_cox(s,"expr_log2",label="continuous_log2_expression"),ct,"per_cancer")
  add(run_cox(s,"med_grp",label="within_cancer_median_split"),ct,"per_cancer")
}
for(nm in c("all_solid","upregulated_fixed","significantly_upregulated_fixed")) {
  s <- switch(nm,all_solid=d,upregulated_fixed=d[cancer %in% UP],significantly_upregulated_fixed=d[cancer %in% SIG_UP]); if(!nrow(s)) next
  s[,pooled_med_grp:=as.integer(expr_log2>=median(expr_log2))]
  s[,within_cancer_grp:=as.integer(group_cancer=="High")]
  add(run_cox(s,"pooled_med_grp",label="descriptive_pooled_median"),"PAN",nm)
  add(run_cox(s,"expr_z",strata="cancer",label="cancer_stratified_per_SD"),"PAN",nm)
  add(run_cox(s,"within_cancer_grp",strata="cancer",label="cancer_stratified_median"),"PAN",nm)
  # KM is reproducible directly from these per-sample group/event/time data.
  s[,`:=`(group_pooled=ifelse(pooled_med_grp==1,"High","Low"),subset=nm)]
  write_table(s[,.(patient,sample15,cancer,subset,time_months,event,expr_log2,group_pooled,group_cancer)],paste0("01_KM_data_",nm,".tsv"))
  # Quartile exploratory contrasts are calculated within each defined subset.
  q <- stats::quantile(s$expr_log2,c(.25,.75),na.rm=TRUE)
  s[,`:=`(q4=as.integer(expr_log2>=q[2]),extreme=ifelse(expr_log2<=q[1],0,ifelse(expr_log2>=q[2],1,NA_real_)))]
  add(run_cox(s,"q4",label="top_quartile_vs_rest"),"PAN",nm)
  add(run_cox(s[!is.na(extreme)],"extreme",label="top_vs_bottom_quartile"),"PAN",nm)
}
r <- data.table::rbindlist(res,fill=TRUE); if(!nrow(r)) stop("No evaluable Cox models")
r[,FDR:=p.adjust(p,"BH"),by=.(cohort,model)]
write_table(r,"01_TCGA_survival_Cox_models.tsv")
