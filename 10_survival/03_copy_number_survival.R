# cBioPortal TCGA copy-number vs OS: continuous and categorical gain models.
# log2CN > 0.2 primary; >0.585 amplification; upper 30% exploratory sensitivity.
# Does NOT mix these log2 thresholds with GISTIC categorical calls.
source("00_config.R"); need(c("data.table","survival"))
cn <- data.table::fread(required(CBIO_CNV_FILE)); cl <- data.table::fread(required(CBIO_CLIN_FILE))
for(k in c("Sample ID","Study ID")) if(!k %in% names(cn)) stop("CNV missing ",k)
value_col <- pick_col(names(cn),c("VPS72: Log2 copy-number values","VPS72"))
ct_id <- pick_col(names(cl),c("Sample ID","SAMPLE_ID"))
ct_study <- pick_col(names(cl),c("Study ID","STUDY_ID"),required=FALSE)
time_col <- pick_col(names(cl),c("Overall Survival (Months)","OS Months","OS_MONTHS"))
status_col <- pick_col(names(cl),c("Overall Survival Status","OS Status","OS_STATUS"))
a <- data.table::data.table(sample=clean_id(cn[["Sample ID"]]),study=as.character(cn[["Study ID"]]),CNV=as.numeric(cn[[value_col]]))
a[,`:=`(patient=tcga12(sample),sample15=tcga15(sample),cancer=toupper(sub("_tcga.*$","",study)))]
a <- a[!is.na(sample15)&substr(sample15,14,15)=="01"&is.finite(CNV)]
a <- a[order(patient,sample)][!duplicated(patient)]
b <- data.table::data.table(patient=tcga12(cl[[ct_id]]),time_months=as.numeric(cl[[time_col]]),event=parse_event(cl[[status_col]]))
b <- b[!is.na(patient)&is.finite(time_months)&time_months>0&!is.na(event)]
if(b[,.(n=uniqueN(round(time_months,5)),e=uniqueN(event)),by=patient][n>1|e>1,.N]) stop("Conflicting clinical OS entries")
b <- unique(b,by="patient"); d <- merge(a,b,by="patient"); d <- d[!cancer %in% LIQUID]
if(!nrow(d)) stop("No matching primary-tumor TCGA CNV and survival rows")
d[,gained:=as.integer(CNV>0.2)]; d[,amplified:=as.integer(CNV>0.585)]
d[,top30:=as.integer(CNV>=as.numeric(quantile(CNV,.70))),by=cancer]
write_table(d,"03_CNV_survival_samples.tsv")
out <- list()
add <- function(result,ca,subset) if(!is.null(result)) out[[length(out)+1L]] <<- cbind(cancer=ca,cohort=subset,result)
for(ct in sort(unique(d$cancer))) for(v in c("CNV","gained","amplified","top30"))
  add(run_cox(d[cancer==ct],v,label=v),ct,"per_cancer")
for(nm in c("all_solid","significantly_upregulated_fixed")) {
  sub <- if(nm=="all_solid") d else d[cancer %in% SIG_UP]; if(!nrow(sub)) next
  for(v in c("CNV","gained","amplified","top30"))
    add(run_cox(sub,v,strata="cancer",label=paste0("cancer_stratified_",v)),"PAN",nm)
}
if(!length(out)) stop("No evaluable CNV-survival models")
r <- data.table::rbindlist(out,fill=TRUE); r[,FDR:=p.adjust(p,"BH"),by=.(cohort,model)]
write_table(r,"03_CNV_survival_Cox_models.tsv")
