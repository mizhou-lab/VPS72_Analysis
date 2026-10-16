# Per-cancer apparent Cox C-index and 1/3/5-year IPCW time-dependent AUC; TSV only.
# Requires 01_tcga_overall_survival.R output; no arbitrary scatter-plot exclusions.
source("00_config.R"); need(c("data.table","survival","timeROC"))
d <- data.table::fread(output("01_TCGA_survival_eligible_samples.tsv"))
if(!nrow(d)) stop("Run 01_tcga_overall_survival.R first")
res <- list(); auc <- list()
for(ct in sort(unique(d$cancer))) {
  s <- d[cancer==ct & is.finite(time_months) & time_months>0 & !is.na(event) & is.finite(expr_log2)]
  if(nrow(s)<30 || sum(s$event)<10 || length(unique(s$expr_log2))<2) next
  fit <- tryCatch(survival::coxph(survival::Surv(time_months,event)~expr_log2,data=s),error=function(e) NULL)
  if(is.null(fit)) next
  z <- summary(fit); concord <- z$concordance
  ci <- unname(concord["C"]); se <- unname(concord["se(C)"])
  res[[length(res)+1L]] <- data.frame(cancer=ct,n=nrow(s),events=sum(s$event),
    c_index=ci,ci_low=max(0,ci-1.96*se),ci_high=min(1,ci+1.96*se),
    cox_beta=unname(stats::coef(fit)[1]),p=unname(z$coefficients[1,"Pr(>|z|)"]))
  # AUC uses fitted Cox log-risk rather than raw expression so that the marker
  # always points toward *higher hazard* (original script did not adjust direction).
  risk <- as.numeric(stats::predict(fit,type="lp"))
  for(t in c(12,36,60)) {
    # Need events before t AND follow-up beyond t for meaningful discrimination.
    if(sum(s$event==1 & s$time_months<=t)<5 || sum(s$time_months>t)<5) next
    roc <- tryCatch(timeROC::timeROC(T=s$time_months,delta=s$event,marker=risk,cause=1,
                                    times=t,iid=FALSE,weighting="marginal"),error=function(e) NULL)
    if(is.null(roc)) next
    ix <- which.min(abs(roc$times-t)); val <- roc$AUC[ix]
    if(length(val)!=1 || !is.finite(val)) next
    auc[[length(auc)+1L]] <- data.frame(cancer=ct,n=nrow(s),events=sum(s$event),year=t/12,
                                      time_months=t,AUC=val,marker="fitted_Cox_risk_score")
  }
}
if(!length(res)) stop("No cancer types meet discrimination criteria")
r <- data.table::rbindlist(res,fill=TRUE); r[,FDR:=p.adjust(p,"BH")]; write_table(r,"02_C_index_by_cancer.tsv")
a <- if(length(auc)) data.table::rbindlist(auc) else data.table::data.table(cancer=character(),year=numeric(),AUC=numeric())
write_table(a,"02_time_dependent_AUC_by_cancer.tsv")
if(nrow(a)) {
  p <- merge(r[,.(cancer,c_index,n)],a[year==3,.(cancer,AUC_3yr=AUC)],by="cancer",all=FALSE)
  write_table(p,"02_C_index_vs_3yr_AUC.tsv")
}
