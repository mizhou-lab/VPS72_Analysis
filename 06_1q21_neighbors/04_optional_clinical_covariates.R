# OPTIONAL sensitivity analyses: survival adjusted for MCL1/SETDB1, and
# proliferation adjusted for regional CN, FGA/aneuploidy and true 1q-arm gain.
# Requires the optional sample-level files listed in 00_config.R.
source("00_config.R")
check_packages(c("data.table","survival"))
library(data.table)
library(survival)
if (!file.exists(out_file("1q21_regional_primary_tumors.rds")))
  stop("Run 03_regional_adjustment.R first")
dat <- readRDS(out_file("1q21_regional_primary_tumors.rds"))
find_col <- function(nms, opts) {
  fmt<-function(x)tolower(gsub("[^a-z0-9]","",x))
  hit<-match(fmt(opts),fmt(nms)); hit<-hit[!is.na(hit)]
  if(!length(hit)) return(NA_character_)
  nms[hit[1]]
}
parse_event <- function(x) {
  if(is.numeric(x)) return(as.integer(x>0))
  z<-toupper(trimws(as.character(x)))
  as.integer(ifelse(z %in% c("1","1:DECEASED","DECEASED","DEAD","EVENT"),1L,
                    ifelse(z %in% c("0","0:LIVING","LIVING","ALIVE","CENSORED"),0L,NA_integer_)))
}
coef_row <- function(model,term,label) {
  s<-summary(model)$coefficients
  if(!term %in% rownames(s)) return(NULL)
  r<-s[term,]
  data.table(model=label,term=term,HR=exp(r["coef"]),CI_low=exp(r["coef"]-1.96*r["se(coef)"]),
             CI_high=exp(r["coef"]+1.96*r["se(coef)"]),p=unname(r["Pr(>|z|)"]),
             n=model$n,events=model$nevent)
}

if(file.exists(FILES_1Q21$clinical)) {
  clin<-fread(FILES_1Q21$clinical)
  idc<-find_col(names(clin),c("PATIENT_ID","patient","bcr_patient_barcode","SAMPLE_ID"))
  tc<-find_col(names(clin),c("OS_MONTHS","OS_TIME","OS.time"))
  ec<-find_col(names(clin),c("OS_STATUS","OS_EVENT","vital_status"))
  if(anyNA(c(idc,tc,ec))) {
    warning("Optional clinical table missing recognized patient/OS time/OS event columns; survival skipped")
  } else {
    clin[,patient:=patient_id(get(idc))]
    clin[,os_time:=suppressWarnings(as.numeric(get(tc)))]
    clin[,os_event:=parse_event(get(ec))]
    clin<-clin[!is.na(patient),.(patient,os_time,os_event)]
    if(anyDuplicated(clin$patient)) stop("Duplicate patients in survival file; resolve rather than selecting arbitrarily")
    s<-merge(dat[,.(patient,cancer,exprz_VPS72,exprz_MCL1,exprz_SETDB1)],
             clin,by="patient",all=FALSE)
    s<-s[complete.cases(s) & os_time>0]
    if(nrow(s)>=100 && sum(s$os_event==1)>=20 && uniqueN(s$cancer)>=2) {
      m1<-survival::coxph(survival::Surv(os_time,os_event)~exprz_VPS72+strata(cancer),data=s)
      m2<-survival::coxph(survival::Surv(os_time,os_event)~exprz_VPS72+exprz_MCL1+exprz_SETDB1+strata(cancer),data=s)
      res<-rbindlist(list(coef_row(m1,"exprz_VPS72","VPS72 + cancer strata"),
                          coef_row(m2,"exprz_VPS72","+ MCL1 + SETDB1")),fill=TRUE)
      fwrite(res,out_file("S1_cox_pancancer.tsv"),sep="\t")
      by_ct<-rbindlist(lapply(sort(unique(s$cancer)),function(ct) {
        d<-s[cancer==ct]
        if(nrow(d)<60L || sum(d$os_event==1)<20L) return(NULL)
        fit<-try(survival::coxph(survival::Surv(os_time,os_event)~exprz_VPS72+exprz_MCL1+exprz_SETDB1,data=d),silent=TRUE)
        if(inherits(fit,"try-error")) return(NULL)
        out<-coef_row(fit,"exprz_VPS72",ct)
        if(is.null(out)) return(NULL)
        out[,cancer:=ct]
        out
      }),fill=TRUE)
      if(nrow(by_ct)) by_ct[,fdr:=p.adjust(p,"BH")]
      fwrite(by_ct,out_file("S1_cox_by_cancer.tsv"),sep="\t")
    } else warning("Too few evaluable events/samples for optional Cox model")
  }
} else message("No optional clinical file; Cox module skipped")

if(file.exists(FILES_1Q21$covariates)) {
  cov<-fread(FILES_1Q21$covariates)
  idc<-find_col(names(cov),c("PATIENT_ID","patient","bcr_patient_barcode","SAMPLE_ID"))
  pc<-find_col(names(cov),c("PROLIFERATION_SCORE","ProliferationScore","cell_cycle_score"))
  if(anyNA(c(idc,pc))) {
    warning("Optional covariates missing patient or proliferation score; module skipped")
  } else {
    cov[,patient:=patient_id(get(idc))]
    cov[,proliferation:=suppressWarnings(as.numeric(get(pc)))]
    ac<-find_col(names(cov),c("ANEUPLOIDY_SCORE","AneuploidyScore"))
    fc<-find_col(names(cov),c("FGA","fraction_genome_altered"))
    qc<-find_col(names(cov),c("ARM_1Q_GAIN","1Q_GAIN","chr1q_gain"))
    cov[,aneuploidy:=if(!is.na(ac)) suppressWarnings(as.numeric(get(ac))) else NA_real_]
    cov[,FGA:=if(!is.na(fc)) suppressWarnings(as.numeric(get(fc))) else NA_real_]
    cov[,arm1q_gain:=if(!is.na(qc)) suppressWarnings(as.numeric(get(qc))) else NA_real_]
    cov<-cov[!is.na(patient),.(patient,proliferation,aneuploidy,FGA,arm1q_gain)]
    if(anyDuplicated(cov$patient)) stop("Duplicate patients in covariate file; resolve upstream")
    d<-merge(dat[,.(patient,cancer,exprz_VPS72,exprz_MCL1,exprz_SETDB1,regional_cn_z)],
             cov,by="patient",all=FALSE)
    d[,proliferation_z:=zsafe(proliferation),by=cancer]
    d[,aneuploidy_z:=zsafe(aneuploidy),by=cancer]
    d[,FGA_z:=zsafe(FGA),by=cancer]
    formulas<-list(
      P1="proliferation_z ~ exprz_VPS72 + factor(cancer)",
      P2="proliferation_z ~ exprz_VPS72 + regional_cn_z + factor(cancer)",
      P3="proliferation_z ~ exprz_VPS72 + regional_cn_z + aneuploidy_z + FGA_z + factor(cancer)",
      P4="proliferation_z ~ exprz_VPS72 + regional_cn_z + arm1q_gain + factor(cancer)",
      P5="proliferation_z ~ exprz_VPS72 + exprz_MCL1 + exprz_SETDB1 + regional_cn_z + factor(cancer)"
    )
    # Each sensitivity model reports its actual complete-case sample size.
    tests<-rbindlist(lapply(names(formulas),function(label) {
      f<-as.formula(formulas[[label]])
      cols<-all.vars(f)
      sub<-d[complete.cases(d[,..cols])]
      if(nrow(sub)<100L || uniqueN(sub$cancer)<2L || sd(sub$proliferation_z)==0) return(NULL)
      fit<-try(lm(f,data=sub),silent=TRUE)
      if(inherits(fit,"try-error")) return(NULL)
      co<-summary(fit)$coefficients
      if(!"exprz_VPS72" %in% rownames(co)) return(NULL)
      a<-co["exprz_VPS72",]
      data.table(model=label,covariates=formulas[[label]],n=nrow(sub),
                 beta=unname(a["Estimate"]),se=unname(a["Std. Error"]),
                 p=unname(a["Pr(>|t|)"]),adj_r2=summary(fit)$adj.r.squared)
    }),fill=TRUE)
    if(nrow(tests)) fwrite(tests,out_file("S2_proliferation_adjusted.tsv"),sep="\t")
  }
} else message("No optional covariate file; proliferation sensitivity module skipped")
message("04 complete: available optional survival / proliferation sensitivity results")
