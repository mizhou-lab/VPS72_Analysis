# Regional 1q21 analyses from raw cBioPortal RSEM and continuous CN values.
# Separate from scripts 01-02, which use cBioPortal expression z-scores.
# Regional 1q21 summaries are NOT chromosome-arm 1q gain, FGA or aneuploidy.
source("00_config.R")
check_packages("data.table")
library(data.table)
check_required(unlist(FILES_1Q21[c("gistic","log2cn","rsem")]))

cg<-read_cbio_wide(FILES_1Q21$gistic,"gistic_")
cc<-read_cbio_wide(FILES_1Q21$log2cn,"cn_")
ce<-read_cbio_wide(FILES_1Q21$rsem,"rsem_")
dat<-merge(merge(cg,cc,by=c("study","sample")),ce,by=c("study","sample"))
if (!nrow(dat)) stop("No shared primary tumour samples in raw RSEM/CN exports")
dat[,cancer:=toupper(sub("_tcga.*$","",study))]
dat[,patient:=patient_id(sample)]
setorder(dat,study,patient,sample)
dat<-dat[!duplicated(dat[,.(study,patient)])]
for(g in GENES_1Q21) {
  rsem<-paste0("rsem_",g)
  if(any(dat[[rsem]] < 0,na.rm=TRUE)) stop("Negative RSEM values for ",g,
                                         "; expected non-negative unlogged values")
  lg<-paste0("logexpr_",g); ze<-paste0("exprz_",g)
  dat[,(lg):=log2(get(rsem)+1)]
  dat[,(ze):=zsafe(get(lg)),by=cancer]
}
regional_cn_cols<-paste0("cn_",NEIGHBOURS)
regional_expr_cols<-paste0("exprz_",NEIGHBOURS)
regional_gistic_cols<-paste0("gistic_",NEIGHBOURS)
dat[,regional_cn_mean:=finite_row_mean(.SD),.SDcols=regional_cn_cols]
dat[,regional_expr_mean:=finite_row_mean(.SD),.SDcols=regional_expr_cols]
dat[,regional_gain_fraction:=finite_row_mean(as.data.frame(lapply(.SD,function(z) as.numeric(z>=1)))),
    .SDcols=regional_gistic_cols]
dat[,regional_amp_fraction:=finite_row_mean(as.data.frame(lapply(.SD,function(z) as.numeric(z==2)))),
    .SDcols=regional_gistic_cols]
dat[,regional_cn_z:=zsafe(regional_cn_mean),by=cancer]
dat[,regional_expr_z:=zsafe(regional_expr_mean),by=cancer]
dat[,VPS72_cn_resid:={
  v<-cn_VPS72; reg<-regional_cn_mean
  ok<-is.finite(v)&is.finite(reg)
  ans<-rep(NA_real_,.N)
  if(sum(ok)>=20L && sd(reg[ok])>0) ans[ok]<-resid(lm(v[ok]~reg[ok]))
  ans
},by=cancer]
dat[,VPS72_cn_resid_z:=zsafe(VPS72_cn_resid),by=cancer]
dat[,VPS72_gain:=gistic_VPS72>=1]
dat[,broad_regional_gain:=regional_gain_fraction>=0.8]
dat[,gain_architecture:=fcase(is.na(VPS72_gain)|is.na(broad_regional_gain),NA_character_,
                               !VPS72_gain,"VPS72 not gained",
                               VPS72_gain & broad_regional_gain,"VPS72 + broad regional co-gain",
                               default="VPS72 gain with <80% regional co-gain")]
saveRDS(dat,out_file("1q21_regional_primary_tumors.rds"))
qc<-data.table(metric=c("primary_tumor_patients","cancer_types","genes_1q21"),
               value=c(nrow(dat),uniqueN(dat$cancer),length(GENES_1Q21)))
fwrite(qc,out_file("R0_regional_QC.tsv"),sep="\t")
arch<-dat[VPS72_gain==TRUE & !is.na(gain_architecture),
          .(n=.N),by=.(cancer,gain_architecture)]
arch[,pct:=100*n/sum(n),by=cancer]
fwrite(arch,out_file("R1_gain_architecture_by_cancer.tsv"),sep="\t")

# Matched-observation nested models, within cancer:
# M0: VPS72 expression ~ regional expression
# M1: M0 + regional CN; M2: M1 + VPS72-local CN residual.
reg_models<-rbindlist(lapply(sort(unique(dat$cancer)),function(ct) {
  d<-dat[cancer==ct,.(exprz_VPS72,regional_expr_z,regional_cn_z,VPS72_cn_resid_z)]
  d<-d[complete.cases(d)]
  if(nrow(d)<50L || any(vapply(d,function(z) sd(z)==0,logical(1)))) return(NULL)
  m0<-lm(exprz_VPS72~regional_expr_z,data=d)
  m1<-lm(exprz_VPS72~regional_expr_z+regional_cn_z,data=d)
  m2<-lm(exprz_VPS72~regional_expr_z+regional_cn_z+VPS72_cn_resid_z,data=d)
  cs<-summary(m2)$coefficients
  if(!"VPS72_cn_resid_z" %in% rownames(cs)) return(NULL)
  row<-cs["VPS72_cn_resid_z",]
  data.table(cancer=ct,n=nrow(d),R2_M0=summary(m0)$r.squared,
             R2_M1=summary(m1)$r.squared,R2_M2=summary(m2)$r.squared,
             delta_R2_regional_CN=summary(m1)$r.squared-summary(m0)$r.squared,
             delta_R2_local_CN=summary(m2)$r.squared-summary(m1)$r.squared,
             beta_local_CN=unname(row["Estimate"]),se_local_CN=unname(row["Std. Error"]),
             p_local_CN=unname(row["Pr(>|t|)"]))
}),fill=TRUE)
if(nrow(reg_models)) reg_models[,fdr_local_CN:=p.adjust(p_local_CN,"BH")]
fwrite(reg_models,out_file("R2_regional_adjustment_by_cancer.tsv"),sep="\t")
d<-dat[,.(cancer,exprz_VPS72,regional_expr_z,regional_cn_z,VPS72_cn_resid_z)]
d<-d[complete.cases(d)]
if(nrow(d)>100L && uniqueN(d$cancer)>1L) {
  pooled<-lm(exprz_VPS72~regional_expr_z+regional_cn_z+VPS72_cn_resid_z+factor(cancer),data=d)
  sm<-summary(pooled)$coefficients
  pooled_res<-data.table(term=rownames(sm),estimate=sm[,1],se=sm[,2],t=sm[,3],p=sm[,4],
                          n=nrow(d),model_r2=summary(pooled)$r.squared)
  fwrite(pooled_res,out_file("R2_regional_adjustment_pancancer.tsv"),sep="\t")
}

# Descriptive expression collinearity relative to other important 1q21 genes.
cors<-rbindlist(lapply(PRIORITY,function(g) {
  v<-dat[[paste0("exprz_",g)]]; y<-dat$exprz_VPS72
  ok<-is.finite(v)&is.finite(y)
  if(sum(ok)<10L || sd(v[ok])==0 || sd(y[ok])==0)
    return(data.table(gene=g,n=sum(ok),rho=NA_real_,p=NA_real_))
  test<-suppressWarnings(cor.test(v[ok],y[ok],method="spearman",exact=FALSE))
  data.table(gene=g,n=sum(ok),rho=unname(test$estimate),p=test$p.value)
}))
cors[,fdr:=p.adjust(p,"BH")]
fwrite(cors,out_file("R3_expression_correlations_priority_genes.tsv"),sep="\t")

vif_cols<-c("exprz_VPS72","exprz_MCL1","exprz_SETDB1","regional_cn_z")
x<-dat[,..vif_cols]; x<-x[complete.cases(x)]
vif<-rbindlist(lapply(vif_cols,function(v) {
  others<-setdiff(vif_cols,v)
  if(nrow(x)<30L || any(vapply(x, function(z) sd(z)==0,logical(1))))
    return(data.table(term=v,VIF=NA_real_,n=nrow(x)))
  mod<-lm(reformulate(others,response=v),data=x)
  r2<-summary(mod)$r.squared
  data.table(term=v,VIF=if(is.finite(r2)&&r2<1)1/(1-r2) else NA_real_,n=nrow(x))
}))
fwrite(vif,out_file("R3_collinearity_VIF.tsv"),sep="\t")
message("03 complete: regional 1q21 analyses; proxy NOT equivalent to whole-arm gain")
