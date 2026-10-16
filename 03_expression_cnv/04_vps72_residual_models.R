# VPS72 association beyond SRCAP/TIP60 modules and optional GISTIC CNV.
# Uses actual aligned primary-tumor samples, not an assumed intermediate file.
# Requires 01_primary_tumor_signatures.R to have completed successfully.
source("00_config.R")

exp_fp <- file.path(output_dir,"model_expression_samples.rds")
cnv_fp <- file.path(output_dir,"model_gistic_samples.rds")
man_fp <- file.path(output_dir,"primary_tumor_manifest.tsv")
expr <- readRDS(require_input(exp_fp))
cnv <- readRDS(require_input(cnv_fp))
manifest <- fread(require_input(man_fp))
needed <- unique(c("VPS72",srcap_genes,tip60_genes))
missing <- setdiff(needed,rownames(expr))
if (length(missing)) stop("Expression matrix lacks residual model members: ",paste(missing,collapse=", "))
if (!"VPS72" %in% rownames(cnv)) stop("GISTIC VPS72 row missing")
if (!identical(colnames(expr),colnames(cnv))) stop("Expression/CNV sample alignment mismatch")
if (!identical(colnames(expr),manifest$sample)) stop("Sample manifest does not match matrices")

# Standardize each gene across primary tumors, then average module genes.
standardize <- function(m) {
  sd <- matrixStats::rowSds(m,na.rm=TRUE)
  sd[sd==0 | !is.finite(sd)] <- NA_real_
  sweep(sweep(m,1L,rowMeans(m,na.rm=TRUE),"-"),1L,sd,"/")
}
z <- standardize(expr)
model <- data.table(sample=colnames(z),cohort=manifest$cohort,
  VPS72=as.numeric(z["VPS72",]),
  SRCAP_score=colMeans(z[srcap_genes,,drop=FALSE],na.rm=TRUE),
  TIP60_score=colMeans(z[tip60_genes,,drop=FALSE],na.rm=TRUE),
  VPS72_cnv=as.numeric(cnv["VPS72",]))
for (nm in c("VPS72","SRCAP_score","TIP60_score","VPS72_cnv"))
  model <- model[is.finite(get(nm))]
if (nrow(model)<25L || uniqueN(model$cohort)<4L)
  stop("Insufficient complete primary tumor data for residual modeling")

# Distinct nested models. Cancer-specific residual means are calculated from
# pan-cancer fits WITHOUT cancer as a predictor; otherwise these means are forced
# close to zero by design.
model_specs <- list(
  complexes_only=VPS72 ~ SRCAP_score + TIP60_score,
  complexes_plus_CNV=VPS72 ~ SRCAP_score + TIP60_score + VPS72_cnv
)
model_fits <- lapply(model_specs,function(f) lm(f,data=model))
coefs <- rbindlist(lapply(names(model_fits),function(nm) {
  fit <- model_fits[[nm]]
  terms <- as.data.table(summary(fit)$coefficients,keep.rownames="term")
  setnames(terms,c("Estimate","Std. Error","t value","Pr(>|t|)"),c("estimate","SE","t_stat","p_value"))
  terms[,`:=`(model=nm,n=nobs(fit),R2=summary(fit)$r.squared)]
  terms
}),fill=TRUE)
write_tsv(coefs,"residual_model_coefficients.tsv")

sample_residuals <- copy(model[,.(sample,cohort,VPS72,VPS72_cnv)])
for (nm in names(model_fits)) sample_residuals[,(paste0("residual_",nm)):=residuals(model_fits[[nm]])]
write_tsv(sample_residuals,"vps72_sample_residuals.tsv")

cancer_results <- rbindlist(lapply(names(model_fits),function(nm) {
  r <- model[,.(cohort)]
  r[,residual:=as.numeric(residuals(model_fits[[nm]]))]
  out <- r[,.(n=.N,mean_residual=mean(residual),median_residual=median(residual),
                SD=sd(residual)),by=cohort]
  out[,`:=`(model=nm,SE=SD/sqrt(n))]
  out[,p_value:=ifelse(n>=3L & is.finite(SE) & SE>0,
                        2*pt(-abs(mean_residual/SE),df=n-1L),NA_real_)]
  out[,FDR:=p.adjust(p_value,method="BH")]
  out
}),fill=TRUE)
write_tsv(cancer_results,"vps72_per_cancer_residuals.tsv")

# A separate descriptive regression on cancer-type means (equal weight/cancer).
cancer_level <- model[,.(VPS72=median(VPS72),SRCAP_score=median(SRCAP_score),
                         TIP60_score=median(TIP60_score),VPS72_cnv=mean(VPS72_cnv)),by=cohort]
if (nrow(cancer_level)>=8L) {
  for (nm in names(model_specs)) {
    fit <- lm(model_specs[[nm]],data=cancer_level)
    cancer_level[,(paste0("cancer_residual_",nm)):=as.numeric(residuals(fit))]
  }
  write_tsv(cancer_level,"vps72_cancer_level_model.tsv")
}
cat("Residual models fitted to",nrow(model),"primary tumors in",uniqueN(model$cohort),"cohorts\n")
