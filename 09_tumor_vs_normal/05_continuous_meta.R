# Continuous VPS72-programme coupling, tumour versus normal:
#   1 per-cancer Spearman correlation within each tissue state
#   2 random-effects meta-analysis of those correlations, by state
#   3 meta-regression with tissue state as moderator (the formal difference)
#   4 per-cancer VPS72 x state interaction
#   5 random-effects meta-analysis of the interaction coefficients
#   6 paired-patient mixed-model sensitivity analysis
source("00_config.R")
suppressPackageStartupMessages({library(metafor)})

d <- fread(require_input(inputs$scores))
d <- d[cancer %chin% eligible & state %chin% c("Normal","Tumor")]
if (!nrow(d)) stop("No eligible samples in the score table")

score_cols <- c(if ("proliferation_score" %chin% names(d)) "proliferation_score",
                grep("^HALLMARK_",names(d),value=TRUE))
if (!length(score_cols)) stop("No pathway-score columns found")
cat("scores:",length(score_cols),"; cancers:",uniqueN(d$cancer),"\n")

# Correlations are pooled on the Fisher z scale, where the sampling variance is
# 1/(n-3) and does not depend on rho.
cor_one <- function(pw,ca,st) {
  z <- d[cancer==ca & state==st,.(VPS72,score=get(pw))]
  z <- z[is.finite(VPS72) & is.finite(score)]
  if (nrow(z)<min_cor_n || sd(z$VPS72)==0 || sd(z$score)==0) return(NULL)
  ct <- suppressWarnings(cor.test(z$VPS72,z$score,method="spearman",exact=FALSE))
  rho <- min(max(as.numeric(ct$estimate),-0.999999),0.999999)
  yi <- atanh(rho); se <- 1/sqrt(nrow(z)-3)
  data.table(pathway=pw,Program=pretty_name(pw),cancer=ca,state=st,
             n=nrow(z),rho=rho,p=ct$p.value,yi=yi,vi=1/(nrow(z)-3),
             ci_low=tanh(yi-1.96*se),ci_high=tanh(yi+1.96*se))
}

grid <- CJ(pw=score_cols,ca=eligible,st=c("Normal","Tumor"),sorted=FALSE)
cors <- rbindlist(Map(cor_one,grid$pw,grid$ca,grid$st))
if (!length(cors) || !nrow(cors)) stop("No cancer/state cell reached ",min_cor_n," samples")
cors[,FDR_within_state:=p.adjust(p,method="BH"),by=.(pathway,state)]
write_tsv(cors,"per_cancer_spearman.tsv")

meta_one <- function(pw,st) {
  z <- cors[pathway==pw & state==st]
  if (nrow(z)<min_meta_k) return(NULL)
  m <- tryCatch(rma.uni(yi=yi,vi=vi,data=z,method="REML",test="knha"),error=function(e) NULL)
  if (is.null(m)) return(NULL)
  pr <- predict(m,transf=tanh)   # back to the correlation scale for reporting
  data.table(pathway=pw,Program=pretty_name(pw),state=st,
             n_cancers=nrow(z),total_n=sum(z$n),
             pooled_rho=as.numeric(pr$pred),CI_low=as.numeric(pr$ci.lb),
             CI_high=as.numeric(pr$ci.ub),p=as.numeric(m$pval),
             tau2=as.numeric(m$tau2),I2=as.numeric(m$I2))
}
g2 <- CJ(pw=unique(cors$pathway),st=c("Normal","Tumor"),sorted=FALSE)
meta_state <- rbindlist(Map(meta_one,g2$pw,g2$st))
if (!length(meta_state) || !nrow(meta_state)) stop("No pathway reached ",min_meta_k," cancer types")
meta_state[,FDR:=p.adjust(p,method="BH"),by=state]
write_tsv(meta_state,"meta_by_state.tsv")

# Tumour and normal estimates from the same cancer are correlated, so the
# moderator model carries a random intercept per cancer rather than treating
# the two states as independent studies.
diff_one <- function(pw) {
  z <- copy(cors[pathway==pw])
  if (uniqueN(z$state)<2L || uniqueN(z$cancer)<min_meta_k) return(NULL)
  z[,state_num:=as.integer(state=="Tumor")]
  fit <- tryCatch(rma.mv(yi=yi,V=vi,mods=~state_num,random=~1|cancer,
                         data=z,method="REML",test="t"),error=function(e) NULL)
  if (is.null(fit)) return(NULL)
  i <- which(names(coef(fit))=="state_num")
  if (!length(i)) return(NULL)
  data.table(pathway=pw,Program=pretty_name(pw),n_cancers=uniqueN(z$cancer),
             fisher_z_difference=as.numeric(coef(fit)[i]),SE=as.numeric(fit$se[i]),
             CI_low=as.numeric(fit$ci.lb[i]),CI_high=as.numeric(fit$ci.ub[i]),
             p=as.numeric(fit$pval[i]))
}
diff_meta <- rbindlist(lapply(unique(cors$pathway),diff_one))
if (!length(diff_meta) || !nrow(diff_meta)) stop("Meta-regression produced no estimate")
diff_meta[,FDR:=p.adjust(p,method="BH")]
setorder(diff_meta,p)
write_tsv(diff_meta,"tumor_vs_normal_meta_regression.tsv")

# Within-state rank standardization makes each slope comparable to a Spearman
# rho, so the interaction term is a difference between two correlations.
int_one <- function(pw,ca) {
  z <- d[cancer==ca,.(state,patient,VPS72,score=get(pw))]
  z <- z[is.finite(VPS72) & is.finite(score)]
  if (uniqueN(z$state)<2L) return(NULL)
  if (any(z[,.N,by=state]$N<min_state_n)) return(NULL)
  z[,`:=`(x=rank_z(VPS72),y=rank_z(score)),by=state]
  z <- z[is.finite(x) & is.finite(y)]
  z[,stateF:=factor(state,levels=c("Normal","Tumor"))]
  fit <- tryCatch(lm(y~x*stateF,data=z),error=function(e) NULL)
  if (is.null(fit)) return(NULL)
  sm <- summary(fit)$coefficients
  rn <- grep("^x:stateF",rownames(sm),value=TRUE)
  if (!length(rn)) return(NULL)
  b <- sm[rn[1L],"Estimate"]; se <- sm[rn[1L],"Std. Error"]
  data.table(pathway=pw,Program=pretty_name(pw),cancer=ca,
             n_normal=sum(z$state=="Normal"),n_tumor=sum(z$state=="Tumor"),
             rho_normal=suppressWarnings(cor(z[state=="Normal"]$VPS72,z[state=="Normal"]$score,method="spearman")),
             rho_tumor=suppressWarnings(cor(z[state=="Tumor"]$VPS72,z[state=="Tumor"]$score,method="spearman")),
             interaction_beta=b,SE=se,CI_low=b-1.96*se,CI_high=b+1.96*se,
             p=sm[rn[1L],"Pr(>|t|)"])
}
g3 <- CJ(pw=score_cols,ca=eligible,sorted=FALSE)
interactions <- rbindlist(Map(int_one,g3$pw,g3$ca))
if (!length(interactions) || !nrow(interactions)) stop("No cancer supported an interaction model")
interactions[,FDR:=p.adjust(p,method="BH"),by=pathway]
write_tsv(interactions,"per_cancer_interactions.tsv")

imeta_one <- function(pw) {
  z <- interactions[pathway==pw & is.finite(interaction_beta) & SE>0]
  if (nrow(z)<min_meta_k) return(NULL)
  fit <- tryCatch(rma.uni(yi=interaction_beta,sei=SE,data=z,method="REML",test="knha"),
                  error=function(e) NULL)
  if (is.null(fit)) return(NULL)
  data.table(pathway=pw,Program=pretty_name(pw),n_cancers=nrow(z),
             pooled_interaction_beta=as.numeric(coef(fit)[1L]),
             CI_low=as.numeric(fit$ci.lb),CI_high=as.numeric(fit$ci.ub),
             p=as.numeric(fit$pval),tau2=as.numeric(fit$tau2),I2=as.numeric(fit$I2))
}
int_meta <- rbindlist(lapply(unique(interactions$pathway),imeta_one))
if (!length(int_meta) || !nrow(int_meta)) stop("Interaction meta-analysis produced no estimate")
int_meta[,FDR:=p.adjust(p,method="BH")]
setorder(int_meta,p)
write_tsv(int_meta,"interaction_meta.tsv")

# Sensitivity analysis restricted to patients contributing both a tumour and an
# adjacent normal, with a random intercept per patient. Skipped rather than
# approximated if lme4/lmerTest are absent.
if (!all(c("lme4","lmerTest") %in% rownames(installed.packages()))) {
  cat("paired mixed-model sensitivity skipped: install lme4 and lmerTest\n")
} else {
  pair_one <- function(pw,ca) {
    z <- d[cancer==ca & !is.na(patient) & paired==TRUE,.(state,patient,VPS72,score=get(pw))]
    z <- z[is.finite(VPS72) & is.finite(score)]
    ids <- z[,.(k=uniqueN(state)),by=patient][k==2L,patient]
    z <- z[patient %chin% ids]
    if (uniqueN(z$patient)<min_pairs) return(NULL)
    z[,`:=`(x=rank_z(VPS72),y=rank_z(score)),by=state]
    z <- z[is.finite(x) & is.finite(y)]
    z[,stateF:=factor(state,levels=c("Normal","Tumor"))]
    fit <- tryCatch(lmerTest::lmer(y~x*stateF+(1|patient),data=z,REML=FALSE),
                    error=function(e) NULL)
    if (is.null(fit)) return(NULL)
    sm <- coef(summary(fit))
    rn <- grep("^x:stateF",rownames(sm),value=TRUE)
    if (!length(rn)) return(NULL)
    data.table(pathway=pw,Program=pretty_name(pw),cancer=ca,
               n_pairs=uniqueN(z$patient),interaction_beta=sm[rn[1L],"Estimate"],
               SE=sm[rn[1L],"Std. Error"],p=sm[rn[1L],"Pr(>|t|)"],
               singular=lme4::isSingular(fit,tol=1e-4))
  }
  paired <- rbindlist(Map(pair_one,g3$pw,g3$ca))
  if (length(paired) && nrow(paired)) {
    paired[,FDR:=p.adjust(p,method="BH"),by=pathway]
    write_tsv(paired,"paired_patient_sensitivity.tsv")
  } else cat("no cancer reached",min_pairs,"paired patients\n")
}

cat("Composite proliferation, pooled rho by state:\n")
print(meta_state[pathway=="proliferation_score",
                 .(state,n_cancers,pooled_rho=round(pooled_rho,3),
                   CI_low=round(CI_low,3),CI_high=round(CI_high,3),I2=round(I2))])
