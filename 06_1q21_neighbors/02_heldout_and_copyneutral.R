# A4b held-out regional expression model and A5b copy-neutral vs normal analysis.
# Replaces invalid in-sample mean-residual tests and all-tumour-reference A5.
source("00_config.R")
check_packages("data.table")
library(data.table)
fp <- out_file("1q21_merged_long.rds")
check_required(fp)
dat <- readRDS(fp)

# Build shared wide sample matrix. The neighbour signature excludes the target gene.
zw <- dcast(dat, study+sample+cancer ~ gene, value.var="z_all")
cn <- dcast(dat, study+sample ~ gene, value.var="gistic")
setnames(cn,GENES_1Q21,paste0("CN_",GENES_1Q21))
zw <- merge(zw,cn,by=c("study","sample"),all=FALSE)

heldout <- rbindlist(lapply(GENES_1Q21,function(g) {
  # For VPS72 use 20 neighbours. For a neighbour use the other 19 neighbours,
  # excluding both itself and VPS72 from the signature.
  others <- setdiff(NEIGHBOURS,g)
  if (!length(others)) return(NULL)
  x <- copy(zw)
  x[,signature:=finite_row_mean(.SD),.SDcols=others]
  x[,n_signature_genes:=rowSums(is.finite(as.matrix(.SD))),.SDcols=others]
  x <- x[is.finite(signature) & n_signature_genes>=10L & is.finite(get(g)) &
           !is.na(get(paste0("CN_",g)))]
  ans <- rbindlist(lapply(sort(unique(x$cancer)),function(ct) {
    d<-x[cancer==ct]
    cn_g<-d[[paste0("CN_",g)]]
    train<-d[cn_g==0]; test<-d[cn_g>=1]
    if (nrow(train)<30L || nrow(test)<20L || sd(train$signature)==0) return(NULL)
    fit<-lm(y ~ signature,data=data.frame(y=train[[g]],signature=train$signature))
    pred<-predict(fit,newdata=data.frame(signature=test$signature))
    resid<-test[[g]]-pred
    resid<-resid[is.finite(resid)]
    if(length(resid)<20L || sd(resid)==0) return(NULL)
    tt<-t.test(resid,mu=0)
    data.table(gene=g,cancer=ct,n_train=nrow(train),n_test=length(resid),
               beta_neighbour=unname(coef(fit)[2]),mean_resid=mean(resid),
               median_resid=median(resid),ci_lo=tt$conf.int[1],ci_hi=tt$conf.int[2],
               p=tt$p.value)
  }),fill=TRUE)
  ans
}),fill=TRUE)
if (!nrow(heldout)) stop("Held-out models produced no valid tests; check sample coverage")
# BH correction across cancer types for each gene (including VPS72).
heldout[,fdr:=p.adjust(p,"BH"),by=gene]
setorder(heldout,gene,cancer)
fwrite(heldout,out_file("A4b_heldout_all_genes_by_cancer.tsv"),sep="\t")
fwrite(heldout[gene==GOI],out_file("A4b_heldout_positional_residuals.tsv"),sep="\t")
a4_rank<-heldout[is.finite(mean_resid),.(n_cancers=.N,median_resid=median(mean_resid),
                                        n_positive_fdr=sum(fdr<0.05 & mean_resid>0,na.rm=TRUE)),
                 by=gene][order(-median_resid)]
a4_rank[,rank:=seq_len(.N)]
fwrite(a4_rank,out_file("A4b_residual_ranking_all_genes.tsv"),sep="\t")

# A5b: in copy-neutral tumours, test expression relative to the provided normal
# reference. The normal-reference export is optional because not all cohorts have
# such scores; do not substitute z_all or z_dip when it is unavailable.
if (!"z_norm" %in% names(dat) || !any(is.finite(dat$z_norm))) {
  warning("A5b SKIPPED: no normal-reference z-scores. Provide z_norm input and rerun 01, then 02.")
} else {
  neutral<-dat[cn_group=="neutral" & is.finite(z_norm)]
  a5ct<-neutral[,{
    if(.N>=20L) {
      p<-tryCatch(suppressWarnings(wilcox.test(z_norm,mu=0,exact=FALSE)$p.value),
                  error=function(e) NA_real_)
      .(n=.N,median_z=median(z_norm),mean_z=mean(z_norm),pct_up=100*mean(z_norm>1),p=p)
    } else .(n=.N,median_z=NA_real_,mean_z=NA_real_,pct_up=NA_real_,p=NA_real_)
  },by=.(cancer,gene)]
  # Family: the 21 genes evaluated within each cancer type.
  a5ct[,fdr:=p.adjust(p,"BH"),by=cancer]
  fwrite(a5ct,out_file("A5b_copyneutral_vs_normal_by_cancer.tsv"),sep="\t")
  # Cancer types are weighted equally, not by their patient counts.
  a5pan<-a5ct[is.finite(median_z),
              .(n_cancers=.N,n_samples=sum(n),median_z=median(median_z),
                mean_z=mean(mean_z),pct_up=mean(pct_up),
                n_positive_fdr=sum(fdr<0.05 & median_z>0,na.rm=TRUE)),
              by=gene][order(-median_z)]
  a5pan[,rank:=seq_len(.N)]
  fwrite(a5pan,out_file("A5b_copyneutral_vs_normal.tsv"),sep="\t")
}
message("02 complete: held-out residuals and available copy-neutral/normal-reference results")
