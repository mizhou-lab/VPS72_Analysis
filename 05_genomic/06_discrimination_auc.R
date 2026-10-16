# ROC and PR AUC for the top-ranked genes from 05.
# In-sample: predictions come from the data the model was fitted on, with no
# resampling. The covariates-only AUC is fitted alongside, because the full
# model contains cancer type and mutation frequency varies sharply by lineage.
source("00_config.R")
suppressPackageStartupMessages({library(pROC);library(PRROC)})

master <- fread(require_input(inputs$master))
mut_mat <- fread(require_input(inputs$mutmat))
assoc <- fread(require_input(inputs$assoc))
if (!identical(master$patient,mut_mat$patient)) stop("Master and mutation tables are misaligned; rerun 01")

top_n <- 20L
genes <- intersect(head(assoc[order(p)]$gene,top_n),names(mut_mat))
if (!length(genes)) stop("No gene from the association table is in the mutation matrix")

auc_pair <- function(truth,preds) {
  roc <- try(pROC::roc(truth,preds,quiet=TRUE),silent=TRUE)
  pr <- try(PRROC::pr.curve(scores.class0=preds[truth==1L],scores.class1=preds[truth==0L]),silent=TRUE)
  c(roc=if (inherits(roc,"try-error")) NA_real_ else as.numeric(pROC::auc(roc)),
    pr=if (inherits(pr,"try-error")) NA_real_ else pr$auc.integral)
}

covars <- master[,.(VPS72_z,age,sex,cancer)]
score_gene <- function(g) {
  d <- data.table(mutated=mut_mat[[g]],covars)
  d <- d[complete.cases(d)]
  if (sum(d$mutated==1L)<min_auc_n) return(NULL)
  full <- try(glm(mutated~VPS72_z+age+sex+factor(cancer),data=d,family=binomial()),silent=TRUE)
  base <- try(glm(mutated~age+sex+factor(cancer),data=d,family=binomial()),silent=TRUE)
  if (inherits(full,"try-error") || inherits(base,"try-error")) return(NULL)
  a <- auc_pair(d$mutated,predict(full,type="response"))
  b <- auc_pair(d$mutated,predict(base,type="response"))
  data.table(gene=g,n_mutated=sum(d$mutated==1L),n_total=nrow(d),
             ROC_AUC=a[["roc"]],PR_AUC=a[["pr"]],
             ROC_AUC_covariates_only=b[["roc"]],
             ROC_AUC_gain_from_VPS72=a[["roc"]]-b[["roc"]])
}

res <- rbindlist(lapply(genes,score_gene))
if (!nrow(res)) stop("No gene reached ",min_auc_n," mutated samples")
setorder(res,-ROC_AUC)

write_tsv(res,"discrimination_auc.tsv")
print(res)
