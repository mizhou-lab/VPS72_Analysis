# Logistic association of somatic mutation status with VPS72 expression:
#   mutated ~ VPS72_z + age + sex + cancer type
# BH correction is applied once across every gene tested, not per batch.
source("00_config.R")

master <- fread(require_input(inputs$master))
mut_mat <- fread(require_input(inputs$mutmat))
if (!identical(master$patient,mut_mat$patient)) stop("Master and mutation tables are misaligned; rerun 01")

all_genes <- setdiff(names(mut_mat),"patient")
if (file.exists(inputs$drivers)) {
  drivers <- toupper(unique(fread(inputs$drivers)[[1L]]))
  candidates <- intersect(all_genes,drivers)
  if (!length(candidates)) stop("No driver-list gene appears in the mutation matrix; check symbol spelling")
  cat("driver list:",length(drivers),"genes,",length(candidates),"present in the data\n")
} else {
  candidates <- all_genes
  cat("no driver list at",inputs$drivers,"- testing all",length(candidates),"genes\n")
}

# A logistic fit on a handful of mutated samples is not interpretable.
carriers <- colSums(mut_mat[,..candidates],na.rm=TRUE)
candidates <- names(carriers)[carriers>=min_mutated]
if (!length(candidates)) stop("No gene reaches ",min_mutated," carriers")

covars <- master[,.(VPS72_z,age,sex,cancer)]
fit_one <- function(g) {
  d <- data.table(mutated=mut_mat[[g]],covars)
  d <- d[complete.cases(d)]
  if (uniqueN(d$mutated)<2L) return(NULL)
  m <- try(glm(mutated~VPS72_z+age+sex+factor(cancer),data=d,family=binomial()),silent=TRUE)
  if (inherits(m,"try-error") || !isTRUE(m$converged)) return(NULL)
  tid <- as.data.table(tidy(m,conf.int=TRUE,exponentiate=TRUE))[term=="VPS72_z"]
  if (!nrow(tid)) return(NULL)
  data.table(gene=g,n_mutated=sum(d$mutated==1L),n_total=nrow(d),
             OR=tid$estimate,CI_low=tid$conf.low,CI_high=tid$conf.high,p=tid$p.value)
}

res <- rbindlist(lapply(candidates,fit_one))
if (!nrow(res)) stop("No model converged")
res[,FDR:=p.adjust(p,method="BH")]
res[,direction:=fifelse(OR>1,"enriched_in_VPS72_high","depleted_in_VPS72_high")]
setorder(res,p)

write_tsv(res,"mutation_association.tsv")
cat("Genes tested:",length(candidates),"; converged:",nrow(res),
    "; at FDR <",fdr_cut,":",sum(res$FDR<fdr_cut,na.rm=TRUE),"\n")
