# TCGA tumor-versus-adjacent-normal expression and SRCAP/TIP60/INO80 modules.
# Uses TPM (not the EB++ primary-tumor matrix); retains tissue-state barcodes.
source("00_config.R")

required <- unique(c("VPS72",unlist(comparison_modules,use.names=FALSE)))
mat <- load_gene_matrix(inputs$tpm_expression,required,"TPM")
mat <- mat[,keep_primary(colnames(mat)) | keep_normal(colnames(mat)),drop=FALSE]
mat <- collapse_samples(mat)
if (!"VPS72" %in% rownames(mat)) stop("VPS72 missing from TPM matrix")
if (any(mat<0,na.rm=TRUE)) stop("Expected untransformed, nonnegative TPM expression")

pheno <- fread(require_input(inputs$tpm_phenotype))
find_pheno <- function(patterns,required=TRUE) {
  for (pat in patterns) {
    n <- grep(pat,names(pheno),ignore.case=TRUE,value=TRUE)
    if (length(n)) return(n[1L])
  }
  if (required) stop("Required phenotype field not found: ",paste(patterns,collapse=" / "))
  NA_character_
}
samp_col <- find_pheno(c("^sample$","^sample_id$","^samples$","_sample$","^sample"))
ca_col <- find_pheno(c("cancer.type.abbreviation","^_primary_disease$","^primary_disease$",
                       "^cancer_type$","^cohort$","^project$"))
stage_col <- find_pheno(c("^pathologic_stage$","^ajcc_pathologic_tumor_stage$","^clinical_stage$","stage"),FALSE)
ids <- tcga15(pheno[[samp_col]])
ph <- data.table(sample=ids,patient=tcga12(pheno[[samp_col]]),
                 cancer=trimws(as.character(pheno[[ca_col]])),
                 stage=if (!is.na(stage_col)) as.character(pheno[[stage_col]]) else NA_character_)
ph <- ph[!is.na(patient) & !is.na(cancer) & nzchar(cancer)]
# Patient cancer type must be unambiguous; map metadata by patient to avoid
# introducing duplicate tissue-state rows from a phenotype export.
patient_cancer <- unique(ph[,.(patient,cancer)])
if (patient_cancer[,uniqueN(cancer),by=patient][V1>1L,.N]>0L)
  stop("Phenotype contains conflicting cancer labels for one patient")
patient_cancer <- unique(patient_cancer,by="patient")
manifest <- data.table(sample=colnames(mat),patient=substr(colnames(mat),1L,12L))
manifest[,cancer:=patient_cancer$cancer[match(patient,patient_cancer$patient)]]
manifest[,state:=fifelse(substr(sample,14L,15L)=="01","Tumor","Normal")]
manifest <- manifest[!is.na(cancer)]
if (!nrow(manifest)) stop("No annotated tumor/normal TCGA samples")
# Use sample-specific stage if available; it is never inferred from another tissue state.
ph_stage <- ph[!is.na(sample) & !is.na(stage) & nzchar(stage),.(sample,stage)]
if (ph_stage[,uniqueN(stage),by=sample][V1>1L,.N]>0L) stop("Conflicting stage annotations")
ph_stage <- unique(ph_stage,by="sample")
manifest[,stage:=ph_stage$stage[match(sample,ph_stage$sample)]]
mat <- mat[,manifest$sample,drop=FALSE]
write_tsv(manifest,"tumor_normal_sample_manifest.tsv")

# Numeric TPM per individual TCGA sample, not log-transformed or collapsed by patient.
targets <- intersect(required,rownames(mat))
long <- data.table(gene=rep(targets,times=ncol(mat)),
                   sample=rep(colnames(mat),each=length(targets)),
                   expression=as.vector(mat[targets,,drop=FALSE]))
long <- merge(long,manifest,by="sample",sort=FALSE)
write_tsv(long[gene=="VPS72"],"vps72_tumor_normal_samples.tsv")

# Comparative log2 fold change of MEAN TPM, with Wilcoxon on untransformed TPM.
compare_states <- function(d,by_cols) {
  out <- d[,{
    t <- expression[state=="Tumor" & is.finite(expression)]
    n <- expression[state=="Normal" & is.finite(expression)]
    p <- if (length(t)>=min_tumor_normal_n && length(n)>=min_tumor_normal_n)
      tryCatch(wilcox.test(t,n,exact=FALSE)$p.value,error=function(e) NA_real_) else NA_real_
    .(n_tumor=length(t),n_normal=length(n),mean_tumor=if (length(t)) mean(t) else NA_real_,
      mean_normal=if (length(n)) mean(n) else NA_real_,p_value=p)
  },by=by_cols]
  out[,log2FC:=log2(mean_tumor+1)-log2(mean_normal+1)]
  out[,FDR:=p.adjust(p_value,method="BH")]
  out
}
gene_tests <- compare_states(long,c("gene","cancer"))
write_tsv(gene_tests,"gene_tumor_normal_tests.tsv")
write_tsv(gene_tests[gene=="VPS72"],"vps72_tumor_normal_tests.tsv")

module_data <- rbindlist(lapply(names(comparison_modules),function(name) {
  genes <- intersect(comparison_modules[[name]],rownames(mat))
  if (length(genes)<3) stop("Module ",name," has <3 measured genes")
  scores <- colMeans(mat[genes,,drop=FALSE],na.rm=TRUE)
  data.table(sample=colnames(mat),module=name,expression=as.numeric(scores),
             genes_found=length(genes))
}))
module_data <- merge(module_data,manifest,by="sample",sort=FALSE)
write_tsv(module_data,"module_expression_samples.tsv")
write_tsv(compare_states(module_data,c("module","cancer")),"module_tumor_normal_tests.tsv")

# Pathologic-stage comparisons are evaluated WITHIN cancer type and in tumors only.
parse_stage <- function(x) {
  x <- toupper(trimws(as.character(x)))
  z <- sub("^STAGE[[:space:]]*","",x)
  z <- sub("[ABC]$","",z)
  z[z %in% c("I","II","III","IV")] <- as.character(match(z[z %in% c("I","II","III","IV")],c("I","II","III","IV")))
  suppressWarnings(as.integer(ifelse(z %in% as.character(1:4),z,NA_character_)))
}
stage_sources <- rbindlist(list(
  long[gene=="VPS72",.(sample,cancer,state,stage,measure="VPS72",value=expression)],
  module_data[,.(sample,cancer,state,stage,measure=module,value=expression)]
))
stage_sources[,stage_int:=parse_stage(stage)]
stage_sources <- stage_sources[state=="Tumor" & is.finite(value) & !is.na(stage_int)]
stage_tests <- if (nrow(stage_sources)) {
  stage_sources[,{
    if (.N < min_stage_n || uniqueN(stage_int)<2L) .(n=.N,rho=NA_real_,p=NA_real_)
    else {r<-safe_spearman(value,stage_int,min_n=min_stage_n); as.list(r)}
  },by=.(cancer,measure)]
} else data.table(cancer=character(),measure=character(),n=integer(),rho=numeric(),p=numeric())
stage_tests[,FDR:=p.adjust(p,method="BH")]
write_tsv(stage_tests,"vps72_module_stage_spearman.tsv")
cat("Tumor/normal samples:",ncol(mat),"; cancers:",uniqueN(manifest$cancer),"\n")
