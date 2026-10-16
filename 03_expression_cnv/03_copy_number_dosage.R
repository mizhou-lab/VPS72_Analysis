# VPS72 categorical copy number and continuous CNV-expression association.
# Reads the three original tab-delimited source exports; generates data tables only.
source("00_config.R")

amp <- fread(require_input(inputs$copy_number),check.names=FALSE)
cna <- fread(require_input(inputs$copy_expression),check.names=FALSE)
meta <- fread(require_input(inputs$copy_metadata),check.names=FALSE)
find_col <- function(d,preferred,pattern=NULL) {
  hit <- preferred[preferred %in% names(d)]
  if (!length(hit) && !is.null(pattern)) hit <- grep(pattern,names(d),ignore.case=TRUE,value=TRUE)
  if (!length(hit)) stop("Missing required column (",paste(preferred,collapse=" / "),")")
  hit[1L]
}
amp_id <- find_col(amp,c("Sample Id","sample","sample_id"))
cna_id <- find_col(cna,c("Sample Id","sample","sample_id"))
meta_id <- find_col(meta,c("sample","Sample Id","sample_id"))
ca_col <- find_col(meta,c("cancer type abbreviation","cancer_type_abbreviation","cancer_type"))
cat_col <- find_col(amp,c("Copy Number Alterations"))
log_col <- find_col(amp,c("VPS72: Log2 copy-number values"),"VPS72.*[Ll]og2")
ex_col <- find_col(cna,c(
  "VPS72: mRNA Expression, RSEM (Batch normalized from Illumina HiSeq_RNASeqV2)",
  "VPS72: mRNA Expression, RSEM (Batch normalized from Illumina HiSeq_RNASeqV2) (log2(value + 1))"
),"VPS72.*mRNA Expression")

# Metadata provides cancer labels. Only include primary tumors if the source
# contains complete TCGA sample-type barcodes. Never silently merge states.
make_keys <- function(x,label) {
  id15 <- tcga15(x)
  if (all(is.na(id15))) stop(label,": full TCGA sample barcodes are required to distinguish tumor and normal")
  if (any(is.na(id15))) warning(label,": discarding ",sum(is.na(id15))," samples with incomplete barcodes")
  id15
}
amp_dt <- data.table(sample=make_keys(amp[[amp_id]],"CNV source"),
                     category=gsub("^VPS72:[[:space:]]*","",trimws(as.character(amp[[cat_col]]))),
                     log2cnv=suppressWarnings(as.numeric(amp[[log_col]])))
cna_dt <- data.table(sample=make_keys(cna[[cna_id]],"Expression source"),
                     expression=suppressWarnings(as.numeric(cna[[ex_col]])))
meta_dt <- data.table(sample=make_keys(meta[[meta_id]],"Cancer metadata"),
                      cancer=trimws(as.character(meta[[ca_col]])))
amp_dt <- amp_dt[keep_primary(sample) & !is.na(category) & nzchar(category)]
cna_dt <- cna_dt[keep_primary(sample)]
meta_dt <- meta_dt[keep_primary(sample) & !is.na(cancer) & nzchar(cancer)]
for (pair in list(list(amp_dt,"CNV"),list(cna_dt,"Expression"),list(meta_dt,"Metadata"))) {
  if (anyDuplicated(pair[[1L]]$sample)) stop("Duplicate primary-tumor samples in ",pair[[2L]]," table; verify or deduplicate source data")
}
if (!nrow(amp_dt) || !nrow(cna_dt) || !nrow(meta_dt)) stop("Insufficient primary-tumor data")
joined <- merge(merge(amp_dt,cna_dt,by="sample"),meta_dt,by="sample")
if (!nrow(joined)) stop("No common primary tumor samples across the three inputs")
known <- c("Amplification","Gain","Diploid","Shallow Deletion","Deep Deletion")
if (any(!joined$category %in% known)) {
  other <- unique(joined$category[!joined$category %in% known])
  warning("Unexpected alteration labels excluded from frequency tables: ",paste(other,collapse=", "))
}
write_tsv(joined[,.(sample,cancer,category,log2cnv,expression)],"cnv_expression_matched_samples.tsv")

counts <- joined[category %in% known,.(n=.N),by=.(cancer,category)]
totals <- joined[category %in% known,.(total_samples=.N),by=cancer]
counts <- merge(counts,totals,by="cancer")
counts[,fraction:=n/total_samples]
setorder(counts,cancer,category)
write_tsv(counts,"per_cancer_cnv_categories.tsv")
pooled <- joined[category %in% known,.(n=.N),by=category]
pooled[,fraction:=n/sum(n)]
write_tsv(pooled,"pan_cancer_cnv_categories.tsv")

dosage <- joined[is.finite(log2cnv) & is.finite(expression), {
  r <- safe_spearman(log2cnv,expression,min_n=min_group_n)
  as.list(r)
},by=cancer]
dosage[,FDR:=p.adjust(p,method="BH")]
setorder(dosage,cancer)
write_tsv(dosage,"per_cancer_cnv_expression_spearman.tsv")
# As in the categorical source analysis, use the pooled expression median to
# dichotomize expression, and contrast diploid vs gain within each cancer.
# This is distinct from the continuous dosage correlation above.
median_all <- median(joined$expression,na.rm=TRUE)
joined[,expression_group:=fifelse(is.finite(expression),
   fifelse(expression>=median_all,"High","Low"),NA_character_)]
chisq_data <- joined[category %in% c("Diploid","Gain") & !is.na(expression_group)]
chisq_res <- if (nrow(chisq_data)) chisq_data[,{
  tab <- table(factor(category,levels=c("Diploid","Gain")),
               factor(expression_group,levels=c("Low","High")))
  valid <- all(rowSums(tab)>0) && all(colSums(tab)>0)
  chi <- if (valid) suppressWarnings(chisq.test(tab,correct=FALSE)) else NULL
  .(n=sum(tab),diploid_low=tab[1L,1L],diploid_high=tab[1L,2L],
    gain_low=tab[2L,1L],gain_high=tab[2L,2L],
    chi_square=if (!is.null(chi)) unname(chi$statistic) else NA_real_,
    p_value=if (!is.null(chi)) chi$p.value else NA_real_,
    low_expected=if (!is.null(chi)) any(chi$expected<5) else NA)
},by=cancer] else data.table(cancer=character(),n=integer(),diploid_low=integer(),
   diploid_high=integer(),gain_low=integer(),gain_high=integer(),
   chi_square=numeric(),p_value=numeric(),low_expected=logical())
chisq_res[,FDR:=p.adjust(p_value,method="BH")]
write_tsv(chisq_res,"per_cancer_cnv_expression_chisquare.tsv")
cat("Matched primary tumors:",nrow(joined),"; cancer groups:",uniqueN(joined$cancer),"\n")
