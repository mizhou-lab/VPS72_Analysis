# 1q21/VPS72 analysis: shared configuration (no figures or spreadsheet exports)
# Paths are relative to the repository working directory; override with env vars.
DATA_DIR <- Sys.getenv("VPS72_1Q21_DATA", unset = "data/1q21")
OUT_DIR <- Sys.getenv("VPS72_1Q21_OUT", unset = "results/1q21")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

GOI <- "VPS72"
NEIGHBOURS <- c("CHD1L","PRKAB2","FMO5","BCL9","H2AC18","ANP32E","MCL1","ECM1",
                "ENSA","ARNT","SETDB1","CERS2","MLLT11","PRUNE1","PIP5K1A",
                "PSMD4","ZNF687","S100A4","ILF2","GATAD2B")
GENES_1Q21 <- c(NEIGHBOURS, GOI)
PRIORITY <- c("MCL1","SETDB1","CHD1L","PSMD4","ILF2","PRUNE1")
MIN_GROUP <- 10L

FILES_1Q21 <- list(
  gistic = file.path(DATA_DIR, "cna.txt"),
  log2cn = file.path(DATA_DIR, "Log2 copy-number values.txt"),
  z_dip = file.path(DATA_DIR, "mRNA expression z-scores relative to diploid samples (RNA Seq V2 RSEM).txt"),
  z_all = file.path(DATA_DIR, "mRNA expression z-scores relative to all samples (log RNA Seq V2 RSEM).txt"),
  z_norm = file.path(DATA_DIR, "mRNA expression z-scores relative to normal samples (log RNA Seq V2 RSEM).txt"),
  rsem = file.path(DATA_DIR, "mRNA Expression, RSEM (Batch normalized from Illumina HiSeq_RNASeqV2) (1).txt"),
  clinical = file.path(DATA_DIR, "TCGA_survival_clinical.tsv"),
  covariates = file.path(DATA_DIR, "TCGA_sample_covariates.tsv")
)
# If cBioPortal downloads use '(1)' suffixes, edit the corresponding FILES_1Q21 entry.

out_file <- function(stem) file.path(OUT_DIR, stem)
check_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required R packages: ", paste(missing, collapse = ", "))
}
check_required <- function(files) {
  missing <- files[!file.exists(files)]
  if (length(missing)) stop("Missing required inputs:\n", paste(missing, collapse = "\n"))
}
# Require an explicit sample-type barcode. No patient-level aggregation of tumour + normal.
is_primary <- function(ids) {
  s <- toupper(gsub("\\.", "-", trimws(as.character(ids))))
  grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-01", s)
}
patient_id <- function(ids) substr(toupper(gsub("\\.", "-", as.character(ids))), 1L, 12L)
zsafe <- function(x) {
  x <- as.numeric(x); n <- sum(is.finite(x)); s <- stats::sd(x, na.rm = TRUE)
  if (n < 3L || !is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE))/s
}
finite_row_mean <- function(x) {
  y <- rowMeans(as.matrix(x), na.rm = TRUE)
  y[!is.finite(y)] <- NA_real_
  y
}
# Stable, self-contained Cliff's delta from the Wilcoxon Mann-Whitney statistic.
cliffs_delta <- function(hi, lo) {
  hi <- hi[is.finite(hi)]; lo <- lo[is.finite(lo)]
  if (!length(hi) || !length(lo)) return(NA_real_)
  w <- suppressWarnings(stats::wilcox.test(hi, lo, exact = FALSE))$statistic
  as.numeric(2 * w / (length(hi) * length(lo)) - 1)
}
read_cbio_long <- function(fp, value_name, genes = GENES_1Q21) {
  check_required(fp)
  dt <- data.table::fread(fp, na.strings = c("NA","NaN","NP","N/A",""), check.names = FALSE)
  if (!all(c("STUDY_ID","SAMPLE_ID") %in% names(dt)))
    stop("Expected STUDY_ID and SAMPLE_ID in ", basename(fp))
  absent <- setdiff(genes, names(dt))
  if (length(absent)) stop(basename(fp), " missing 1q21 genes: ", paste(absent, collapse=", "))
  dt <- dt[is_primary(SAMPLE_ID), c("STUDY_ID","SAMPLE_ID",genes), with = FALSE]
  data.table::setnames(dt, c("STUDY_ID","SAMPLE_ID"), c("study","sample"))
  if (anyDuplicated(dt[, .(study, sample)]))
    stop("Duplicate study + sample keys in ", basename(fp), "; resolve upstream")
  out <- data.table::melt(dt, id.vars = c("study","sample"), measure.vars = genes,
                          variable.name = "gene", value.name = value_name,
                          variable.factor = FALSE)
  out[, (value_name) := suppressWarnings(as.numeric(get(value_name)))]
  out
}
read_cbio_wide <- function(fp, prefix, genes = GENES_1Q21) {
  check_required(fp)
  dt <- data.table::fread(fp, na.strings = c("NA","NaN","NP","N/A",""), check.names = FALSE)
  if (!all(c("STUDY_ID","SAMPLE_ID") %in% names(dt))) stop("Missing cBioPortal keys in ", fp)
  absent <- setdiff(genes, names(dt))
  if (length(absent)) stop(basename(fp), " missing genes: ", paste(absent, collapse=", "))
  dt <- dt[is_primary(SAMPLE_ID), c("STUDY_ID","SAMPLE_ID",genes), with = FALSE]
  data.table::setnames(dt, c("STUDY_ID","SAMPLE_ID"), c("study","sample"))
  if (anyDuplicated(dt[, .(study, sample)])) stop("Duplicate study/sample entries in ", fp)
  dt[, (genes) := lapply(.SD, function(x) suppressWarnings(as.numeric(x))), .SDcols = genes]
  data.table::setnames(dt, genes, paste0(prefix, genes))
  dt
}
