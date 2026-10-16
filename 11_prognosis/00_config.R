# Run scripts from the 11_prognosis directory.
PROJ_DIR <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
CBIO_DIR <- Sys.getenv("VPS72_PROGNOSIS_CBIO", unset = file.path(PROJ_DIR, "data", "cbioportal"))
TCGA_DIR <- Sys.getenv("VPS72_PROGNOSIS_TCGA", unset = file.path(PROJ_DIR, "data", "TCGA"))
OUT_DIR <- Sys.getenv("VPS72_PROGNOSIS_OUTPUT", unset = file.path(PROJ_DIR, "results"))
F_EXPR <- file.path(CBIO_DIR, "mRNA expression z-scores relative to all samples (log RNA Seq V2 RSEM).txt")
F_CNV <- file.path(CBIO_DIR, "Log2 copy-number values.txt")
F_CLINICAL <- file.path(TCGA_DIR, "combined_study_clinical_data (4).tsv")
F_PROLIF <- file.path(TCGA_DIR, "Scores_160_Signatures.tsv") # Optional, for model M8.
NEIGHBORS <- c("CHD1L", "PRKAB2", "FMO5", "BCL9", "H2AC18", "ANP32E", "MCL1", "ECM1",
  "ENSA", "ARNT", "SETDB1", "CERS2", "MLLT11", "PRUNE1", "PIP5K1A", "PSMD4",
  "ZNF687", "S100A4", "ILF2", "GATAD2B")
CN_GAIN_CUT <- 0.20
CN_LOSS_CUT <- -0.20
MIN_MODEL_N <- 50L
MIN_MODEL_EVENTS <- 20L
suppressPackageStartupMessages({ library(data.table); library(survival) })
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
write_result <- function(x, name) {
  path <- file.path(OUT_DIR, name)
  data.table::fwrite(x, path, sep = "\t", na = "NA")
  message("Wrote: ", path)
}
