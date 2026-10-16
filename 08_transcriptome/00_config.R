# Run scripts from the 08_transcriptome directory.
PROJ_DIR <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
DATA_DIR <- Sys.getenv("VPS72_TRANSCRIPTOME_INPUT", unset = file.path(PROJ_DIR, "data", "Transcriptome"))
OUT_DIR <- Sys.getenv("VPS72_TRANSCRIPTOME_OUTPUT", unset = file.path(PROJ_DIR, "results"))
CANCERS <- c("BRCA", "CRC", "GBM", "HNSC", "LUAD", "LUSC")
SIGNATURES <- c("VPS72", "SRCAP", "SRCAP_resid")
FDR_CUT <- 0.05
DEG_LFC_CUT <- 0.50  # For descriptive DEG counts; fgsea uses all genes.
FGSEA_PERMUTATIONS <- 10000L
FGSEA_SEED <- 2026L
FGSEA_MIN_SIZE <- 15L
FGSEA_MAX_SIZE <- 500L
# For stable gene sets, set this to a versioned MSigDB Hallmark GMT file.
HALLMARK_GMT <- Sys.getenv("VPS72_HALLMARK_GMT", unset = "")
# Filenames follow the supplied per-cancer DEG exports; edit here if different.
DEG_FILES <- list(
  VPS72 = setNames(file.path(DATA_DIR, "VPS72", CANCERS,
    c("BRCA.tsv", "CRC.tsv", "table (22).tsv", "HNSC.tsv", "LUAD.tsv", "LUSC.tsv")), CANCERS),
  SRCAP = setNames(file.path(DATA_DIR, "SRCAP", CANCERS,
    paste0(c("table (34)", "table (32)", "table (30)", "table (28)", "table (37)", "table (36)"), ".tsv")), CANCERS),
  SRCAP_resid = setNames(file.path(DATA_DIR, "SRCAP_resid", CANCERS,
    paste0(c("table (35)", "table (33)", "table (31)", "table (29)", "table (39)", "table (38)"), ".tsv")), CANCERS)
)
suppressPackageStartupMessages(library(data.table))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
write_result <- function(x, name) {
  path <- file.path(OUT_DIR, name)
  data.table::fwrite(x, path, sep = "\t", na = "NA")
  message("Wrote: ", path)
}
