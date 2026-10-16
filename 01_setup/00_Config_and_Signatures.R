# Defines all project paths, output directories, analysis parameters, and loads MSigDB gene signatures for the VPS72 pan-cancer dotblot pipeline.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(matrixStats)
  library(dplyr)
})

# ---- Base path (update only this) ----
base_dir <- "path/to/ATAC_2"
data_dir <- file.path(base_dir, "Data")

# ---- Signature files ----
sig_files <- c(
  MYC_V1 = file.path(data_dir, "HALLMARK_MYC_TARGETS_V1.v2025.1.Hs.tsv"),
  MYC_V2 = file.path(data_dir, "HALLMARK_MYC_TARGETS_V2.v2025.1.Hs.tsv"),
  P53    = file.path(data_dir, "KEGG_P53_SIGNALING_PATHWAY.v2025.1.Hs.tsv"),
  CCYCLE = file.path(data_dir, "REACTOME_CELL_CYCLE.v2025.1.Hs.tsv"),
  PROLIF = file.path(data_dir, "CELL_PROLIFERATION_GO_0008283.v2025.1.Hs.tsv"),
  E2F    = file.path(data_dir, "HALLMARK_E2F_TARGETS.v2025.1.Hs.tsv"),
  G2M    = file.path(data_dir, "HALLMARK_G2M_CHECKPOINT.v2025.1.Hs.tsv")
)

# ---- Data files ----
expr_fp <- file.path(data_dir, "expression.xena.gz")
cnv_fp  <- file.path(data_dir, "TCGA.PANCAN.sampleMap_Gistic2_CopyNumber_Gistic2_all_data_by_genes (1).gz")
clin_fp <- file.path(data_dir, "combined_study_clinical_data (8).tsv")

# ---- Output directories ----
out_root       <- file.path(base_dir, "derived", "VPS72_signature_dotblot")
dir_fig_expr   <- file.path(out_root, "figures_expr")
dir_fig_cnv    <- file.path(out_root, "figures_cnv")
dir_fig_violin <- file.path(out_root, "figures_violin")
dir_tables     <- file.path(out_root, "tables")

for (d in c(dir_fig_expr, dir_fig_cnv, dir_fig_violin, dir_tables))
  dir.create(d, recursive=TRUE, showWarnings=FALSE)

# ---- Analysis parameters ----
vps72_gene          <- "VPS72"
score_method_expr   <- "zmean"   # z-score then mean
score_method_cnv    <- "mean"
fraction_threshold  <- 0.0
min_group_n         <- 10
dotplot_point_min   <- 1.5
dotplot_point_max   <- 8.0
col_low             <- "#2B6CB0"
col_mid             <- "white"
col_high            <- "#8B0000"

# ---- Load MSigDB signatures (GENE_SYMBOLS row, comma-separated) ----
load_msigdb_tsv <- function(fp) {
  if (!file.exists(fp)) { warning("Missing: ", fp); return(character(0)) }
  dt <- fread(fp, sep="\t", header=FALSE)
  gs <- dt[V1=="GENE_SYMBOLS", V2]
  if (length(gs)==0) return(character(0))
  unique(toupper(trimws(unlist(strsplit(gs, ",")))))
}

sig_list  <- lapply(sig_files, load_msigdb_tsv)
sig_names <- names(sig_list)

# Save gene lists
sig_long <- rbindlist(lapply(sig_names, function(nm)
  data.table(signature=nm, gene=sig_list[[nm]])))
fwrite(sig_long, file.path(dir_tables,"signature_genes_long.tsv"), sep="\t")
saveRDS(sig_list, file.path(dir_tables,"signature_genes_list.rds"))

message("Config loaded. Signatures: ", paste(sig_names, collapse=", "))
