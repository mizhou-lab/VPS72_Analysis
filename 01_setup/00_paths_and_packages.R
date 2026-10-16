# Set project paths and load the packages used by the Figure 2 workflow.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readxl)
  library(openxlsx)
  library(pheatmap)
  library(ggrepel)
  library(ggpubr)
  library(purrr)
  library(tibble)
  library(stringr)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

project_root <- "."
data_dir <- file.path(project_root, "Data", "TCGA")
figure_dir <- file.path(project_root, "Figures", "Fig2_Expression")
supp_dir <- file.path(project_root, "Figures", "Supplementary")
table_dir <- file.path(project_root, "Results", "Tables")
plot_dir <- file.path(project_root, "Results", "Plots")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

expr_file <- file.path(data_dir, "TcgaTargetGtex_rsem_gene_tpm.gz")
pheno_file <- file.path(data_dir, "TcgaTargetGTEX_phenotype.txt.gz")
expr_sym_file <- file.path(data_dir, "TcgaTargetGtex_rsem_gene_tpm_annotated_bySymbol.csv.gz")
master_data_file <- file.path(table_dir, "Figure2_MasterData.xlsx")
