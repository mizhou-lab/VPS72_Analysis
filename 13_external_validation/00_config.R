# Paths, contrasts, gene sets and shared helpers for the external validation
# of VPS72/YL1 in the Jostes et al. perturbation data.
# Set PROJ_DIR; every path below is built from it.
suppressPackageStartupMessages({library(data.table)})

PROJ_DIR <- "."
OUT_DIR <- file.path(PROJ_DIR,"output")
dir.create(OUT_DIR,recursive=TRUE,showWarnings=FALSE)

# Jostes et al. 2024, Genes & Development. SK-MEL-147 melanoma cells.
#   GSE242227  RNA-seq, shRNA knockdowns vs scrambled control
#   GSE246121  H4ac ChIP-seq, control vs YL1 knockdown
gse_rnaseq <- "GSE242227"
gse_chip   <- "GSE246121"

inputs <- list(
  geo_dir  = file.path(PROJ_DIR,"data/GEO"),
  genesets = file.path(PROJ_DIR,"data/gene_sets"),
  deg      = file.path(OUT_DIR,"deg_standardized.tsv"),
  sets     = file.path(OUT_DIR,"regulated_gene_sets.tsv"),
  bigwig   = file.path(OUT_DIR,"h4ac_bigwig_samples.tsv"),
  fgsea    = file.path(OUT_DIR,"fgsea_hallmark.tsv")
)

require_input <- function(p) {
  if (!file.exists(p)) stop("Cannot see required input: ",p,
                            "\nDownload it and place it under data/ (see README), or run the script that writes it.")
  p
}

write_tsv <- function(x,file) {
  fp <- file.path(OUT_DIR,file)
  fwrite(x,fp,sep="\t")
  cat("wrote",fp,"(",nrow(x),"rows )\n")
  invisible(fp)
}

find_col <- function(d,options) {
  z <- options[options %in% names(d)]
  if (length(z)) return(z[1L])
  z <- names(d)[grepl(paste(options,collapse="|"),names(d),ignore.case=TRUE)]
  if (!length(z)) stop("Required field missing: ",paste(options,collapse=", "),
                       "\nAvailable: ",paste(names(d),collapse=", "))
  z[1L]
}

# Author-provided DESeq2 result tables, day 6 contrasts. These are Jostes et
# al.'s own differential testing, not a re-analysis from raw counts.
deg_files <- c(VPS72="GSE242227_shYL1_vs_SCR_d6.csv",
               SRCAP="GSE242227_shSRCAP_vs_SCR_d6.csv",
               H2AZ1="GSE242227_shH2AZ1_vs_SCR_d6.csv")

# H4ac ChIP-seq conditions: two independent controls, two independent YL1
# hairpins. Matched on the file-name tokens GEO ships.
chip_conditions <- c(shSCR="scr-h4ac", shSCR2="scr_h4ac",
                     YL1sh30="sh30", YL1sh84="sh84")
chip_controls <- c("shSCR","shSCR2")
chip_kd       <- c("YL1sh30","YL1sh84")

signature_files <- c(
  KEGG_P53_SIGNALING_PATHWAY    = "KEGG_P53_SIGNALING_PATHWAY.v2025.1.Hs.tsv",
  REACTOME_CELL_CYCLE           = "REACTOME_CELL_CYCLE.v2025.1.Hs.tsv",
  CELL_PROLIFERATION_GO_0008283 = "CELL_PROLIFERATION_GO_0008283.v2025.1.Hs.tsv",
  HALLMARK_E2F_TARGETS          = "HALLMARK_E2F_TARGETS.v2025.1.Hs.tsv",
  HALLMARK_G2M_CHECKPOINT       = "HALLMARK_G2M_CHECKPOINT.v2025.1.Hs.tsv",
  HALLMARK_MYC_TARGETS_V1       = "HALLMARK_MYC_TARGETS_V1.v2025.1.Hs.tsv",
  HALLMARK_MYC_TARGETS_V2       = "HALLMARK_MYC_TARGETS_V2.v2025.1.Hs.tsv")

padj_cut  <- 0.05
lfc_cut   <- 0.5
tss_flank <- 500    # half-width of the promoter window for the H4ac ranking
min_set_n <- 10     # genes required before a partition is tested for enrichment
chip_genome <- "hg19"   # the ChIP data are hg19; the pan-cancer ATAC work is hg38

# MSigDB TSV exports carry members on a GENE_SYMBOLS row.
read_geneset <- function(fp) {
  raw <- readLines(require_input(fp),warn=FALSE)
  row <- grep("GENE_SYMBOLS",raw,value=TRUE)[1L]
  if (is.na(row)) stop("No GENE_SYMBOLS row in ",fp)
  g <- toupper(trimws(strsplit(sub("^[^\t]*\t","",row),",")[[1L]]))
  g[nzchar(g)]
}

load_signatures <- function() {
  s <- lapply(names(signature_files),
              function(nm) read_geneset(file.path(inputs$genesets,signature_files[[nm]])))
  names(s) <- names(signature_files)
  s
}

# Seven mutually exclusive categories over three perturbations: which programmes
# VPS72 shares with SRCAP, with H2A.Z.1, with both, or with neither.
make_partitions <- function(v,s,h) {
  list(VPS72_only  = setdiff(setdiff(v,s),h),
       SRCAP_only  = setdiff(setdiff(s,v),h),
       H2AZ1_only  = setdiff(setdiff(h,v),s),
       VPS72_SRCAP = setdiff(intersect(v,s),h),
       VPS72_H2AZ1 = setdiff(intersect(v,h),s),
       SRCAP_H2AZ1 = setdiff(intersect(s,h),v),
       CORE        = Reduce(intersect,list(v,s,h)))
}
