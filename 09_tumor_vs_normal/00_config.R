# Paths, cohorts, gene sets and shared helpers for the VPS72 tumour-vs-normal
# analyses. Set PROJ_DIR; every path below is built from it.
suppressPackageStartupMessages({library(data.table)})

PROJ_DIR <- "."
OUT_DIR <- file.path(PROJ_DIR,"output")
dir.create(OUT_DIR,recursive=TRUE,showWarnings=FALSE)

# GEO series are downloaded by 01. TCGA/Xena and MSigDB files are not
# distributed with this repository; see the README.
inputs <- list(
  geo_dir   = file.path(PROJ_DIR,"data/GEO"),
  xena_expr = file.path(PROJ_DIR,"data/TCGA/TcgaTargetGtex_rsem_gene_tpm.gz"),
  xena_pheno= file.path(PROJ_DIR,"data/TCGA/TcgaTargetGTEX_phenotype.txt.gz"),
  genesets  = file.path(PROJ_DIR,"data/gene_sets"),
  geo_expr  = file.path(OUT_DIR,"geo_normal_expression.tsv.gz"),
  geo_meta  = file.path(OUT_DIR,"geo_normal_metadata.tsv"),
  tcga      = file.path(OUT_DIR,"tcga_tumor_normal.tsv"),
  master    = file.path(OUT_DIR,"master_sample_table.tsv"),
  scores    = file.path(OUT_DIR,"pathway_scores_per_sample.tsv.gz"),
  highlow   = file.path(OUT_DIR,"pathway_high_low_results.tsv"),
  spearman  = file.path(OUT_DIR,"per_cancer_spearman.tsv"),
  meta_state= file.path(OUT_DIR,"meta_by_state.tsv"),
  interact  = file.path(OUT_DIR,"per_cancer_interactions.tsv")
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

tcga12 <- function(x) substr(as.character(x),1,12)
tcga15 <- function(x) substr(as.character(x),1,15)

find_col <- function(d,options) {
  z <- options[options %in% names(d)]
  if (length(z)) return(z[1L])
  z <- names(d)[grepl(paste(options,collapse="|"),names(d),ignore.case=TRUE)]
  if (!length(z)) stop("Required field missing: ",paste(options,collapse=", "),
                       "\nAvailable: ",paste(names(d),collapse=", "))
  z[1L]
}

# Normal-tissue bulk RNA-seq references. Within each series only the untreated
# normal samples are kept; 01 records which were dropped and why.
geo_series <- c(GSE190615="28 normal primary human cell types",
                GSE117271="normal prostate, FACS-sorted populations",
                GSE161892="normal breast epithelial populations",
                GSE227502="primary bronchial epithelial cells, normoxia arm only")

# TCGA projects with enough adjacent-normal samples to support a paired contrast.
eligible <- c("BLCA","BRCA","COAD","ESCA","HNSC","KICH","KIRC","KIRP",
              "LIHC","LUAD","LUSC","PRAD","READ","STAD","THCA","UCEC")

tissue_map <- data.table(
  cancer=eligible,
  tissue=c("Bladder","Breast","Colon","Esophagus","Head and neck","Kidney",
           "Kidney","Kidney","Liver","Lung","Lung","Prostate","Rectum",
           "Stomach","Thyroid","Uterus"))

# MSigDB Hallmark sets scored per sample. The composite proliferation score is
# the mean of prolif_sets; the remaining sets are reported individually.
hallmark_files <- c(
  HALLMARK_E2F_TARGETS                = "HALLMARK_E2F_TARGETS.v2025.1.Hs.tsv",
  HALLMARK_G2M_CHECKPOINT             = "HALLMARK_G2M_CHECKPOINT.v2025.1.Hs.tsv",
  HALLMARK_MYC_TARGETS_V1             = "HALLMARK_MYC_TARGETS_V1.v2025.1.Hs.tsv",
  HALLMARK_MYC_TARGETS_V2             = "HALLMARK_MYC_TARGETS_V2.v2025.1.Hs.tsv",
  HALLMARK_MITOTIC_SPINDLE            = "HALLMARK_MITOTIC_SPINDLE.v2025.1.Hs.tsv",
  HALLMARK_MTORC1_SIGNALING           = "HALLMARK_MTORC1_SIGNALING.v2025.1.Hs.tsv",
  HALLMARK_DNA_REPAIR                 = "HALLMARK_DNA_REPAIR.v2025.1.Hs.tsv",
  HALLMARK_OXIDATIVE_PHOSPHORYLATION  = "HALLMARK_OXIDATIVE_PHOSPHORYLATION.v2025.1.Hs.tsv")

prolif_sets <- c("HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT",
                 "HALLMARK_MYC_TARGETS_V1","HALLMARK_MYC_TARGETS_V2",
                 "HALLMARK_MITOTIC_SPINDLE")

program_labels <- c(
  proliferation_score                 = "Composite proliferation score",
  HALLMARK_E2F_TARGETS                = "E2F targets",
  HALLMARK_G2M_CHECKPOINT             = "G2M checkpoint",
  HALLMARK_MYC_TARGETS_V1             = "MYC targets V1",
  HALLMARK_MYC_TARGETS_V2             = "MYC targets V2",
  HALLMARK_MITOTIC_SPINDLE            = "Mitotic spindle",
  HALLMARK_MTORC1_SIGNALING           = "mTORC1 signaling",
  HALLMARK_DNA_REPAIR                 = "DNA repair",
  HALLMARK_OXIDATIVE_PHOSPHORYLATION  = "Oxidative phosphorylation")

pretty_name <- function(x) {
  lab <- unname(program_labels[x])
  fifelse(is.na(lab),gsub("_"," ",sub("^HALLMARK_","",x)),lab)
}

vps72_ensg <- "ENSG00000163159"
fdr_cut    <- 0.05
min_cor_n  <- 6    # samples in a cancer/state cell before a correlation is run
min_state_n<- 5    # samples per state before an interaction model is fitted
min_pairs  <- 10   # paired patients before the mixed-model sensitivity is fitted
min_meta_k <- 3    # cancer types before a random-effects meta-analysis is run

# Ranks are standardized within tissue state, so a within-state slope
# corresponds closely to Spearman rho and the interaction term is a difference
# of correlations rather than of raw scales.
rank_z <- function(x) {
  r <- rank(x,ties.method="average",na.last="keep")
  if (sum(!is.na(r))<3L || sd(r,na.rm=TRUE)==0) return(rep(NA_real_,length(x)))
  as.numeric(scale(r))
}

# MSigDB TSV exports carry members on a GENE_SYMBOLS row.
read_geneset <- function(fp) {
  raw <- readLines(require_input(fp),warn=FALSE)
  row <- grep("GENE_SYMBOLS",raw,value=TRUE)[1L]
  if (is.na(row)) stop("No GENE_SYMBOLS row in ",fp)
  toupper(trimws(strsplit(sub("^[^\t]*\t","",row),",")[[1L]]))
}
