# Paths, constants and shared helpers for the VPS72 TCGA ATAC-seq analyses.
# Set PROJ_DIR; every path below is built from it. Nothing else needs editing.
suppressPackageStartupMessages({
  library(data.table);library(matrixStats);library(limma)
  library(GenomicRanges);library(IRanges)
})

PROJ_DIR <- "."
OUT_DIR <- file.path(PROJ_DIR,"output")
CACHE_DIR <- file.path(OUT_DIR,"cache")
for (d in c(OUT_DIR,CACHE_DIR)) dir.create(d,recursive=TRUE,showWarnings=FALSE)

# TCGA pan-cancer ATAC-seq (Corces et al., Science 2018) and UCSC Xena TOIL
# expression. Not distributed with this repository.
inputs <- list(
  atac     = file.path(PROJ_DIR,"data/ATAC/TCGA_ATAC_peak_Log2Counts_dedup_sample.gz"),
  idmap    = file.path(PROJ_DIR,"data/ATAC/TCGA_identifier_mapping.txt"),
  probemap = file.path(PROJ_DIR,"data/ATAC/TCGA_ATAC_peak.all.probeMap"),
  rna      = file.path(PROJ_DIR,"data/TCGA/TcgaTargetGtex_rsem_gene_tpm.gz"),
  genesets = file.path(PROJ_DIR,"data/gene_sets"),
  coords   = file.path(OUT_DIR,"atac_matrix_with_coords.tsv"),
  groups   = file.path(OUT_DIR,"sample_groups.tsv"),
  da       = file.path(OUT_DIR,"differential_accessibility.tsv")
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

find_col <- function(d,options) {
  z <- options[options %in% names(d)]
  if (length(z)) return(z[1L])
  z <- names(d)[grepl(paste(options,collapse="|"),names(d),ignore.case=TRUE)]
  if (!length(z)) stop("Required field missing: ",paste(options,collapse=", "),
                       "\nAvailable: ",paste(names(d),collapse=", "))
  z[1L]
}

# Ensembl coordinates are cached on first use, so results do not shift with
# Ensembl releases and later runs work offline.
ensembl_cached <- function(name,attributes,filters,values) {
  fp <- file.path(CACHE_DIR,paste0(name,".tsv"))
  if (file.exists(fp)) return(fread(fp))
  mart <- biomaRt::useEnsembl("genes",dataset="hsapiens_gene_ensembl")
  res <- biomaRt::getBM(attributes,filters=filters,values=values,mart=mart)
  if (!nrow(res)) stop("Ensembl returned nothing for ",name)
  fwrite(res,fp,sep="\t")
  as.data.table(res)
}

vps72_ensg <- "ENSG00000163159"
fdr_cut    <- 0.05
lfc_annot  <- 0.25   # effect-size floor for peak annotation and TSS profiling
tss_window <- 5000
tss_bin    <- 250
flank_bp   <- 50000  # flanking sequence added to TF gene bodies
min_locus_n <- 3     # samples required before a locus mean is reported

# Include cancer type in the differential model.
#
# FALSE is the published analysis: the pan-cancer contrast reported in the
# manuscript. Setting it TRUE adds a lineage term, which absorbs baseline
# differences in accessibility between tumour types and makes the VPS72
# coefficient a within-lineage contrast. That is a useful sensitivity check
# for anyone adapting this pipeline, but it does not reproduce the published
# estimates. Either way the script reports which model it fitted.
adjust_for_cancer_type <- FALSE

tf_genes <- c("MYC",paste0("E2F",1:8),"JUN","FOS","HIF1A","NFE2L2","STAT3","RELA",
              "TP53",paste0("TEAD",1:4),"FOXM1")

# Loci for the per-cancer accessibility means (hg38).
loci <- list(MYC =list(chrom="chr8", start=127000000,end=128800000),
             E2F1=list(chrom="chr20",start=33600000, end=34000000),
             E2F3=list(chrom="chr6", start=20400000, end=21000000))

geneset_files <- c(
  CELL_PROLIFERATION_GO   = "CELL_PROLIFERATION_GO_0008283.v2025.1.Hs.tsv",
  HALLMARK_E2F_TARGETS    = "HALLMARK_E2F_TARGETS.v2025.1.Hs.tsv",
  HALLMARK_G2M_CHECKPOINT = "HALLMARK_G2M_CHECKPOINT.v2025.1.Hs.tsv",
  HALLMARK_MYC_TARGETS_V1 = "HALLMARK_MYC_TARGETS_V1.v2025.1.Hs.tsv",
  HALLMARK_MYC_TARGETS_V2 = "HALLMARK_MYC_TARGETS_V2.v2025.1.Hs.tsv",
  REACTOME_CELL_CYCLE     = "REACTOME_CELL_CYCLE.v2025.1.Hs.tsv",
  KEGG_P53_SIGNALING      = "KEGG_P53_SIGNALING_PATHWAY.v2025.1.Hs.tsv")
