# VPS72 pan-cancer expression/CNV: data-analysis configuration.
# Run from this directory, or set VPS72_DATA_DIR and VPS72_OUT_DIR.
# Files and values below identify DATASETS, not any investigator's local directories.
suppressPackageStartupMessages({library(data.table); library(matrixStats)})

input_dir <- Sys.getenv("VPS72_DATA_DIR", unset=file.path(getwd(), "data"))
output_dir <- Sys.getenv("VPS72_OUT_DIR", unset=file.path(getwd(), "results"))
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

inputs <- list(
  tcga_expression = file.path(input_dir, "expression.xena.gz"),
  tcga_gistic = file.path(input_dir, "TCGA.PANCAN.sampleMap_Gistic2_CopyNumber_Gistic2_all_data_by_genes (1).gz"),
  tcga_clinical = file.path(input_dir, "combined_study_clinical_data (8).tsv"),
  tpm_expression = file.path(input_dir, "TcgaTargetGtex_rsem_gene_tpm.gz"),
  tpm_phenotype = file.path(input_dir, "TcgaTargetGTEX_phenotype.txt.gz"),
  copy_number = file.path(input_dir, "amplificaion.txt"),
  copy_expression = file.path(input_dir, "cna.txt"),
  copy_metadata = file.path(input_dir, "full dataset.tsv"),
  xena_merged = file.path(input_dir, "xena_merge.tsv"),
  xena_signatures = file.path(input_dir, "Scores_160_Signatures.tsv")
)

# Exact seven gene sets used in the original expression/CNV signature workflow.
signature_files <- c(
  MYC_V1="HALLMARK_MYC_TARGETS_V1.v2025.1.Hs.tsv",
  MYC_V2="HALLMARK_MYC_TARGETS_V2.v2025.1.Hs.tsv",
  P53="KEGG_P53_SIGNALING_PATHWAY.v2025.1.Hs.tsv",
  CCYCLE="REACTOME_CELL_CYCLE.v2025.1.Hs.tsv",
  PROLIF="CELL_PROLIFERATION_GO_0008283.v2025.1.Hs.tsv",
  E2F="HALLMARK_E2F_TARGETS.v2025.1.Hs.tsv",
  G2M="HALLMARK_G2M_CHECKPOINT.v2025.1.Hs.tsv"
)
min_group_n <- 10L
min_tumor_normal_n <- 3L # sensitivity/coverage rule for the broad comparison
min_stage_n <- 10L
srcap_genes <- c("SRCAP","DMAP1","EP400","ACTL6A","ACTL6B","RUVBL1","RUVBL2","MRGBP","YEATS4")
tip60_genes <- c("KAT5","TRRAP","ING3","EP400","DMAP1","RUVBL1","RUVBL2","ACTL6A","ACTL6B","MRGBP","YEATS4")
ino80_genes <- c("INO80","ACTR8","ACTL6A","CHRAC1","DMAP1","BAZ1A","BAZ1B","ANP32E")
# These displayed gene modules are intentionally distinct from the VPS72-excluded
# SRCAP/TIP60 modules used to regress VPS72 expression.
comparison_modules <- list(
  SRCAP=c("SRCAP","VPS72","RUVBL1","RUVBL2","YEATS4","TRRAP"),
  TIP60=c("KAT5","EP400","EPC1","EPC2","ING3","YEATS4","H2AFZ","H2AFV"),
  INO80=ino80_genes
)

require_input <- function(fp) {
  if (!file.exists(fp)) stop("Missing input: ", fp, call.=FALSE)
  fp
}
tcga15 <- function(x) {
  x <- toupper(gsub("\\.", "-", trimws(as.character(x))))
  ok <- grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}", x)
  ifelse(ok, substr(x, 1L, 15L), NA_character_)
}
tcga12 <- function(x) {
  x <- toupper(gsub("\\.", "-", trimws(as.character(x))))
  ifelse(grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", x), substr(x,1L,12L), NA_character_)
}
# Use proper TCGA sample type; a 12-character barcode cannot distinguish tissue states.
keep_primary <- function(x) !is.na(x) & substr(x,14L,15L)=="01"
keep_normal <- function(x) !is.na(x) & substr(x,14L,15L)=="11"

# Allow SYMBOL, SYMBOL|ENTREZ, or Ensembl identifiers. Map Ensembl IDs only
# when required; duplicate symbols choose the row with highest mean expression.
resolve_symbols <- function(ids) {
  ids <- toupper(trimws(as.character(ids)))
  sym <- sub("\\|.*$", "", ids)
  ens <- grepl("^ENSG[0-9]+(\\.[0-9]+)?$", sym)
  if (any(ens)) {
    if (!requireNamespace("AnnotationDbi", quietly=TRUE) ||
        !requireNamespace("org.Hs.eg.db", quietly=TRUE))
      stop("Ensembl IDs detected; install AnnotationDbi and org.Hs.eg.db")
    keys <- unique(sub("\\.[0-9]+$", "", sym[ens]))
    map <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys=keys,
                                 columns="SYMBOL", keytype="ENSEMBL")
    map <- map[!is.na(map$SYMBOL) & nzchar(map$SYMBOL), , drop=FALSE]
    map <- map[order(map$ENSEMBL,map$SYMBOL),,drop=FALSE]
    map <- map[!duplicated(map$ENSEMBL),,drop=FALSE]
    lookup <- setNames(toupper(map$SYMBOL),map$ENSEMBL)
    sym[ens] <- unname(lookup[sub("\\.[0-9]+$", "",sym[ens])])
  }
  sym
}
load_gene_matrix <- function(fp, required_genes, label) {
  dt <- data.table::fread(require_input(fp), check.names=FALSE)
  if (ncol(dt)<2) stop(label, ": expected gene ID plus sample columns")
  symbols <- resolve_symbols(dt[[1L]])
  idx <- which(!is.na(symbols) & symbols %in% toupper(required_genes))
  if (!length(idx)) stop(label, ": none of the requested gene symbols found; verify gene ID format")
  raw <- as.matrix(dt[idx, -1L, with=FALSE]); storage.mode(raw) <- "double"
  rownames(raw) <- symbols[idx]
  ranks <- rowMeans(raw, na.rm=TRUE); ranks[!is.finite(ranks)] <- -Inf
  ord <- order(rownames(raw), -ranks, seq_along(ranks))
  raw <- raw[ord,,drop=FALSE]
  raw <- raw[!duplicated(rownames(raw)),,drop=FALSE]
  colnames(raw) <- tcga15(colnames(raw))
  raw
}
# Collapse aliquots for the SAME TCGA sample/tissue state, never across 01/11.
collapse_samples <- function(mat) {
  good <- !is.na(colnames(mat))
  mat <- mat[,good,drop=FALSE]
  if (!anyDuplicated(colnames(mat))) return(mat)
  labs <- unique(colnames(mat))
  out <- vapply(labs, function(x) {
    m <- mat[,colnames(mat)==x,drop=FALSE]
    if (ncol(m)==1L) as.numeric(m[,1L]) else rowMeans(m,na.rm=TRUE)
  }, numeric(nrow(mat)))
  if (is.null(dim(out))) out <- matrix(out,nrow=nrow(mat))
  rownames(out) <- rownames(mat); colnames(out) <- labs
  out[!is.finite(out)] <- NA_real_
  out
}
write_tsv <- function(dat,name) {
  data.table::fwrite(as.data.table(dat),file.path(output_dir,name),sep="\t",na="NA")
}
safe_wilcox <- function(x, group, min_n=min_group_n) {
  ok <- is.finite(x) & !is.na(group)
  x <- x[ok]; group <- as.character(group[ok])
  if (sum(group=="High") < min_n || sum(group=="Low") < min_n) return(NA_real_)
  tryCatch(suppressWarnings(wilcox.test(x[group=="High"],x[group=="Low"],exact=FALSE)$p.value),
           error=function(e) NA_real_)
}
safe_spearman <- function(x,y,min_n=10L) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok)<min_n || length(unique(x[ok]))<2L || length(unique(y[ok]))<2L)
    return(c(n=sum(ok),rho=NA_real_,p=NA_real_))
  test <- suppressWarnings(cor.test(x[ok],y[ok],method="spearman",exact=FALSE))
  c(n=sum(ok),rho=unname(test$estimate),p=test$p.value)
}
