# Paths, gene sets, thresholds and shared helpers for the VPS72 genomic analyses.
# Set PROJ_DIR; every path below is built from it. Nothing else needs editing.
suppressPackageStartupMessages({library(data.table);library(broom)})

PROJ_DIR <- "."
OUT_DIR <- file.path(PROJ_DIR,"output")
dir.create(OUT_DIR,recursive=TRUE,showWarnings=FALSE)

# Open-access TCGA/cBioPortal releases; not distributed with this repository.
inputs <- list(
  expr    = file.path(PROJ_DIR,"data/TCGA/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv"),
  mc3     = file.path(PROJ_DIR,"data/TCGA/mc3.v0.2.8.PUBLIC.toil.xena.gz"),
  clin    = file.path(PROJ_DIR,"data/TCGA/combined_study_clinical_data.tsv"),
  gistic  = file.path(PROJ_DIR,"data/TCGA/GDC-PANCAN.gistic.tsv"),
  annot   = file.path(PROJ_DIR,"data/TCGA/merged_sample_quality_annotations.tsv"),
  drivers = file.path(PROJ_DIR,"data/TCGA/cancer_driver_genes.csv"),
  cohorts = file.path(PROJ_DIR,"data/cbioportal"),
  master  = file.path(OUT_DIR,"master_sample_table.tsv"),
  mutmat  = file.path(OUT_DIR,"mutation_matrix.tsv"),
  expr_mat= file.path(OUT_DIR,"complex_expression_medianZ.tsv"),
  cnv_mat = file.path(OUT_DIR,"complex_cnv_mean.tsv"),
  amp_mat = file.path(OUT_DIR,"complex_ampfreq.tsv"),
  assoc   = file.path(OUT_DIR,"mutation_association.tsv")
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

tcga12 <- function(x) substr(as.character(x),1,12)   # patient
tcga16 <- function(x) substr(as.character(x),1,16)   # aliquot

# Clinical and annotation exports rename columns between releases, so fields are
# located by exact name first and then by pattern.
find_col <- function(d,options) {
  z <- options[options %in% names(d)]
  if (length(z)) return(z[1L])
  z <- names(d)[grepl(paste(options,collapse="|"),names(d),ignore.case=TRUE)]
  if (!length(z)) stop("Required field missing: ",paste(options,collapse=", "),
                       "\nAvailable: ",paste(names(d),collapse=", "))
  z[1L]
}

# SRCAP / TIP60(EP400) / INO80 members. VPS72 is in the panel but is kept out of
# both complex scores below, since those scores are used to predict it.
panel_genes <- c("SRCAP","KAT5","INO80","VPS72","DMAP1","RUVBL1","RUVBL2","ACTL6A",
                 "ACTL6B","YEATS4","EP400","TRRAP","ING3","MRGBP","ACTR5")
srcap_genes <- c("SRCAP","DMAP1","EP400","ACTL6A","ACTL6B","RUVBL1","RUVBL2","MRGBP","YEATS4")
tip60_genes <- c("KAT5","TRRAP","ING3","EP400","DMAP1","RUVBL1","RUVBL2","ACTL6A","ACTL6B","MRGBP","YEATS4")
if ("VPS72" %in% c(srcap_genes,tip60_genes)) stop("VPS72 must not appear in its own predictor scores")

common_drivers <- c("TP53","PIK3CA","KRAS","APC","ARID1A","PTEN","CTNNB1","RB1","SMAD4",
                    "ATM","NF1","EGFR","IDH1","NFE2L2","CDKN2A","CDH1","NOTCH1","TERT",
                    "BRAF","FGFR3","KEAP1","STK11","NRAS","MTOR","GNAS","FBXW7","AKT1",
                    "MAP2K1","MLH1","VHL")

cohorts <- c(BRCA="brca_tcga_pan_can_atlas_2018", GBM="gbm_tcga_pan_can_atlas_2018",
             HNSC="hnsc_tcga_pan_can_atlas_2018", LUAD="luad_tcga_pan_can_atlas_2018",
             LUSC="lusc_tcga_pan_can_atlas_2018", SKCM="skcm_tcga_pan_can_atlas_2018",
             CRC="coadread_tcga_pan_can_atlas_2018")

nonsilent <- c("Missense_Mutation","Nonsense_Mutation","Splice_Site","Frame_Shift_Del",
               "Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Translation_Start_Site",
               "Nonstop_Mutation")

fdr_cut       <- 0.05
min_mutated   <- 10   # carriers before a gene or variant is modelled
min_auc_n     <- 5    # mutated samples before an AUC is computed
amp_gistic    <- 2    # GISTIC2 call counted as amplification
high_quantile <- 0.75 # top quartile defines VPS72-high in the co-mutation test
enrich_fdr    <- 0.10 # looser cut feeding over-representation testing

# Non-silent calls, one row per patient and gene. Used by 01 and 07.
read_mc3 <- function() {
  m <- fread(require_input(inputs$mc3))
  setnames(m,c(find_col(m,c("sample","Tumor_Sample_Barcode")),
               find_col(m,c("gene","Hugo_Symbol")),
               find_col(m,c("effect","Variant_Classification"))),
           c("sample","gene","effect"))
  m <- m[toupper(effect) %chin% toupper(nonsilent)]
  if (!nrow(m)) stop("No non-silent calls after filtering; check the effect column vocabulary")
  m[,`:=`(patient=tcga12(sample),gene=toupper(gene))]
  m
}

# Cancer type per aliquot, from the merged sample quality annotation. Used by 02 and 04.
read_annot <- function() {
  a <- fread(require_input(inputs$annot))
  a <- unique(data.table(sample16=tcga16(a[[find_col(a,c("aliquot_barcode","sample"))]]),
                         cancer=as.character(a[[find_col(a,c("cancer type","cancer_type"))]])))
  a[!is.na(cancer) & nzchar(cancer)]
}
