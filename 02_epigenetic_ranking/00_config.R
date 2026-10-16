# Paths, the epigenetic regulator panel and shared helpers for the pan-cancer
# disruption ranking. Set PROJ_DIR; every path below is built from it.
suppressPackageStartupMessages({library(data.table);library(matrixStats)})

PROJ_DIR <- "."
DATA_DIR <- file.path(PROJ_DIR,"data","TCGA")
OUT_DIR  <- file.path(PROJ_DIR,"output")
CACHE_DIR<- file.path(OUT_DIR,"cache")
for (d in c(OUT_DIR,CACHE_DIR)) dir.create(d,recursive=TRUE,showWarnings=FALSE)

inputs <- list(
  expr   = file.path(DATA_DIR,"EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv"),
  gistic = file.path(DATA_DIR,"GDC-PANCAN.gistic.tsv"),
  annot  = file.path(DATA_DIR,"merged_sample_quality_annotations.tsv"))

require_input <- function(p) {
  if (!file.exists(p)) stop("Cannot see required input: ",p,
                            "\nDownload it and place it under data/TCGA/ (see README).")
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

tcga12 <- function(x) toupper(substr(as.character(x),1,12))
tcga16 <- function(x) toupper(substr(as.character(x),1,16))

# Cancer type comes from the TCGA sample annotation, keyed on barcode. It is
# NOT derived from the tissue source site (barcode positions 6-7): a source
# site can contribute to several projects, so that mapping is ambiguous by
# construction. See the README.
read_sample_annotation <- function() {
  a <- fread(require_input(inputs$annot))
  a <- data.table(patient=tcga12(a[[find_col(a,c("patient_barcode","aliquot_barcode","sample"))]]),
                  cancer=toupper(as.character(a[[find_col(a,c("cancer type","cancer_type"))]])))
  a <- a[nzchar(patient) & nzchar(cancer) & cancer!="NA"]
  a <- a[,.SD[1L],by=patient]
  if (!nrow(a)) stop("No usable rows in the sample annotation file")
  a
}

# Ensembl to symbol, cached on first use so later runs are offline and results
# do not shift with annotation releases.
ensembl_to_symbol <- function(ids) {
  fp <- file.path(CACHE_DIR,"ensembl_to_symbol.tsv")
  if (file.exists(fp)) return(fread(fp))
  suppressPackageStartupMessages({library(org.Hs.eg.db);library(AnnotationDbi)})
  res <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db,keys=unique(ids),
                                                keytype="ENSEMBL",columns="SYMBOL"))
  m <- as.data.table(res)[!is.na(SYMBOL) & SYMBOL!=""]
  setnames(m,c("ENSEMBL","SYMBOL"),c("ensembl","symbol"))
  m[,symbol:=toupper(symbol)]
  m <- m[,.SD[1L],by=ensembl]
  if (!nrow(m)) stop("No Ensembl ID mapped to a symbol")
  fwrite(m,fp,sep="\t")
  m
}

# Duplicate symbols are collapsed to the highest-variance row, so the retained
# probe is the informative one rather than whichever came first.
collapse_by_variance <- function(mat) {
  if (!anyDuplicated(rownames(mat))) return(mat)
  rv <- rowVars(mat,na.rm=TRUE)
  keep <- tapply(seq_along(rv),rownames(mat),function(i) i[which.max(rv[i])])
  mat[unlist(keep),,drop=FALSE]
}

amp_gistic <- 2     # GISTIC2 call counted as high-level amplification
min_cancers <- 5    # cancer types a gene needs before it is ranked

# Epigenetic regulators curated from EpiFactors, plus the H2A.Z chaperone and
# SRCAP/TIP60/INO80 complex members this study is concerned with.
epi_genes <- toupper(unique(c(
  "VPS72","SRCAP","EP400","RUVBL1","RUVBL2","ACTR5","ACTR6","ACTR8",
  "INO80","INO80B","INO80C","INO80D","INO80E","ZNHIT1","ANP32E",
  "H2AFZ","H2AFV","DMAP1","MRGBP","TRRAP","MEAF6","YEATS4","KAT5",
  "EPC1","EPC2","ING3","ING4","ING5","MBTD1","BRD8","ACTL6A","ACTL6B",
  "EZH1","EZH2","NSD1","NSD2","NSD3","SETD1A","SETD1B","SETD2","SETD3",
  "SETD4","SETD5","SETD6","SETD7","SETD8","SETDB1","SETDB2","SMYD1",
  "SMYD2","SMYD3","SMYD4","SMYD5","SUV39H1","SUV39H2","EHMT1","EHMT2",
  "DOT1L","KMT2A","KMT2B","KMT2C","KMT2D","KMT2E","ASH1L",
  "PRDM1","PRDM2","PRDM6","PRDM9","PRDM14","PRDM16",
  "KDM1A","KDM1B","KDM2A","KDM2B","KDM3A","KDM3B","KDM4A","KDM4B",
  "KDM4C","KDM4D","KDM5A","KDM5B","KDM5C","KDM5D","KDM6A","KDM6B",
  "KDM7A","KDM8","PHF2","PHF8","JHDM1D","JMJD3","RIOX1","RIOX2",
  "KAT2A","KAT2B","KAT6A","KAT6B","KAT7","KAT8","EP300",
  "CREBBP","KAT14","ELP3","HAT1",
  "HDAC1","HDAC2","HDAC3","HDAC4","HDAC5","HDAC6","HDAC7","HDAC8",
  "HDAC9","HDAC10","HDAC11","SIRT1","SIRT2","SIRT3","SIRT4","SIRT5",
  "SIRT6","SIRT7",
  "DNMT1","DNMT3A","DNMT3B","DNMT3L","TET1","TET2","TET3",
  "SMARCA2","SMARCA4","SMARCB1","SMARCC1","SMARCC2","SMARCD1","SMARCD2",
  "SMARCD3","SMARCE1","ARID1A","ARID1B","ARID2","PBRM1","BRD7","BRD9",
  "PHF10","DPF1","DPF2","DPF3","BCL7A","BCL7B","BCL7C",
  "MBD2","MBD3","CHD3","CHD4","CHD5","MTA1","MTA2","MTA3",
  "GATAD2A","GATAD2B","RBBP4","RBBP7",
  "CHD1","CHD2","CHD6","CHD7","CHD8","CHD9",
  "SMARCA1","SMARCA5","BPTF","BAZ1A","BAZ1B","BAZ2A","BAZ2B","CECR2","RSF1",
  "BMI1","RING1","RNF2","PCGF1","PCGF2","PCGF4","PCGF6",
  "CBX2","CBX4","CBX6","CBX7","CBX8","PHC1","PHC2","PHC3",
  "EED","SUZ12","JARID2","AEBP2","PHF1","PHF19","MTF2",
  "BRD2","BRD3","BRD4","BRDT","BRD1","BRD6","ATAD2","ATAD2B",
  "BRPF1","BRPF3","BRWD1","BRWD3","TAF1","ZMYND8","ZMYND11",
  "TRIM24","TRIM28","TRIM33",
  "NAP1L1","NAP1L2","NAP1L4","CHAF1A","CHAF1B","ASF1A","ASF1B",
  "HIRA","CABIN1","UBN1","DAXX","ATRX","DEK","SSRP1","SUPT16H",
  "RNF8","RNF168","USP3","USP16","USP22","BAP1","ASXL1","ASXL2","ASXL3",
  "MBD1","MBD4","MECP2","UHRF1","UHRF2")))

# Members of the complexes VPS72 sits in, flagged in the output so the ranking
# can be read against them.
complex_members <- c("VPS72","SRCAP","EP400","RUVBL1","RUVBL2","ACTR5","ACTR6",
                     "ACTR8","KAT5","DMAP1","MRGBP","TRRAP","YEATS4","ACTL6A",
                     "ACTL6B","INO80","INO80B","INO80C","INO80D","INO80E",
                     "H2AFZ","H2AFV","EPC1","EPC2","MBTD1","ZNHIT1","ANP32E","MEAF6")
