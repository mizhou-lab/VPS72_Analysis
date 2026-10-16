# Paths, thresholds and shared helpers for the DepMap 25Q2 VPS72 analyses.
# Set PROJ_DIR; every path below is built from it.
suppressPackageStartupMessages({library(data.table)})

PROJ_DIR <- "."
DATA_DIR <- file.path(PROJ_DIR,"data","DepMap")
OUT_DIR  <- file.path(PROJ_DIR,"output")
dir.create(OUT_DIR,recursive=TRUE,showWarnings=FALSE)

# DepMap 25Q2. Portal gene-level exports and the full omics matrices are both
# used; file names are the portal's own and are not distributed here.
inputs <- list(
  model        = file.path(DATA_DIR,"Model.csv"),
  profiles     = file.path(DATA_DIR,"OmicsProfiles.csv"),
  global_sig   = file.path(DATA_DIR,"OmicsGlobalSignaturesProfile.csv"),
  expr_gene    = file.path(DATA_DIR,"VPS72_Expression_Public_25Q2.csv"),
  cnv_gene     = file.path(DATA_DIR,"VPS72_CopyNumber_Public_25Q2_Log2.csv"),
  crispr       = file.path(DATA_DIR,"VPS72_CRISPR_Chronos_25Q2.csv"),
  rnai         = file.path(DATA_DIR,"VPS72_RNAi_DEMETER2.csv"),
  crispr_codep = file.path(DATA_DIR,"VPS72_Top100_Codependencies_CRISPR.csv"),
  rnai_codep   = file.path(DATA_DIR,"VPS72_Top100_Codependencies_RNAi.csv"),
  expr_matrix  = file.path(DATA_DIR,"OmicsExpressionProteinCodingGenesTPMLogp1.csv"),
  cnv_matrix   = file.path(DATA_DIR,"OmicsCNGeneWGS.csv"),
  crispr_matrix= file.path(DATA_DIR,"CRISPRGeneEffect.csv"),
  dep          = file.path(OUT_DIR,"vps72_dependency.tsv"),
  scores       = file.path(OUT_DIR,"hallmark_program_scores.tsv")
)

require_input <- function(p) {
  if (!file.exists(p)) stop("Cannot see required input: ",p,
                            "\nDownload it from the DepMap portal into data/DepMap/ (see README), ",
                            "or run the script that writes it.")
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

clean_names <- function(d) { setnames(d,names(d),gsub("[^A-Za-z0-9_]","_",names(d))); d }

# Portal single-gene exports carry the model ID first and the value second,
# under a header that includes the release name.
read_gene_export <- function(p,value_name) {
  d <- clean_names(fread(require_input(p)))
  if (ncol(d)<2L) stop("Expected at least two columns in ",p)
  setnames(d,1:2,c("ModelID",value_name))
  d <- d[,.SD,.SDcols=c("ModelID",value_name)]
  d[[value_name]] <- suppressWarnings(as.numeric(d[[value_name]]))
  unique(d[!is.na(ModelID) & nzchar(ModelID)])
}

read_model <- function() {
  m <- clean_names(fread(require_input(inputs$model)))
  out <- data.table(
    ModelID = as.character(m[[find_col(m,c("ModelID","DepMap_ID","model_id"))]]),
    name    = as.character(m[[find_col(m,c("CellLineName","StrippedCellLineName"))]]),
    lineage = as.character(m[[find_col(m,c("OncotreeLineage","Lineage"))]]),
    disease = as.character(m[[find_col(m,c("OncotreePrimaryDisease","PrimaryDisease"))]]))
  out[,is_cancer:=fifelse(disease=="Non-Cancerous","Normal","Cancer")]
  out[nzchar(ModelID)][,.SD[1L],by=ModelID]
}

# Gene symbols in the omics matrices are "SYMBOL (ENTREZ)".
strip_entrez <- function(x) toupper(trimws(sub("\\s*\\(.*$","",x)))

# Dependency thresholds. DepMap's own convention calls a Chronos score of -1
# the median of known common-essential genes, so -1 marks "as essential as a
# core essential gene", while -0.5 is the looser cut often used for "depleted".
# Both are reported everywhere, because the two give very different pictures
# for VPS72 and the choice of line on a density plot is not cosmetic.
dep_strong <- -1.0
dep_loose  <- -0.5

codep_cut   <- 0.3   # |correlation| called a co-dependency in the published figure
fdr_cut     <- 0.05
min_lineage_n <- 10  # cell lines in a lineage before it is correlated
min_set_n   <- 5     # member genes present before a Hallmark set is scored

instability_metrics <- c("CIN","Aneuploidy","Ploidy","LoHFraction","MSIScore")

# Members of the complexes VPS72 sits in, flagged in the co-dependency output.
complex_members <- c("SRCAP","EP400","KAT5","TRRAP","DMAP1","RUVBL1","RUVBL2",
                     "ACTL6A","ACTL6B","YEATS4","ING3","MRGBP","ACTR5","ACTR6",
                     "INO80","H2AZ1","H2AZ2","H2AFZ","H2AFV")

spearman <- function(x,y,min_n=10L) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok)<min_n) return(NULL)
  ct <- suppressWarnings(cor.test(x[ok],y[ok],method="spearman",exact=FALSE))
  data.table(n=sum(ok),rho=as.numeric(ct$estimate),p=ct$p.value)
}
