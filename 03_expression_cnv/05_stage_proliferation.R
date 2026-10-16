# Independent TCGA/Xena VPS72 stage and proliferation correlations.
# Uses Module11_Prolif_score explicitly; never substitutes another signature.
source("00_config.R")

xena <- fread(require_input(inputs$xena_merged))
score <- fread(require_input(inputs$xena_signatures))
find_col <- function(d,options) {
  z <- options[options %in% names(d)]
  if (!length(z)) stop("Required field missing: ",paste(options,collapse=", "))
  z[1L]
}
# Original Xena merged data may have underscore or mixed-case headers.
setnames(xena,names(xena),tolower(gsub("[^[:alnum:]]+","_",names(xena))))
samp_col <- find_col(xena,c("sample","samples"))
vps_col <- find_col(xena,c("vps72","vps72_average"))
ca_col <- find_col(xena,c("cancer_type_abbreviation","cancer"))
type_col <- find_col(xena,c("sample_type"))
stage_col <- find_col(xena,c("clinical_stage","ajcc_pathologic_tumor_stage"))
if (!"SetName" %in% names(score)) stop("Signature table must contain SetName")
if (!"Module11_Prolif_score" %in% score$SetName) stop("Module11_Prolif_score not found")
sig <- score[SetName=="Module11_Prolif_score"]
if (nrow(sig)!=1L) stop("Expected exactly one Module11_Prolif_score row")
sig_cols <- setdiff(names(sig),c("SetName","Source"))
if (!length(sig_cols)) stop("Proliferation signature table has no sample columns")
prolif <- data.table(sample=tcga15(sig_cols),prolif=as.numeric(unlist(sig[,..sig_cols],use.names=FALSE)))
prolif <- prolif[!is.na(sample)]
if (anyDuplicated(prolif$sample)) stop("Duplicate TCGA sample barcodes in proliferation scores")

parse_stage <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("^stage[[:space:]]*","",x)
  x <- gsub("[[:space:]]+","",x)
  x <- sub("[abc]$","",x)
  roman <- c(i=1L,ii=2L,iii=3L,iv=4L)
  as.integer(unname(roman[x]))
}
dat <- data.table(sample=tcga15(xena[[samp_col]]),
                  cancer=as.character(xena[[ca_col]]),
                  type=as.character(xena[[type_col]]),
                  stage=parse_stage(xena[[stage_col]]),
                  VPS72=suppressWarnings(as.numeric(xena[[vps_col]])))
dat <- dat[!is.na(sample) & !is.na(cancer) & nzchar(cancer)]
if (anyDuplicated(dat$sample)) stop("Duplicate TCGA samples in merged Xena data")
dat[,state:=fcase(type=="Primary Tumor","Primary",
                  type=="Solid Tissue Normal","Normal",
                  type=="Metastatic","Metastatic",
                  type=="Recurrent Tumor","Recurrent",
                  default="Other")]
dat <- merge(dat,prolif,by="sample",all.x=TRUE,sort=FALSE)

class_counts <- dat[is.finite(VPS72) & state!="Other",
                    .(n=.N,median_VPS72=median(VPS72),mean_VPS72=mean(VPS72)),by=state]
write_tsv(class_counts,"sample_class_vps72_summary.tsv")

# Per-cancer Spearman tests, with BH correction across cancer types separately
# for stage and proliferation endpoints.
primary <- dat[state=="Primary" & is.finite(VPS72)]
if (!nrow(primary)) stop("No primary tumor samples with VPS72")
stage <- if (any(!is.na(primary$stage))) primary[!is.na(stage),{
  as.list(safe_spearman(VPS72,as.numeric(stage),min_n=min_stage_n))
},by=cancer] else data.table(cancer=character(),n=integer(),rho=numeric(),p=numeric())
stage[,FDR:=p.adjust(p,method="BH")]
write_tsv(stage,"per_cancer_vps72_stage_spearman.tsv")

prolif_tests <- primary[is.finite(prolif),{
  as.list(safe_spearman(VPS72,prolif,min_n=min_stage_n))
},by=cancer]
prolif_tests[,FDR:=p.adjust(p,method="BH")]
write_tsv(prolif_tests,"per_cancer_vps72_proliferation_spearman.tsv")

pooled <- rbindlist(list(
  as.list(safe_spearman(primary$VPS72,as.numeric(primary$stage),min_stage_n)),
  as.list(safe_spearman(primary$VPS72,primary$prolif,min_stage_n))
))
pooled[,endpoint:=c("stage","Module11_Prolif_score")]
write_tsv(pooled,"pooled_vps72_associations.tsv")
cat("Xena primary tumors:",nrow(primary),"; cancer groups:",uniqueN(primary$cancer),"\n")
