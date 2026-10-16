# VPS72 survival analysis: common configuration and helpers (R >= 4.2).
# Set INPUT_DIR to a directory containing the public-source files. No personal paths required.
INPUT_DIR <- Sys.getenv("VPS72_SURVIVAL_INPUT", unset="data")
OUTPUT_DIR <- Sys.getenv("VPS72_SURVIVAL_OUTPUT", unset="results")
dir.create(OUTPUT_DIR, recursive=TRUE, showWarnings=FALSE)
input <- function(name) file.path(INPUT_DIR, name)
output <- function(name) file.path(OUTPUT_DIR, name)
required <- function(name) { p <- input(name); if (!file.exists(p)) stop("Missing input: ",p,call.=FALSE); p }
need <- function(pkgs) { absent <- pkgs[!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)]; if(length(absent)) stop("Install required R packages: ",paste(absent,collapse=", "),call.=FALSE) }
write_table <- function(x,name) { data.table::fwrite(x, output(name),sep="\t",na="NA"); message("Wrote ",output(name)) }
clean_id <- function(x) toupper(gsub("\\.","-",trimws(as.character(x))))
tcga12 <- function(x) { z <- clean_id(x); y <- rep(NA_character_,length(z)); ok <- grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}",z); y[ok] <- substr(z[ok],1,12); y }
tcga15 <- function(x) { z <- clean_id(x); y <- rep(NA_character_,length(z)); ok <- grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}",z); y[ok] <- substr(z[ok],1,15); y }
pick_col <- function(nms,variants,required=TRUE) {
  normalize <- function(z) tolower(gsub("[^a-z0-9]","",z))
  ix <- match(normalize(variants),normalize(nms)); ix <- ix[!is.na(ix)]
  if (!length(ix)) { if(required) stop("Missing column (expected one of): ",paste(variants,collapse=", ")); return(NA_character_) }
  nms[ix[1]]
}
parse_event <- function(x) {
  y <- toupper(trimws(as.character(x)))
  out <- rep(NA_integer_,length(y))
  out[grepl("^(1($|:)|DECEASED$|DEAD$|EVENT$)",y)] <- 1L
  out[grepl("^(0($|:)|LIVING$|ALIVE$|CENSORED$)",y)] <- 0L
  out
}
first_nonmissing <- function(x) { x <- x[!is.na(x) & trimws(as.character(x))!=""]; if(length(x)) x[1] else NA }
run_cox <- function(d, predictor, covariates=character(), strata=NULL, label="") {
  vars <- unique(c("time_months","event",predictor,covariates,strata))
  d <- as.data.frame(d)
  d <- d[stats::complete.cases(d[,vars,drop=FALSE]),,drop=FALSE]
  d <- d[is.finite(d$time_months) & d$time_months>0 & d$event %in% 0:1,,drop=FALSE]
  if(nrow(d)<30 || sum(d$event)<10 || length(unique(d[[predictor]]))<2) return(NULL)
  if(length(covariates) && any(vapply(d[,covariates,drop=FALSE],function(x) length(unique(x))<2,logical(1)))) return(NULL)
  f <- stats::as.formula(paste("survival::Surv(time_months,event) ~",paste(c(predictor,covariates,if(!is.null(strata)) paste0("strata(",strata,")")),collapse=" + ")))
  fit <- tryCatch(survival::coxph(f,data=d,x=TRUE),error=function(e) NULL)
  if(is.null(fit)) return(NULL)
  z <- summary(fit)
  if(!predictor %in% rownames(z$coefficients)) return(NULL)
  q <- z$coefficients[predictor,,drop=FALSE]; ci <- z$conf.int[predictor,,drop=FALSE]
  data.frame(model=label,term=predictor,n=nrow(d),events=sum(d$event),HR=unname(ci[1,"exp(coef)"]),
     ci_low=unname(ci[1,"lower .95"]),ci_high=unname(ci[1,"upper .95"]),p=unname(q[1,"Pr(>|z|)"]),
     formula=paste(deparse(f),collapse=""),stringsAsFactors=FALSE)
}
safe_wilcox <- function(x,y) { x <- x[is.finite(x)]; y <- y[is.finite(y)]; if(length(x)<3||length(y)<3) return(NA_real_); suppressWarnings(stats::wilcox.test(x,y,exact=FALSE)$p.value) }
# Cohort lists come from the uploaded survival scripts; they are *fixed exploratory subsets*.
SIG_UP <- c("BRCA","LUAD","LUSC","LIHC","COAD","UCEC","HNSC","THCA","CHOL","READ","ESCA","PRAD","KIRC","GBM","BLCA","STAD","KIRP")
UP <- c(SIG_UP,"CESC","PCPG","PAAD","SARC")
LIQUID <- c("LAML","DLBC")
# Input file names are configurable here only; other scripts source this file.
XENA_EXPR_FILE <- "VPS72expression pan cancer.tsv"
XENA_CLIN_FILE <- "Clinical.tsv"
CBIO_CLIN_FILE <- "combined_study_clinical_data.tsv"
CBIO_EXPR_FILE <- "VPS72__mRNA_Expression,_RSEM_(Batch_normalized_from_Illumina_HiSeq_RNASeqV2) (1).txt"
CBIO_CNV_FILE <- "VPS72__Log2_copy-number_values.txt"
TCGA_TPM_MATRIX_FILE <- "TcgaTargetGtex_rsem_gene_tpm_annotated_bySymbol.csv.gz"
PCAWG_CLIN_FILE <- "pancan_pcawg_2020_clinical_data.tsv"
PCAWG_EXPR_FILE <- "VPS72__mRNA_expression_(FPKM_UQ).txt"
GLIOMA_CLIN_FILE <- "combined_study_clinical_data (2).tsv"
GLIOMA_MUT_FILES <- c("plot (9).txt","plot (10).txt")
GLIOMA_PROLIF_FILE <- "Scores_160_Signatures.tsv"
