# TCGA PanCancer Atlas cBioPortal clinical + TCGA/GTEx TPM expression signatures.
# Computes per-cancer continuous and median-split Cox estimates, without plots/Excel.
# Preserves gene sets and arithmetic-mean TPM scores from the uploaded source.
source("00_config.R"); need(c("data.table","survival"))
f <- required(TCGA_TPM_MATRIX_FILE); cl <- data.table::fread(required(CBIO_CLIN_FILE))
sig <- list(
 SRCAP=c("SRCAP","VPS72","RUVBL1","RUVBL2","DMAP1","ACTR5","ACTR6","ACTL6A","YEATS4","ZNHIT1","ARID1A"),
 TIP60=c("KAT5","EPC1","EPC2","VPS72","DMAP1","RUVBL1","RUVBL2","ING3","TRRAP"),
 INO80=c("INO80","RUVBL1","RUVBL2","ACTR5","ACTR6","MCRS1","NFRKB","UCHL5","YY1"))
all_genes <- unique(c(unlist(sig),"VPS72"))
expr <- data.table::fread(f)
if(ncol(expr)<3) stop("TPM matrix must be genes x samples")
gc <- names(expr)[1]; expr[,GENE_SYMBOL:=toupper(trimws(as.character(get(gc))))]
col_samples <- setdiff(names(expr),c(gc,"GENE_SYMBOL"))
col_samples <- col_samples[!is.na(tcga15(col_samples)) & substr(tcga15(col_samples),14,15)=="01"]
expr <- expr[GENE_SYMBOL %in% all_genes]
if(!"VPS72" %in% expr$GENE_SYMBOL) stop("No VPS72 gene row in TPM matrix")
if(anyDuplicated(expr$GENE_SYMBOL)) stop("Duplicated TPM gene symbols; resolve mapping before scoring")
if(!length(col_samples)) stop("No TCGA primary-tumor sample columns in TPM matrix")
mat <- as.matrix(expr[,..col_samples]); storage.mode(mat) <- "double"
d <- data.table::as.data.table(t(mat)); data.table::setnames(d,expr$GENE_SYMBOL)
d[,sample:=col_samples]; d[,`:=`(patient=tcga12(sample),sample15=tcga15(sample))]
d <- d[!is.na(patient)][order(patient,sample)][!duplicated(patient)]
# cBioPortal merged clinical tables can contain multiple studies per patient;
# use PanCancer Atlas records when available instead of joining all exports blindly.
id <- pick_col(names(cl),c("Sample ID","Patient ID")); study <- pick_col(names(cl),c("Study ID"))
time <- pick_col(names(cl),c("Overall Survival (Months)","OS Months","OS_MONTHS"))
event <- pick_col(names(cl),c("Overall Survival Status","OS Status","OS_STATUS"))
cl <- cl[grepl("_tcga_pan_can_atlas_2018$",as.character(cl[[study]])),]
if(!nrow(cl)) stop("No TCGA PanCancer Atlas 2018 clinical records; check clinical export")
c <- data.table::data.table(patient=tcga12(cl[[id]]),cancer=toupper(sub("_tcga.*$","",as.character(cl[[study]]))),
                            time_months=as.numeric(cl[[time]]),event=parse_event(cl[[event]]))
c <- c[!is.na(patient)&is.finite(time_months)&time_months>0&!is.na(event)]
if(c[,.(n_cancer=uniqueN(cancer),n_time=uniqueN(round(time_months,5)),n_event=uniqueN(event)),by=patient][n_cancer>1|n_time>1|n_event>1,.N]) stop("Conflicting clinical study/time/event records within patients")
c <- unique(c,by="patient"); d <- merge(d,c,by="patient",all=FALSE)
if(!nrow(d)) stop("No matched samples in TPM and cBioPortal clinical files")
for(s in names(sig)) {
  avail <- intersect(sig[[s]],names(d)); if(length(avail)<3) stop(s," signature has fewer than 3 mapped genes")
  d[[paste0(s,"_score")]] <- rowMeans(as.matrix(d[,..avail]),na.rm=TRUE)
  message(s,": ",length(avail),"/",length(sig[[s]])," genes mapped")
}
if("VPS72" %in% names(d)) {
  avail <- setdiff(intersect(sig$SRCAP,names(d)),"VPS72")
  d[,SRCAP_noVPS72:=rowMeans(as.matrix(d[,..avail]),na.rm=TRUE)]
  d[,SRCAP_minus_VPS72:=SRCAP_noVPS72 - VPS72]
}
write_table(d,"04_TCGA_signature_survival_samples.tsv")
out <- list(); markers <- c("VPS72","SRCAP_score","TIP60_score","INO80_score","SRCAP_minus_VPS72")
for(ca in sort(unique(d$cancer))) for(v in markers) {
  z <- d[cancer==ca & is.finite(get(v))]; if(!nrow(z)) next
  a <- run_cox(z,v,label="continuous_original_TPM_units")
  if(!is.null(a)) out[[length(out)+1]] <- cbind(cancer=ca,signature=v,a)
  z[,median_group:=as.integer(get(v)>=median(get(v)))]; a <- run_cox(z,"median_group",label="within_cancer_median_split")
  if(!is.null(a)) out[[length(out)+1]] <- cbind(cancer=ca,signature=v,a)
}
if(!length(out)) stop("No evaluable per-cancer signature-survival models")
r <- data.table::rbindlist(out,fill=TRUE); r[,FDR:=p.adjust(p,"BH"),by=.(signature,model)]
write_table(r,"04_TCGA_signature_survival_Cox_models.tsv")
