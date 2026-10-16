# Sample-level VPS72 expression, clinical covariates and non-silent mutation calls.
# Written once here so every model in 05-07 is fitted on the same sample set.
source("00_config.R")

expr <- fread(require_input(inputs$expr))
setnames(expr,1,"gene_id")
vps <- expr[grepl("^VPS72\\|",gene_id,ignore.case=TRUE)]
if (nrow(vps)!=1L) stop("Expected exactly one VPS72 row in the expression matrix, found ",nrow(vps))
expr_cols <- setdiff(names(vps),"gene_id")
# Several aliquots can map to one patient; averaged.
vps_expr <- data.table(patient=tcga12(expr_cols),VPS72=as.numeric(vps[1L,..expr_cols]))
vps_expr <- vps_expr[is.finite(VPS72),.(VPS72=mean(VPS72)),by=patient]
rm(expr); gc()

clin <- fread(require_input(inputs$clin))
clin <- unique(data.table(
  patient = tcga12(clin[[find_col(clin,c("Patient ID","patient_id","bcr_patient_barcode"))]]),
  age     = suppressWarnings(as.numeric(clin[[find_col(clin,c("Diagnosis Age","age_at_initial_pathologic_diagnosis","age"))]])),
  sex     = as.character(clin[[find_col(clin,c("Sex","gender"))]]),
  cancer  = as.character(clin[[find_col(clin,c("Cancer Type","cancer_type","TCGA PanCanAtlas Cancer Type Acronym"))]])))
clin <- clin[!is.na(patient) & nzchar(patient)]
clin <- clin[,.SD[1L],by=patient]

mc3 <- read_mc3()
muts <- unique(mc3[,.(patient,gene)])
muts[,mutated:=1L]
mut_mat <- dcast(muts,patient~gene,value.var="mutated",fill=0L)
rm(mc3,muts); gc()

master <- merge(vps_expr,clin,by="patient")
master <- master[!is.na(age) & !is.na(sex) & nzchar(sex) & !is.na(cancer) & nzchar(cancer)]
if (!nrow(master)) stop("No patients survived the covariate join; check the clinical export columns")
master[,VPS72_z:=as.numeric(scale(VPS72))]

# Both tables are restricted to the same patients and ordered identically, so
# downstream scripts can bind columns without re-matching.
keep <- intersect(master$patient,mut_mat$patient)
master <- master[patient %chin% keep]; mut_mat <- mut_mat[patient %chin% keep]
setorder(master,patient); setorder(mut_mat,patient)
if (!identical(master$patient,mut_mat$patient)) stop("Master and mutation tables are misaligned")

write_tsv(master,"master_sample_table.tsv")
write_tsv(mut_mat,"mutation_matrix.tsv")
cat("Patients:",nrow(master),"; cancer types:",uniqueN(master$cancer),
    "; genes mutated in at least one patient:",ncol(mut_mat)-1L,"\n")
