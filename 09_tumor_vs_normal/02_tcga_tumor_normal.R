# TCGA primary tumour and adjacent solid-tissue normal samples, with VPS72
# expression and the patient pairing used by the sensitivity analysis.
source("00_config.R")

pheno <- fread(require_input(inputs$xena_pheno))
setnames(pheno,c(find_col(pheno,c("sample","sampleID")),
                 find_col(pheno,c("_sample_type","sample_type")),
                 find_col(pheno,c("_primary_site","primary_site","detailed_category"))),
         c("sample","sample_type","site"))
proj_col <- find_col(pheno,c("_study","project","study"))

# Xena carries TCGA, TARGET and GTEx in one phenotype table.
pheno <- pheno[get(proj_col)=="TCGA" | grepl("^TCGA",sample)]
pheno[,state:=fcase(sample_type=="Primary Tumor","Tumor",
                    sample_type=="Solid Tissue Normal","Normal",
                    default=NA_character_)]
pheno <- pheno[!is.na(state)]
if (!nrow(pheno)) stop("No primary tumour or adjacent normal samples in the phenotype table")

# Cancer code from the detailed category when present, otherwise from the
# cohort column; samples that cannot be mapped to an eligible project are cut.
code_col <- find_col(pheno,c("cancer type abbreviation","_cohort","cohort","disease_code"))
pheno[,cancer:=toupper(gsub(".*\\(|\\).*","",as.character(get(code_col))))]
pheno[!cancer %chin% eligible,
      cancer:=toupper(sub("^TCGA[ _-]*","",as.character(get(code_col))))]
pheno <- pheno[cancer %chin% eligible]
if (!nrow(pheno)) stop("No sample mapped to an eligible TCGA project; check the cohort column spelling")

pheno[,patient:=tcga12(sample)]
pheno <- unique(pheno[,.(sample=tcga15(sample),patient,cancer,state,site)])
pheno <- pheno[,.SD[1L],by=sample]

# Only the VPS72 row is read out of the expression matrix: the full matrix is
# large and nothing else is needed here.
hdr <- names(fread(require_input(inputs$xena_expr),nrows=0L))
rows <- fread(inputs$xena_expr,select=1L)
setnames(rows,1,"gene")
idx <- grep(paste0("^",vps72_ensg),rows$gene)
if (length(idx)!=1L) stop("Expected exactly one VPS72 row in the expression matrix, found ",length(idx))
vps <- fread(inputs$xena_expr,skip=idx,nrows=1L,header=FALSE)
setnames(vps,hdr)
setnames(vps,1,"gene")

cols <- intersect(setdiff(hdr,hdr[1L]),names(vps))
expr <- data.table(sample=tcga15(cols),VPS72=as.numeric(vps[1L,..cols]))
expr <- expr[is.finite(VPS72),.(VPS72=mean(VPS72)),by=sample]

tcga <- merge(pheno,expr,by="sample")
if (!nrow(tcga)) stop("No sample overlap between the phenotype and expression tables")

# Patients contributing both a tumour and an adjacent normal.
paired <- tcga[,.(states=uniqueN(state)),by=patient][states==2L,patient]
tcga[,paired:=patient %chin% paired]

# One sample per patient and state, so patient-level summaries are not weighted
# by patients who happen to have several aliquots.
setorder(tcga,patient,state,sample)
tcga[,primary_record:=seq_len(.N)==1L,by=.(patient,state)]

write_tsv(tcga,"tcga_tumor_normal.tsv")
cat("Samples:",nrow(tcga),"; patients:",uniqueN(tcga$patient),
    "; paired patients:",length(paired),"\n")
print(dcast(tcga[,.N,by=.(cancer,state)],cancer~state,value.var="N",fill=0L))
