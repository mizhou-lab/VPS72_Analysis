# Collapse technical replicates to one column per patient, attach peak
# coordinates, and assign VPS72-high/low groups from matched RNA-seq.
source("00_config.R")

# The published matrix carries one column per sequencing run ("stan8" IDs, e.g.
# BRCA_A1AS); the identifier mapping links each run to a patient barcode.
atac <- fread(require_input(inputs$atac))
setnames(atac,1,"peak_id")

idmap <- fread(require_input(inputs$idmap))
setnames(idmap,c(find_col(idmap,c("Case_ID","case_id","patient")),
                 find_col(idmap,c("bam_prefix","stan8","run"))),c("case_id","run_id"))
idmap[,patient:=tcga12(case_id)]

run_cols <- setdiff(names(atac),"peak_id")
# Run IDs in the matrix header may carry suffixes; match on the leading token.
map <- idmap[match(sub("_.*$","",run_cols),sub("_.*$","",run_id))]
if (all(is.na(map$patient))) stop("No matrix column matched the identifier mapping; check the stan8 format")
cat("unmatched run columns:",sum(is.na(map$patient)),"of",length(run_cols),"\n")

mat <- as.matrix(atac[,..run_cols])
rownames(mat) <- atac$peak_id
mode(mat) <- "numeric"

# Replicate runs are averaged on the log2 scale, not re-summed from counts.
patients <- unique(na.omit(map$patient))
collapsed <- vapply(patients,function(p) rowMeans(mat[,which(map$patient==p),drop=FALSE],na.rm=TRUE),
                    numeric(nrow(mat)))
colnames(collapsed) <- patients
cat("runs:",length(run_cols),"-> patients:",ncol(collapsed),"\n")

# Coordinates travel with the peaks from here on, so later scripts never
# re-join the probeMap.
pm <- fread(require_input(inputs$probemap))
setnames(pm,1:4,c("peak_id","chrom","start","end"))
pm <- unique(pm[,.(peak_id,chrom,start=as.integer(start),end=as.integer(end))])

out <- merge(pm,data.table(peak_id=rownames(collapsed),as.data.table(collapsed)),by="peak_id")
if (!nrow(out)) stop("No peak ID matched between the matrix and the probeMap")
setcolorder(out,c("peak_id","chrom","start","end"))
write_tsv(out,"atac_matrix_with_coords.tsv")

rna <- fread(require_input(inputs$rna))
setnames(rna,1,"gene")
vps <- rna[grepl(paste0("^",vps72_ensg),gene)]
if (nrow(vps)!=1L) stop("Expected exactly one VPS72 row in the expression matrix, found ",nrow(vps))
rna_cols <- grep("^TCGA-",names(vps),value=TRUE)
if (!length(rna_cols)) stop("No TCGA sample columns in the expression matrix")
expr <- data.table(patient=tcga12(rna_cols),VPS72=as.numeric(vps[1L,..rna_cols]))
expr <- expr[is.finite(VPS72),.(VPS72=mean(VPS72)),by=patient]

samples <- merge(data.table(patient=patients),expr,by="patient",all.x=TRUE)
# Lineage comes from the stan8 prefix, e.g. "BRCA_A1AS" -> "BRCA".
samples <- merge(samples,unique(map[!is.na(patient),.(patient,cancer=sub("_.*$","",run_id))]),
                 by="patient",all.x=TRUE)

# Median split across the whole ATAC cohort: lineage is handled in the model
# (00_config.R: adjust_for_cancer_type) rather than in the grouping.
samples[,group:=fifelse(VPS72>=median(VPS72,na.rm=TRUE),"High","Low")]
samples[!is.finite(VPS72),group:=NA_character_]
if (!any(!is.na(samples$group))) stop("No patient has a VPS72 value; check barcode overlap")

write_tsv(samples,"sample_groups.tsv")
print(table(samples$group,useNA="ifany"))
