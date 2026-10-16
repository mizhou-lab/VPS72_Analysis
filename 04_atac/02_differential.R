# Differential accessibility, VPS72-high vs VPS72-low, by limma.
# Zero-variance peaks are dropped; p-values adjusted by Benjamini-Hochberg.
source("00_config.R")

atac <- fread(require_input(inputs$coords))
coords <- atac[,.(peak_id,chrom,start,end)]

samp <- fread(require_input(inputs$groups))[!is.na(group)]
sample_cols <- intersect(samp$patient,names(atac))
if (!length(sample_cols)) stop("No sample column in the matrix matches the group table")
samp <- samp[match(sample_cols,patient)]
samp[,group:=factor(group,levels=c("Low","High"))]

mat <- as.matrix(atac[,..sample_cols])
rownames(mat) <- atac$peak_id
keep <- rowVars(mat,na.rm=TRUE)>0
mat <- mat[keep,]
if (!nrow(mat)) stop("Every peak has zero variance")

# Lineages represented by a single sample carry no within-group information and
# are dropped when the lineage term is in the model.
adjust <- adjust_for_cancer_type && "cancer" %in% names(samp) && uniqueN(na.omit(samp$cancer))>1L
if (adjust) {
  ok <- samp$cancer %chin% names(which(table(samp$cancer)>1L))
  samp <- samp[ok]; mat <- mat[,ok,drop=FALSE]
  design <- model.matrix(~group+factor(cancer),data=samp)
  cat("model: ~ group + cancer (",uniqueN(samp$cancer),"lineages, n =",ncol(mat),")\n")
} else {
  design <- model.matrix(~group,data=samp)
  cat("model: ~ group (unadjusted, n =",ncol(mat),")\n")
}
if (!"groupHigh" %in% colnames(design)) stop("Design matrix has no groupHigh column; check group levels")

fit <- eBayes(lmFit(mat,design))
res <- as.data.table(topTable(fit,coef="groupHigh",number=Inf,sort.by="none"),keep.rownames="peak_id")
setnames(res,c("logFC","P.Value","adj.P.Val"),c("log2FC","pvalue","padj"))
res <- merge(coords,res,by="peak_id")
res[,direction:=fifelse(log2FC>0,"up_in_high","down_in_high")]
setorder(res,pvalue)

write_tsv(res,"differential_accessibility.tsv")
cat("Peaks tested:",nrow(res),"; at FDR <",fdr_cut,":",sum(res$padj<fdr_cut,na.rm=TRUE),"\n")
