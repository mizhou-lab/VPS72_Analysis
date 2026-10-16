# Pan-cancer disruption ranking of the epigenetic regulator panel.
# Composite score = z(mean expression) + z(mean copy number) + z(amplification
# frequency), each z-scored within the panel, every cancer type weighted equally.
source("00_config.R")

annot <- read_sample_annotation()
cat("annotated patients:",nrow(annot),"across",uniqueN(annot$cancer),"cancer types\n")

## ---- Copy number ------------------------------------------------------------
gistic <- fread(require_input(inputs$gistic),header=TRUE)
setnames(gistic,1,"gene_id")
ids <- gsub("\\.[0-9]+$","",trimws(as.character(gistic$gene_id)))
map <- ensembl_to_symbol(ids)
syms <- map$symbol[match(ids,map$ensembl)]
idx <- which(!is.na(syms) & syms %chin% epi_genes)
if (!length(idx)) stop("No panel gene mapped into the GISTIC table")
cat("panel genes in GISTIC:",uniqueN(syms[idx]),"of",length(epi_genes),"\n")

g_mat <- as.matrix(gistic[idx,-1L,with=FALSE])
mode(g_mat) <- "numeric"
rownames(g_mat) <- syms[idx]
colnames(g_mat) <- tcga12(colnames(g_mat))
g_mat <- collapse_by_variance(g_mat)
rm(gistic); gc()

keep <- intersect(colnames(g_mat),annot$patient)
if (!length(keep)) stop("No GISTIC barcode matched the sample annotation")
g_mat <- g_mat[,keep,drop=FALSE]
cat("GISTIC samples retained:",length(keep),"\n")

g_long <- data.table(gene=rep(rownames(g_mat),length(keep)),
                     patient=rep(colnames(g_mat),each=nrow(g_mat)),
                     gistic=as.vector(g_mat))
g_long <- merge(g_long,annot,by="patient")[is.finite(gistic)]
rm(g_mat); gc()

# Per gene and cancer type first, then averaged over cancer types, so a large
# cohort does not dominate the pan-cancer figure.
cnv_ct <- g_long[,.(cnv_mean_ct=mean(gistic),
                    amp_freq_ct=mean(gistic>=amp_gistic)),by=.(gene,cancer)]
cnv <- cnv_ct[,.(cnv_mean=mean(cnv_mean_ct),amp_mean=mean(amp_freq_ct),
                 n_cancers_cnv=.N),by=gene]
write_tsv(cnv_ct,"cnv_by_cancer_type.tsv")
rm(g_long); gc()

## ---- Expression -------------------------------------------------------------
expr <- fread(require_input(inputs$expr),header=TRUE)
setnames(expr,1,"gene_id")
e_syms <- toupper(gsub("\\|.*$","",trimws(as.character(expr$gene_id))))
eidx <- which(e_syms %chin% epi_genes)
if (!length(eidx)) stop("No panel gene found in the expression matrix")
cat("panel genes in expression:",uniqueN(e_syms[eidx]),"of",length(epi_genes),"\n")

e_mat <- as.matrix(expr[eidx,-1L,with=FALSE])
mode(e_mat) <- "numeric"
rownames(e_mat) <- e_syms[eidx]
colnames(e_mat) <- tcga12(colnames(e_mat))
e_mat <- collapse_by_variance(e_mat)
rm(expr); gc()

# z-scored per gene across all samples, so the matrix is comparable between
# genes before any per-cancer summary is taken.
e_z <- (e_mat-rowMeans(e_mat,na.rm=TRUE))/rowSds(e_mat,na.rm=TRUE)
keep_e <- intersect(colnames(e_z),annot$patient)
if (!length(keep_e)) stop("No expression barcode matched the sample annotation")
e_z <- e_z[,keep_e,drop=FALSE]
cat("expression samples retained:",length(keep_e),"\n")

e_long <- data.table(gene=rep(rownames(e_z),length(keep_e)),
                     patient=rep(colnames(e_z),each=nrow(e_z)),
                     expr_z=as.vector(e_z))
e_long <- merge(e_long,annot,by="patient")[is.finite(expr_z)]
rm(e_mat,e_z); gc()

expr_ct <- e_long[,.(expr_median_z=median(expr_z)),by=.(gene,cancer)]
expr <- expr_ct[,.(expr_mean=mean(expr_median_z),n_cancers_expr=.N),by=gene]
write_tsv(expr_ct,"expression_by_cancer_type.tsv")
rm(e_long); gc()

## ---- Composite --------------------------------------------------------------
dt <- merge(expr,cnv,by="gene",all=TRUE)

# A gene measured in too few cancer types has an unstable mean, so it is held
# out of the ranking rather than scored on thin data.
dt[,n_cancers_min:=pmin(n_cancers_expr,n_cancers_cnv,na.rm=TRUE)]
dropped <- dt[is.na(n_cancers_min) | n_cancers_min<min_cancers]
if (nrow(dropped)) {
  write_tsv(dropped[,.(gene,n_cancers_expr,n_cancers_cnv)],"genes_excluded.tsv")
  cat("excluded for fewer than",min_cancers,"cancer types:",nrow(dropped),"genes\n")
}
dt <- dt[!is.na(n_cancers_min) & n_cancers_min>=min_cancers]
if (nrow(dt)<20L) stop("Only ",nrow(dt)," genes survived; check the barcode overlap")

dt[,`:=`(expr_s=as.numeric(scale(expr_mean)),
         cnv_s=as.numeric(scale(cnv_mean)),
         amp_s=as.numeric(scale(amp_mean)))]
dt[,n_components:=(!is.na(expr_s))+(!is.na(cnv_s))+(!is.na(amp_s))]

# A missing component is substituted with zero, which in z-space is the panel
# mean. That imputes "average" rather than excluding the gene, so a gene scored
# on fewer components is not directly comparable; n_components records which.
dt[,disruption_score:=rowSums(cbind(fifelse(is.na(expr_s),0,expr_s),
                                    fifelse(is.na(cnv_s),0,cnv_s),
                                    fifelse(is.na(amp_s),0,amp_s)))]
dt[,rank:=rank(-disruption_score,ties.method="min")]
dt[,percentile:=1-(rank-1)/.N]
# Fraction of panel genes scoring at or above each gene: a position within this
# panel, not a test against a null model.
dt[,empirical_p:=vapply(disruption_score,
                        function(x) mean(disruption_score>=x,na.rm=TRUE),numeric(1))]
dt[,complex_member:=gene %chin% complex_members]
setorder(dt,rank)

write_tsv(dt,"epi_regulator_disruption_ranking.tsv")

## ---- Report -----------------------------------------------------------------
cat("\nGenes ranked:",nrow(dt),"| complex members among them:",sum(dt$complex_member),
    "| scored on all three components:",sum(dt$n_components==3L),"\n\n")

v <- dt[gene=="VPS72"]
if (nrow(v)==1L) {
  cat(sprintf("VPS72: rank %d of %d (%.0fth percentile), score %.3f [expr %.3f, cnv %.3f, amp %.3f]\n\n",
              v$rank,nrow(dt),v$percentile*100,v$disruption_score,v$expr_s,v$cnv_s,v$amp_s))
} else cat("VPS72 is not in the ranking.\n\n")

print(dt[1:20,.(rank,gene,disruption_score=round(disruption_score,3),
                expr_s=round(expr_s,3),cnv_s=round(cnv_s,3),amp_s=round(amp_s,3),
                empirical_p=round(empirical_p,4),complex_member)])

cat("\nComplex members, ranked:\n")
print(dt[complex_member==TRUE,.(rank,gene,disruption_score=round(disruption_score,3),
                                percentile=round(percentile*100))])
