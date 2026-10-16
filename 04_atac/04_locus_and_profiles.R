# Three summary tables over the differential results: locus means by cancer
# type and VPS72 group, binned effect size around TSSs, and per-gene-set
# directionality with exact binomial tests.
source("00_config.R")
suppressPackageStartupMessages({library(TxDb.Hsapiens.UCSC.hg38.knownGene)})

atac <- fread(require_input(inputs$coords))
samp <- fread(require_input(inputs$groups))[!is.na(group)]
sample_cols <- intersect(samp$patient,names(atac))
if (!length(sample_cols)) stop("No sample column in the matrix matches the group table")
samp <- samp[match(sample_cols,patient)]

locus_means <- rbindlist(lapply(names(loci),function(nm) {
  L <- loci[[nm]]
  sub <- atac[chrom==L$chrom & start>=L$start & end<=L$end]
  if (!nrow(sub)) { cat("  no peaks in locus:",nm,"\n"); return(NULL) }
  per_sample <- colMeans(as.matrix(sub[,..sample_cols]),na.rm=TRUE)
  data.table(locus=nm,cancer=samp$cancer,group=samp$group,acc=per_sample)[
    ,.(mean_accessibility=mean(acc,na.rm=TRUE),sd_accessibility=sd(acc,na.rm=TRUE),n_samples=.N),
    by=.(locus,cancer,group)]
}))
if (!length(locus_means) || !nrow(locus_means)) stop("No locus had peaks in the matrix")
# Cells below min_locus_n give unstable means; dropped rather than reported.
locus_means <- locus_means[n_samples>=min_locus_n]
setorder(locus_means,locus,cancer,group)
write_tsv(locus_means,"locus_means.tsv")

da <- fread(require_input(inputs$da))
sig <- da[padj<fdr_cut & abs(log2FC)>lfc_annot]
if (!nrow(sig)) stop("No peak passes FDR < ",fdr_cut," and |log2FC| > ",lfc_annot)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
tss <- suppressWarnings(resize(transcripts(txdb),width=1,fix="start"))
mid <- GRanges(sig$chrom,IRanges(round((sig$start+sig$end)/2),width=1))

# Binning is the analysis here; the profile plot is presentation.
prof <- data.table(distance=mcols(distanceToNearest(mid,tss))$distance,log2FC=sig$log2FC)
prof <- prof[!is.na(distance) & distance<=tss_window]
if (!nrow(prof)) stop("No differential peak lies within ",tss_window," bp of a TSS")
prof[,bin:=floor(distance/tss_bin)*tss_bin]
tss_profile <- prof[,.(n_peaks=.N,mean_log2FC=mean(log2FC),sem_log2FC=sd(log2FC)/sqrt(.N)),
                    by=.(distance_to_tss_bp=bin)][order(distance_to_tss_bp)]
write_tsv(tss_profile,"tss_binned_profile.tsv")

genes_gr <- suppressWarnings(genes(txdb))
sym <- as.data.table(ensembl_cached("entrez_to_symbol",
  attributes=c("entrezgene_id","hgnc_symbol"),filters="entrezgene_id",values=names(genes_gr)))
sym <- sym[nzchar(hgnc_symbol)]
sig[,entrez:=names(genes_gr)[nearest(mid,genes_gr)]]
sig <- merge(sig,sym,by.x="entrez",by.y="entrezgene_id",all.x=TRUE)

# Peaks nearest a set's member genes, counted by direction and tested against a
# 50:50 split. This asks whether accessibility shifts coherently with VPS72
# status, which the per-peak FDR does not answer.
geneset <- rbindlist(lapply(names(geneset_files),function(nm) {
  fp <- file.path(inputs$genesets,geneset_files[[nm]])
  if (!file.exists(fp)) { cat("  skipped, not found:",nm,"\n"); return(NULL) }
  raw <- readLines(fp,warn=FALSE)
  row <- grep("GENE_SYMBOLS",raw,value=TRUE)[1L]
  if (is.na(row)) stop("No GENE_SYMBOLS row in ",fp)   # MSigDB TSV export format
  members <- toupper(trimws(strsplit(sub("^[^\t]*\t","",row),",")[[1L]]))
  hit <- sig[toupper(hgnc_symbol) %chin% members]
  if (!nrow(hit)) return(NULL)
  n_up <- sum(hit$log2FC>0); n_tot <- nrow(hit)
  data.table(gene_set=nm,n_genes_in_set=length(members),n_peaks=n_tot,
             n_up_in_high=n_up,n_down_in_high=n_tot-n_up,frac_up=n_up/n_tot,
             median_log2FC=median(hit$log2FC),binom_p=binom.test(n_up,n_tot,0.5)$p.value)
}))
if (!length(geneset) || !nrow(geneset)) stop("No gene set matched any differential peak")
geneset[,binom_FDR:=p.adjust(binom_p,method="BH")]
setorder(geneset,binom_p)

write_tsv(geneset,"geneset_directionality.tsv")
print(geneset)
