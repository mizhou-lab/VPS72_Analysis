# Genomic category of each differential peak, and the subset at the
# transcription-factor loci in tf_genes.
source("00_config.R")
suppressPackageStartupMessages({library(TxDb.Hsapiens.UCSC.hg38.knownGene)})

da <- fread(require_input(inputs$da))
sig <- da[padj<fdr_cut & abs(log2FC)>lfc_annot]
if (!nrow(sig)) stop("No peak passes FDR < ",fdr_cut," and |log2FC| > ",lfc_annot)
cat("peaks carried forward:",nrow(sig),"\n")

gr <- GRanges(sig$chrom,IRanges(sig$start,sig$end))
seqlevelsStyle(gr) <- "UCSC"

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
prom <- suppressWarnings(promoters(genes(txdb),upstream=2000,downstream=2000))
exon <- exons(txdb)
intr <- unlist(intronsByTranscript(txdb))

# Precedence promoter > exon > intron > intergenic, so each peak carries
# exactly one label.
sig[,genomic_category:="intergenic"]
sig[overlapsAny(gr,intr),genomic_category:="intron"]
sig[overlapsAny(gr,exon),genomic_category:="exon"]
sig[overlapsAny(gr,prom),genomic_category:="promoter"]

write_tsv(sig,"peaks_annotated.tsv")
print(sig[,.N,by=.(direction,genomic_category)][order(direction,-N)])

# Gene bodies plus flank_bp, to catch upstream regulatory peaks.
ann <- as.data.table(ensembl_cached("tf_gene_coords",
  attributes=c("hgnc_symbol","chromosome_name","start_position","end_position"),
  filters="hgnc_symbol",values=tf_genes))
ann <- ann[chromosome_name %chin% as.character(c(1:22,"X","Y")) & nzchar(hgnc_symbol)]
if (!nrow(ann)) stop("No TF gene coordinates on a primary chromosome")

tf_gr <- GRanges(paste0("chr",ann$chromosome_name),
                 IRanges(pmax(1L,ann$start_position-flank_bp),ann$end_position+flank_bp),
                 gene=ann$hgnc_symbol)
hits <- findOverlaps(gr,tf_gr)
if (!length(hits)) stop("No differential peak falls at a TF locus")
tf_peaks <- cbind(sig[queryHits(hits)],data.table(tf_gene=mcols(tf_gr)$gene[subjectHits(hits)]))

write_tsv(tf_peaks,"peaks_tf_loci.tsv")
print(tf_peaks[,.N,by=.(tf_gene,direction)][order(tf_gene)])
