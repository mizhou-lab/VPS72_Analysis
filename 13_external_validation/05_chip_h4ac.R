# Promoter H4 acetylation around the TSS of signature genes, control versus
# YL1 knockdown. The bigWigs are hg19 (chip_genome in 00_config.R).
source("00_config.R")
suppressPackageStartupMessages({
  library(rtracklayer);library(GenomicRanges);library(GenomeInfoDb)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene);library(org.Hs.eg.db);library(AnnotationDbi)
})

bw_files <- list.files(file.path(inputs$geo_dir,gse_chip),pattern="\\.bw$",
                       full.names=TRUE,recursive=TRUE)
bw_files <- bw_files[grepl("h4ac",tolower(basename(bw_files)))]
if (!length(bw_files)) stop("No H4ac bigWig found under ",file.path(inputs$geo_dir,gse_chip),
                            "\nRun 01 first; the bigWigs arrive inside the RAW tarball.")

# Conditions are matched on the file-name tokens GEO ships; the two control
# files differ only by a hyphen versus an underscore, so order matters here.
assign_condition <- function(f) {
  lf <- tolower(basename(f))
  hit <- names(chip_conditions)[vapply(chip_conditions,grepl,logical(1),x=lf,fixed=TRUE)]
  if (!length(hit)) NA_character_ else hit[1L]
}
meta <- data.table(file=bw_files,filename=basename(bw_files))
meta[,condition:=vapply(file,assign_condition,character(1))]
if (anyNA(meta$condition)) stop("Unmatched bigWig: ",
  paste(meta[is.na(condition)]$filename,collapse=", "),
  "\nAdd its token to chip_conditions in 00_config.R")
meta <- meta[,.SD[1L],by=condition]
if (!all(c(chip_controls,chip_kd) %chin% meta$condition))
  stop("Missing condition: ",paste(setdiff(c(chip_controls,chip_kd),meta$condition),collapse=", "))
write_tsv(meta[,.(condition,filename)],"h4ac_bigwig_samples.tsv")

bw <- setNames(meta$file,meta$condition)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
gene_gr <- genes(txdb)

# One TSS window per symbol, from the first gene model Entrez returns.
sig_genes <- unique(unlist(load_signatures(),use.names=FALSE))
map <- suppressMessages(as.data.table(
  AnnotationDbi::select(org.Hs.eg.db,keys=sig_genes,keytype="SYMBOL",columns="ENTREZID")))
map <- map[!is.na(ENTREZID)][,.SD[1L],by=SYMBOL]
map <- map[ENTREZID %chin% names(gene_gr)]
if (!nrow(map)) stop("No signature gene mapped to an hg19 gene model")

g <- gene_gr[map$ENTREZID]
tss <- ifelse(as.character(strand(g))=="-",end(g),start(g))
win <- GRanges(seqnames(g),IRanges(pmax(1L,tss-tss_flank),tss+tss_flank))
seqlevelsStyle(win) <- "UCSC"
mcols(win)$gene <- map$SYMBOL
cat("promoter windows:",length(win),"of",length(sig_genes),"signature genes\n")

# Mean coverage in the window, read per bigWig in one pass over all windows.
mean_in_windows <- function(f) {
  bwf <- BigWigFile(f)
  common <- intersect(seqlevels(win),seqlevels(seqinfo(bwf)))
  if (!length(common)) stop("No shared chromosome names between ",basename(f),
                            " and the hg19 gene models")
  w <- win[as.character(seqnames(win)) %chin% common]
  cov <- rtracklayer::import(bwf,which=w,as="NumericList")
  data.table(gene=mcols(w)$gene,value=vapply(cov,function(v) {
    v <- v[is.finite(v)]; if (!length(v)) NA_real_ else mean(v)
  },numeric(1)))
}

vals <- rbindlist(lapply(names(bw),function(nm) {
  d <- mean_in_windows(bw[[nm]]); d[,condition:=nm]; d[]
}))
wide <- dcast(vals,gene~condition,value.var="value")

ctrl_ok <- intersect(chip_controls,names(wide)); kd_ok <- intersect(chip_kd,names(wide))
wide[,control_mean:=rowMeans(.SD,na.rm=TRUE),.SDcols=ctrl_ok]
wide[,kd_mean:=rowMeans(.SD,na.rm=TRUE),.SDcols=kd_ok]
wide <- wide[is.finite(control_mean) & is.finite(kd_mean)]
if (!nrow(wide)) stop("No gene had coverage in both arms")

# Delta is the quantity ranked; the log2 ratio is reported alongside it and
# carries a pseudocount, so it is unstable where coverage is near zero.
wide[,`:=`(delta=kd_mean-control_mean,
           log2_ratio=log2((kd_mean+1e-6)/(control_mean+1e-6)))]
setorder(wide,delta)
wide[,rank:=.I]

# Which signature each gene belongs to, so the ranking can be read per programme.
sigs <- load_signatures()
memb <- rbindlist(lapply(names(sigs),function(nm)
  data.table(signature=nm,gene=sigs[[nm]])))[gene %chin% wide$gene]
wide <- merge(wide,memb[,.(signatures=paste(sort(unique(signature)),collapse=";")),by=gene],
              by="gene",all.x=TRUE)
setorder(wide,rank)

write_tsv(wide,"promoter_h4ac_ranking.tsv")
cat("Genes ranked:",nrow(wide),"; losing H4ac on knockdown:",sum(wide$delta<0),"\n")
print(head(wide[,.(rank,gene,control_mean=round(control_mean,3),
                   kd_mean=round(kd_mean,3),delta=round(delta,3))],10))
