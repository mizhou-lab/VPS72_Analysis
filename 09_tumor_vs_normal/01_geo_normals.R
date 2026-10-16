# Download the four normal-tissue bulk RNA-seq series and harmonize them into
# one expression table. Treated arms and tumour-derived samples are dropped here.
source("00_config.R")
suppressPackageStartupMessages({library(GEOquery)})

dir.create(inputs$geo_dir,recursive=TRUE,showWarnings=FALSE)

# Supplementary files only; the series matrices carry no counts for these.
for (gse in names(geo_series)) {
  dest <- file.path(inputs$geo_dir,gse)
  if (dir.exists(dest) && length(list.files(dest))) { cat("  already present:",gse,"\n"); next }
  cat("  downloading:",gse,"-",geo_series[[gse]],"\n")
  getGEOSuppFiles(gse,baseDir=inputs$geo_dir,makeDirectory=TRUE)
}

# Each series ships its counts under its own column naming, so the matrix is
# located by extension and the gene column by position.
read_series <- function(gse) {
  fs <- list.files(file.path(inputs$geo_dir,gse),full.names=TRUE,
                   pattern="\\.(txt|tsv|csv)(\\.gz)?$")
  fs <- fs[!grepl("filelist|README",fs,ignore.case=TRUE)]
  if (!length(fs)) stop("No count matrix found for ",gse,"; inspect ",file.path(inputs$geo_dir,gse))
  m <- fread(fs[1L])
  setnames(m,1,"gene")
  m[,gene:=toupper(sub("\\..*$","",as.character(gene)))]   # strip Ensembl version suffixes
  long <- melt(m,id.vars="gene",variable.name="sample",value.name="expr")
  long[,`:=`(dataset=gse,expr=suppressWarnings(as.numeric(expr)))]
  long[is.finite(expr)]
}

# Sample-level exclusions, stated per series rather than applied silently.
keep_samples <- function(x) {
  switch(unique(x$dataset),
    # Primary-cell preparations, not donors; all are normal.
    GSE190615 = x,
    GSE117271 = x,
    GSE161892 = x[!grepl("tumor|tumour|cancer|BRCA1|mutant",sample,ignore.case=TRUE)],
    # Hypoxia-treated samples are excluded from the normal reference.
    GSE227502 = x[grepl("^N0[1-6]|normox",sample,ignore.case=TRUE)],
    x)
}

parts <- lapply(names(geo_series),function(g) keep_samples(read_series(g)))
geo <- rbindlist(parts,fill=TRUE)
if (!nrow(geo)) stop("No GEO samples survived the exclusions")

# Standardized within each series, because library preparation and platform
# differ between them and a pooled z-score would carry that batch structure.
geo[,expr_z:=as.numeric(scale(expr)),by=.(dataset,gene)]

meta <- unique(geo[,.(dataset,sample,source="GEO",source_type="normal_reference",
                      state="Normal",description=geo_series[dataset])])
write_tsv(meta,"geo_normal_metadata.tsv")

fp <- file.path(OUT_DIR,"geo_normal_expression.tsv.gz")
fwrite(geo[,.(dataset,sample,gene,expr,expr_z)],fp,sep="\t")
cat("wrote",fp,"(",nrow(geo),"rows )\n")
cat("Series:",uniqueN(geo$dataset),"; samples:",uniqueN(geo$sample),
    "; genes:",uniqueN(geo$gene),"\n")
print(geo[,.(samples=uniqueN(sample)),by=dataset])
