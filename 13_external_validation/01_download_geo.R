# Download the two Jostes et al. GEO series and unpack their supplementary
# files: GSE242227 (RNA-seq DESeq2 tables) and GSE246121 (H4ac ChIP bigWigs).
source("00_config.R")
suppressPackageStartupMessages({library(GEOquery);library(R.utils)})

dir.create(inputs$geo_dir,recursive=TRUE,showWarnings=FALSE)

fetch <- function(gse) {
  od <- file.path(inputs$geo_dir,gse)
  if (dir.exists(od) && length(list.files(od,recursive=TRUE))) {
    cat("  already present:",gse,"\n"); return(invisible(od))
  }
  dir.create(od,recursive=TRUE,showWarnings=FALSE)
  cat("  downloading:",gse,"\n")
  # The GEOquery argument name changed across versions.
  ok <- tryCatch({getGEOSuppFiles(GSE=gse,makeDirectory=FALSE,baseDir=od); TRUE},
                 error=function(e) tryCatch({getGEOSuppFiles(GEO=gse,makeDirectory=FALSE,baseDir=od); TRUE},
                                            error=function(e2) FALSE))
  if (!ok) stop("Download failed for ",gse,
                "\nFetch the supplementary files by hand into ",od)

  # RAW tarballs hold the bigWigs; .gz tables are expanded but kept.
  for (tf in list.files(od,pattern="\\.tar$",full.names=TRUE,recursive=TRUE))
    untar(tf,exdir=file.path(od,"untar"))
  gz <- setdiff(list.files(od,pattern="\\.gz$",full.names=TRUE,recursive=TRUE),
                list.files(od,pattern="\\.tar\\.gz$",full.names=TRUE,recursive=TRUE))
  for (f in gz) gunzip(f,remove=FALSE,overwrite=TRUE)
  invisible(od)
}

for (g in c(gse_rnaseq,gse_chip)) fetch(g)

# Report what landed, so a missing file is obvious before the next script.
found <- rbindlist(lapply(c(gse_rnaseq,gse_chip),function(g) {
  fs <- list.files(file.path(inputs$geo_dir,g),recursive=TRUE)
  if (!length(fs)) stop("Nothing downloaded for ",g)
  data.table(gse=g,file=fs)
}))
write_tsv(found,"geo_file_manifest.tsv")

missing_deg <- deg_files[!file.exists(file.path(inputs$geo_dir,gse_rnaseq,deg_files))]
if (length(missing_deg))
  cat("expected DEG table not at the top level of",gse_rnaseq,":",
      paste(missing_deg,collapse=", "),
      "\ncheck the manifest for the names GEO actually used\n")

cat("Files per series:\n"); print(found[,.N,by=gse])
