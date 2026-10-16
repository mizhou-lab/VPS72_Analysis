# Load the three author-provided DESeq2 tables and put them on one schema with
# a single significance rule: padj < padj_cut and |log2FC| > lfc_cut.
source("00_config.R")

# The gene identifier sits in an unnamed first column whose header differs
# between the three files, so it is taken by position rather than by name.
read_deg <- function(contrast,file) {
  d <- fread(require_input(file.path(inputs$geo_dir,gse_rnaseq,file)))
  if (!nrow(d)) stop("Empty DEG table: ",file)
  if (names(d)[1L] %chin% c("","V1")) setnames(d,1,"gene")
  gcol <- if ("gene" %chin% names(d)) "gene" else names(d)[1L]
  out <- data.table(
    contrast  = contrast,
    gene      = toupper(as.character(d[[gcol]])),
    base_mean = suppressWarnings(as.numeric(d[[find_col(d,c("baseMean","base_mean"))]])),
    log2fc    = suppressWarnings(as.numeric(d[[find_col(d,c("log2FoldChange","log2fold_change","log2fc"))]])),
    lfc_se    = suppressWarnings(as.numeric(d[[find_col(d,c("lfcSE","lfc_se"))]])),
    stat      = suppressWarnings(as.numeric(d[[find_col(d,c("stat"))]])),
    pvalue    = suppressWarnings(as.numeric(d[[find_col(d,c("pvalue","p_value"))]])),
    padj      = suppressWarnings(as.numeric(d[[find_col(d,c("padj","adj_p","FDR"))]])))
  out[nzchar(gene) & gene!="NA"]
}

deg <- rbindlist(Map(read_deg,names(deg_files),deg_files))

# A gene tested more than once in a table keeps its most significant row.
setorder(deg,contrast,gene,padj)
deg <- deg[,.SD[1L],by=.(contrast,gene)]

deg[,sig:=!is.na(padj) & padj<padj_cut & !is.na(log2fc) & abs(log2fc)>lfc_cut]
deg[,direction:=fcase(sig & log2fc>0,"Up",
                      sig & log2fc<0,"Down",
                      default="NS")]

write_tsv(deg,"deg_standardized.tsv")

summ <- deg[,.(n_tested=.N,n_sig=sum(sig),
               n_up=sum(direction=="Up"),n_down=sum(direction=="Down"),
               median_abs_lfc_sig=median(abs(log2fc[sig]),na.rm=TRUE)),by=contrast]
write_tsv(summ,"deg_summary.tsv")
print(summ)
cat("Significance rule: padj <",padj_cut,"and |log2FC| >",lfc_cut,"\n")
