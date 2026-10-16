# Per-cancer Hallmark GSEA and six-cancer pathway recurrence from existing DEG tables.
source("00_config.R")
if (!requireNamespace("fgsea", quietly = TRUE)) stop("Required package missing: fgsea")
if (!nzchar(HALLMARK_GMT) && !requireNamespace("msigdbr", quietly = TRUE))
  stop("Required package missing: msigdbr (or set a local Hallmark GMT)")
files <- unlist(DEG_FILES, use.names = TRUE)
missing_files <- files[!file.exists(files)]
if (length(missing_files)) stop("Missing DEG input files:\n", paste(missing_files, collapse = "\n"))

load_hallmarks <- function() {
  if (nzchar(HALLMARK_GMT)) {
    if (!file.exists(HALLMARK_GMT)) stop("Hallmark GMT not found: ", HALLMARK_GMT)
    return(fgsea::gmtPathways(HALLMARK_GMT))
  }
  warning("Using installed msigdbr Hallmark collection. Pin a GMT file to reproduce a specific MSigDB release.")
  sets <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
  split(sets$gene_symbol, sets$gs_name)
}
pathways <- load_hallmarks()

read_deg <- function(path) {
  raw <- data.table::fread(path, sep = "\t", check.names = FALSE)
  keys <- gsub("[^a-z0-9]", "", tolower(names(raw)))
  pick <- function(candidates, required = TRUE) {
    j <- match(candidates, keys, nomatch = 0L)
    j <- j[j > 0L]
    if (length(j)) return(names(raw)[j[1]])
    if (required) stop("Required DEG column not found in ", path, ": ", paste(candidates, collapse = "/"))
    NULL
  }
  g <- pick(c("gene", "genesymbol", "symbol", "genename"))
  l <- pick(c("log2fc", "log2foldchange", "avglog2fc", "avglogfc"))
  p <- pick(c("pval", "pvalue", "p", "pvalueunadjusted"), FALSE)
  q <- pick(c("qval", "qvalue", "padj", "fdr", "pvaladj", "adjustedpvalue"), FALSE)
  if (is.null(q)) stop("Adjusted P/FDR column required in ", path)
  x <- data.table(Gene = trimws(as.character(raw[[g]])),
    log2FC = suppressWarnings(as.numeric(raw[[l]])),
    pval = if (is.null(p)) NA_real_ else suppressWarnings(as.numeric(raw[[p]])),
    qval = suppressWarnings(as.numeric(raw[[q]])))
  x <- x[!is.na(Gene) & nzchar(Gene) & is.finite(log2FC)]
  if (any(x$qval[is.finite(x$qval)] < 0 | x$qval[is.finite(x$qval)] > 1))
    stop("FDR values outside [0,1] in ", path)
  # Deterministic one-row-per-symbol mapping: lowest FDR, then largest |log2FC|.
  x[, abs_fc := abs(log2FC)]
  setorder(x, Gene, qval, -abs_fc)
  x <- x[!duplicated(Gene)]
  x[, abs_fc := NULL]
  x
}

all_deg <- list(); all_fg <- list(); summaries <- list(); manifest <- list()
for (sig in SIGNATURES) for (ct in CANCERS) {
  fp <- DEG_FILES[[sig]][[ct]]
  d <- read_deg(fp)
  if (nrow(d) == 0L) stop("No usable genes in ", fp)
  d[, `:=`(Signature = sig, Cancer = ct)]
  all_deg[[paste(sig, ct, sep = "_")]] <- d
  summaries[[paste(sig, ct, sep = "_")]] <- data.table(Signature = sig, Cancer = ct,
    Genes_tested = nrow(d), Significant = sum(d$qval < FDR_CUT & abs(d$log2FC) > DEG_LFC_CUT, na.rm = TRUE),
    Up = sum(d$qval < FDR_CUT & d$log2FC > DEG_LFC_CUT, na.rm = TRUE),
    Down = sum(d$qval < FDR_CUT & d$log2FC < -DEG_LFC_CUT, na.rm = TRUE))
  manifest[[paste(sig, ct, sep = "_")]] <- data.table(Signature = sig, Cancer = ct,
    File = fp, Genes_used = nrow(d))
  ranks <- setNames(d$log2FC, d$Gene)
  ranks <- sort(ranks, decreasing = TRUE)
  set.seed(FGSEA_SEED)
  fg <- as.data.table(fgsea::fgsea(pathways = pathways, stats = ranks,
    nperm = FGSEA_PERMUTATIONS, minSize = FGSEA_MIN_SIZE, maxSize = FGSEA_MAX_SIZE))
  if (nrow(fg)) {
    fg[, leadingEdge := vapply(leadingEdge, function(z) paste(z, collapse = ";"), "")]
    fg[, `:=`(Signature = sig, Cancer = ct, Direction = fifelse(NES > 0, "Up", "Down"))]
    setorder(fg, padj, -NES)
  }
  all_fg[[paste(sig, ct, sep = "_")]] <- fg
  message(sig, " / ", ct, ": ", nrow(d), " genes, ", nrow(fg), " Hallmark results")
}
deg <- rbindlist(all_deg, fill = TRUE)
fg <- rbindlist(all_fg, fill = TRUE)
if (anyDuplicated(fg[, .(Signature, Cancer, pathway)])) stop("Duplicate signature/cancer/pathway fgsea rows")
write_result(deg, "deg_standardized.tsv")
write_result(rbindlist(summaries), "deg_summary.tsv")
write_result(rbindlist(manifest), "deg_input_manifest.tsv")
write_result(fg, "hallmark_fgsea_all.tsv")
rec <- fg[is.finite(padj) & padj < FDR_CUT,
  .(Cancer_Count = uniqueN(Cancer), Cancers = paste(sort(unique(Cancer)), collapse = ";"),
    Mean_NES = mean(NES, na.rm = TRUE), Mean_minusLog10FDR = mean(-log10(pmax(padj, .Machine$double.xmin)))),
  by = .(Signature, pathway, Direction)]
setorder(rec, Signature, Direction, -Cancer_Count, pathway)
write_result(rec, "hallmark_recurrence.tsv")
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
