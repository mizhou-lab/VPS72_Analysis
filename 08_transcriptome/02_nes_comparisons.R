# Numerical pathway-level NES comparisons; no plotting or Excel workbooks.
source("00_config.R")
fp <- file.path(OUT_DIR, "hallmark_fgsea_all.tsv")
if (!file.exists(fp)) stop("Run 01_hallmark_fgsea.R first. Missing: ", fp)
f <- fread(fp)
need <- c("Signature", "Cancer", "pathway", "NES", "padj")
if (!all(need %in% names(f))) stop("Unexpected fgsea output schema")
if (anyDuplicated(f[, .(Signature, Cancer, pathway)])) stop("Duplicate pathway keys")

# Preserve cancer in the join key. Do not replace absent NES values with zero.
wide <- dcast(f[, ..need], Cancer + pathway ~ Signature, value.var = c("NES", "padj"))
for (nm in c("NES_VPS72", "NES_SRCAP", "NES_SRCAP_resid", "padj_VPS72", "padj_SRCAP", "padj_SRCAP_resid"))
  if (!nm %in% names(wide)) wide[, (nm) := NA_real_]
wide[, `:=`(
  Delta_NES_SRCAP_minus_VPS72 = NES_SRCAP - NES_VPS72,
  Delta_NES_SRCAP_minus_resid = NES_SRCAP - NES_SRCAP_resid,
  Concordant_SRCAP_VPS72 = fifelse(is.finite(NES_SRCAP) & is.finite(NES_VPS72) &
    NES_SRCAP != 0 & NES_VPS72 != 0, sign(NES_SRCAP) == sign(NES_VPS72), NA),
  Both_FDR05_SRCAP_VPS72 = !is.na(padj_SRCAP) & !is.na(padj_VPS72) &
    padj_SRCAP < FDR_CUT & padj_VPS72 < FDR_CUT,
  Both_FDR05_SRCAP_resid = !is.na(padj_SRCAP) & !is.na(padj_SRCAP_resid) &
    padj_SRCAP < FDR_CUT & padj_SRCAP_resid < FDR_CUT
)]
setorder(wide, Cancer, pathway)
write_result(wide, "nes_comparisons_by_cancer.tsv")
write_result(wide[Cancer == "BRCA"], "nes_comparisons_BRCA.tsv")

# Descriptive equal-cancer-weighted NES means, matching the 6-cohort plot inputs.
means <- f[, .(Mean_NES = mean(NES, na.rm = TRUE),
    Cancers_measured = uniqueN(Cancer), FDR05_Cancer_Count = sum(padj < FDR_CUT, na.rm = TRUE)),
  by = .(Signature, pathway)]
setorder(means, Signature, pathway)
write_result(means, "nes_means_six_cancers.tsv")
pan <- dcast(means, pathway ~ Signature, value.var = c("Mean_NES", "Cancers_measured", "FDR05_Cancer_Count"))
for (nm in c("Mean_NES_SRCAP", "Mean_NES_VPS72", "Mean_NES_SRCAP_resid"))
  if (!nm %in% names(pan)) pan[, (nm) := NA_real_]
pan[, `:=`(Delta_SRCAP_minus_VPS72 = Mean_NES_SRCAP - Mean_NES_VPS72,
            Delta_SRCAP_minus_resid = Mean_NES_SRCAP - Mean_NES_SRCAP_resid)]
setorder(pan, pathway)
write_result(pan, "nes_pancancer_comparisons.tsv")
