# VPS72 Transcriptome and Hallmark GSEA

R scripts supporting the VPS72, SRCAP and SRCAP-residual transcriptomic analyses in *Multi-omic analyses reveal VPS72 associations with chromatin remodeling, proliferation, and stress-response programs across human cancers*.

## Scripts

| Script | Analysis |
|---|---|
| `00_config.R` | Input file mapping, cohorts, gene-set source and thresholds |
| `01_hallmark_fgsea.R` | Standardizes existing DEG tables, runs Hallmark fgsea, and summarizes significant pathway recurrence across cancers |
| `02_nes_comparisons.R` | Per-cancer, BRCA and equal-cancer-weighted comparisons of VPS72, SRCAP and SRCAP-residual normalized enrichment scores (NES) |

## Inputs

Precomputed gene-level differential-expression TSV files (gene, log2 fold change, and adjusted P value) for VPS72, SRCAP and SRCAP residual in BRCA, CRC, GBM, HNSC, LUAD and LUSC. Exact filenames are defined in `00_config.R`. The scripts analyze these upstream DEG results; they **do not** recalculate differential expression from RNA-seq counts.

Hallmark gene sets are loaded from a versioned MSigDB GMT specified via `VPS72_HALLMARK_GMT`; if unset, the installed `msigdbr` Hallmark collection is used. Pin a GMT/release for stable reproduction.

## Run

From this folder, set input/output paths in `00_config.R`, then run:

```r
source("01_hallmark_fgsea.R")
source("02_nes_comparisons.R")
```

## Methods and outputs

Hallmark fgsea uses each gene's unflipped log2 fold change as the ranking statistic (10,000 permutations; seed 2026). Duplicate gene symbols are resolved by smallest adjusted P value, then largest absolute fold change. Significant pathways use BH-adjusted P < 0.05. Recurrence is counted across six cancer types for each signature and enrichment direction; pan-cancer mean NES weights cancer types equally. NES comparison tables join by both cancer and pathway, without treating missing values as zero.

Outputs in `results/`: `deg_standardized.tsv`, `deg_summary.tsv`, `hallmark_fgsea_all.tsv`, `hallmark_recurrence.tsv`, `nes_comparisons_by_cancer.tsv`, `nes_comparisons_BRCA.tsv`, `nes_means_six_cancers.tsv`, and `nes_pancancer_comparisons.tsv`, plus an input manifest and R session record.

**Requirements:** R >= 4.2; `data.table`, `fgsea`, `msigdbr` (the last is optional when using a local GMT).
