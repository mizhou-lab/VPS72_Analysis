# TCGA ATAC-seq Analysis of VPS72

R scripts supporting the analysis of chromatin accessibility associated with VPS72 expression across TCGA cancer types.

These analyses accompany *Multi-omic analyses reveal VPS72 associations with chromatin remodeling, proliferation, and stress-response programs across human cancers*.

## Analysis Scripts

| Script | Analysis |
|---|---|
| `00_config.R` | Input configuration, genomic annotations, analysis parameters, and gene sets |
| `01_prepare.R` | ATAC-seq and RNA-seq integration, sample matching, and VPS72-high/low stratification |
| `02_differential.R` | Genome-wide differential chromatin accessibility analysis using limma |
| `03_annotate_peaks.R` | Genomic annotation of differentially accessible peaks and identification of transcription-factor loci |
| `04_locus_and_profiles.R` | Cancer-specific locus accessibility, transcription start site profiles, and gene-set directionality analysis |

## Data Requirements

- **TCGA ATAC-seq:** Pan-cancer accessibility matrix, sample identifier mapping, and peak coordinates (Corces et al., *Science*, 2018).
- **TCGA RNA-seq:** UCSC Xena TOIL RSEM gene TPM expression matrix.
- **Gene sets:** MSigDB gene-set files specified in `00_config.R`.

Input datasets are available from TCGA-associated public repositories and UCSC Xena. Required filenames and directory settings are defined in `00_config.R`.

## Execution

Configure the input directories in `00_config.R`, then execute the scripts sequentially:

```r
source("01_prepare.R")
source("02_differential.R")
source("03_annotate_peaks.R")
source("04_locus_and_profiles.R")
```

## Statistical Analysis

VPS72-high and VPS72-low samples are defined by median VPS72 expression across the matched ATAC cohort. Differential accessibility is assessed using limma, with Benjamini–Hochberg adjustment for multiple testing.

Significant peaks are defined by FDR <escape><</escape> 0.05 and absolute log2 fold change > 0.25 for downstream annotation. The primary pan-cancer model does not adjust for cancer type (`adjust_for_cancer_type = FALSE`).

## Outputs

Scripts generate tab-delimited files containing matched sample classifications, genome-wide differential accessibility results, genomic peak annotations, cancer-specific locus summaries, TSS accessibility profiles, and gene-set enrichment directionality statistics.

**Requirements:** R ≥ 4.2; `data.table`, `matrixStats`, `limma`, `GenomicRanges`, `IRanges`, `TxDb.Hsapiens.UCSC.hg38.knownGene`, and `biomaRt`.