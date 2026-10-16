# VPS72 Tumor–Normal Expression and Proliferation Analysis

R scripts examining associations between VPS72 expression and transcriptional programs in TCGA tumors and adjacent normal tissues, including cancer-specific comparisons and pan-cancer meta-analysis.

These analyses accompany *Multi-omic analyses reveal VPS72 associations with chromatin remodeling, proliferation, and stress-response programs across human cancers*.

## Analysis Scripts

| Script | Analysis |
|---|---|
| `00_config.R` | Dataset configuration, Hallmark gene sets, and statistical parameters |
| `01_geo_normals.R` | Retrieval and processing of independent normal epithelial RNA-seq datasets |
| `02_tcga_tumor_normal.R` | TCGA tumor and adjacent-normal sample selection and VPS72 expression |
| `03_master_table.R` | Integration of sample metadata and normalized expression |
| `04_pathway_scores.R` | Hallmark pathway scoring and VPS72-high versus VPS72-low comparisons |
| `05_continuous_meta.R` | Cancer-specific Spearman correlations, random-effects meta-analysis, and tumor–normal interaction tests |
| `06_summary_tables.R` | Cancer-specific and pooled statistical summaries |

## Data Requirements

- **TCGA RNA-seq:** UCSC Xena TOIL RSEM gene TPM expression matrix and corresponding phenotype annotations.
- **Normal epithelial datasets:** GEO accessions GSE190615, GSE117271, GSE161892, and GSE227502.
- **Gene sets:** MSigDB Hallmark collection, as specified in `00_config.R`.

Input paths, cancer types, and analysis parameters are defined in `00_config.R`.

## Execution

Configure the input and output directories in `00_config.R`, then execute sequentially:

```r
source("01_geo_normals.R")
source("02_tcga_tumor_normal.R")
source("03_master_table.R")
source("04_pathway_scores.R")
source("05_continuous_meta.R")
source("06_summary_tables.R")
```

The GEO download step requires an internet connection. Subsequent analyses use the generated expression and metadata tables.

## Statistical Methods

Within each cancer type, tumor and normal samples are independently stratified by median VPS72 expression. Hallmark pathway scores are compared between VPS72-high and VPS72-low groups using two-sided Wilcoxon rank-sum tests.

Continuous associations are evaluated using Spearman correlations, followed by random-effects meta-analysis. Tumor-versus-normal differences are assessed through meta-regression and interaction models.

The proliferation composite includes E2F Targets, G2M Checkpoint, MYC Targets V1/V2, and Mitotic Spindle. Multiple-testing correction uses the Benjamini–Hochberg procedure.

## Outputs

Scripts generate tab-delimited tables containing sample-level pathway scores, high-versus-low comparisons, cancer-specific correlations, pooled meta-analysis estimates, tumor–normal interaction statistics, and summary results.

**Requirements:** R ≥ 4.2; `data.table`, `metafor`, `GEOquery`, with `lme4` and `lmerTest` for optional paired-patient sensitivity analyses.