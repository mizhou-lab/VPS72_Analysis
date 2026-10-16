# VPS72 DepMap Dependency and Expression Analysis

R scripts examining VPS72 genetic dependency, expression, copy-number alterations, genomic instability, and associated transcriptional programs across cancer cell lines.

These analyses accompany *Multi-omic analyses reveal VPS72 associations with chromatin remodeling, proliferation, and stress-response programs across human cancers*.

## Analysis Scripts

| Script | Analysis |
|---|---|
| `00_config.R` | Input configuration, gene sets, dependency thresholds, and analysis parameters |
| `01_dependency.R` | VPS72 CRISPR and RNAi dependency, lineage-specific summaries, and SRCAP co-dependency |
| `02_expression_cnv.R` | VPS72 expression in cancer and non-cancerous cell lines and copy-number–expression associations |
| `03_genomic_instability.R` | Associations with chromosomal instability, aneuploidy, ploidy, loss of heterozygosity, and microsatellite instability |
| `04_codependency.R` | Identification of VPS72 co-dependent genes across CRISPR and RNAi screens |
| `05_hallmark_programs.R` | Associations between VPS72 expression, copy number, and MSigDB Hallmark transcriptional programs |

## Data Requirements

Analyses use the **DepMap Public 25Q2 release**, including:

- CRISPR Chronos and RNAi DEMETER2 dependency data.
- VPS72 expression and gene-level copy-number data.
- Cell-line metadata and genomic instability signatures.
- Genome-wide gene expression matrices and MSigDB Hallmark gene sets.
- VPS72 co-dependency datasets from CRISPR and RNAi screens.

Data are available through the [DepMap Portal](https://depmap.org/portal/). Required filenames and input directories are specified in `00_config.R`.

## Execution

Configure the input and output directories in `00_config.R`, then run:

```r
source("01_dependency.R")
source("02_expression_cnv.R")
source("03_genomic_instability.R")
source("04_codependency.R")
source("05_hallmark_programs.R")
```

Each analysis can be executed independently.

## Statistical Methods

VPS72 dependency is evaluated using CRISPR Chronos and RNAi DEMETER2 scores, with thresholds of ≤ −0.5 and ≤ −1.0 reported separately.

Associations between VPS72 expression, copy number, genomic instability, co-dependency, and transcriptional programs are evaluated using correlation analyses and group comparisons.

Hallmark pathway scores are calculated from standardized gene expression, with cancer-lineage-specific analyses performed where sufficient cell lines are available. Multiple-testing correction uses the Benjamini–Hochberg procedure.

## Outputs

Scripts generate tab-delimited results containing dependency distributions, lineage-specific summaries, expression and copy-number associations, genomic instability correlations, co-dependency statistics, and Hallmark pathway associations.

**Requirements:** R ≥ 4.2; `data.table` and `msigdbr`.