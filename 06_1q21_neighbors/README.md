# VPS72 and Neighboring 1q21 Genes

R scripts supporting the pan-cancer analysis of VPS72 and 20 neighboring genes at chromosome 1q21, including copy-number alterations, expression dosage, regional regulation, and clinical associations.

These analyses accompany *Multi-omic analyses reveal VPS72 associations with chromatin remodeling, proliferation, and stress-response programs across human cancers*.

## Analysis Scripts

| Script | Analysis |
|---|---|
| `00_config.R` | Input configuration, 21-gene panel, and statistical parameters |
| `01_gain_and_dosage.R` | Copy-number gain and amplification frequencies, CNV–expression associations, and dosage responsiveness |
| `02_heldout_and_copyneutral.R` | Copy-neutral reference models and expression associations in held-out gained tumors |
| `03_regional_adjustment.R` | VPS72 expression adjusted for neighboring 1q21 genes and regional copy-number alterations |
| `04_optional_clinical_covariates.R` | Optional survival and proliferation models adjusted for neighboring genes and genomic covariates |

## Data Requirements

Analyses use TCGA PanCancer Atlas datasets obtained through cBioPortal, including:

- GISTIC2 copy-number alteration calls and continuous log2 copy-number values.
- RNA-seq expression z-scores relative to diploid, all-tumor, and, where available, normal-reference samples.
- Batch-normalized RSEM expression data for regional adjustment.
- Optional TCGA clinical survival and proliferation data.

Required filenames and directories are specified in `00_config.R`.

## Execution

Configure the input and output paths in `00_config.R`, then run:

```r
source("01_gain_and_dosage.R")
source("02_heldout_and_copyneutral.R")
source("03_regional_adjustment.R")
```

Optional clinical analyses can be run separately:

```r
source("04_optional_clinical_covariates.R")
```

Scripts 01–02 share the cBioPortal expression z-score workflow. Script 03 independently analyzes raw RSEM expression; script 04 uses the regional analysis output.

## Statistical Methods

Copy-number gains are defined as GISTIC2 +1 or +2, amplifications as +2, and copy-neutral samples as 0.

Within-cancer analyses evaluate CNV–expression associations, gained-versus-neutral expression differences, and dosage responsiveness. Regional models assess VPS72 expression after accounting for neighboring 1q21 genes and local copy-number variation.

Statistical tests include Spearman correlations, Wilcoxon tests, multivariable regression, and optional Cox proportional hazards models. Multiple-testing correction uses the Benjamini–Hochberg procedure.

## Outputs

Scripts generate tab-delimited tables and RDS files containing copy-number frequencies, dosage associations, expression effect sizes, regional regression estimates, and optional clinical associations.

**Requirements:** R ≥ 4.1; `data.table` and `survival` for optional clinical analyses.