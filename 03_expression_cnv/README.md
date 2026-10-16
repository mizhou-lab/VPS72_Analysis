# VPS72 Expression and Copy-Number Analysis

R scripts supporting the pan-cancer analysis of VPS72 expression, copy-number alterations, chromatin-remodeling signatures, and clinical associations.

These analyses accompany *Multi-omic analyses reveal VPS72 associations with chromatin remodeling, proliferation, and stress-response programs across human cancers*.

## Analysis Scripts

| Script | Analysis |
|---|---|
| `00_config.R` | Input configuration, gene signatures, and analysis parameters |
| `01_primary_tumor_signatures.R` | TCGA expression and CNV integration, signature scoring, and VPS72-associated pathway analysis |
| `02_tumor_normal_expression.R` | Tumor-versus-normal expression comparisons and chromatin-complex signatures |
| `03_copy_number_dosage.R` | Copy-number alteration frequencies and CNV–expression associations |
| `04_vps72_residual_models.R` | VPS72 expression residuals adjusted for SRCAP/TIP60 signatures and copy number |
| `05_stage_proliferation.R` | Associations of VPS72 with clinical stage and proliferation |

## Data Requirements

Analyses use publicly available TCGA and UCSC Xena expression, copy-number, clinical, and pathway-signature datasets. Required input filenames are specified in `00_config.R`.

## Execution

Configure the input and output directories in `00_config.R`, then run the relevant scripts in R.

```r
source("01_primary_tumor_signatures.R")
source("02_tumor_normal_expression.R")
source("03_copy_number_dosage.R")
source("04_vps72_residual_models.R")
source("05_stage_proliferation.R")
```

Run script `01` before `04`, as the residual analysis uses its generated sample-level matrices. Other analyses can be executed independently.

**Requirements:** R ≥ 4.1, `data.table`, `matrixStats`, and, where Ensembl identifier mapping is required, `AnnotationDbi` and `org.Hs.eg.db`.

## Outputs

Scripts generate tab-delimited results and RDS files containing sample-level data, gene-signature scores, statistical comparisons, regression coefficients, and cancer-specific summaries.

The analyses use distinct source datasets and expression normalization methods, as specified in the individual scripts.