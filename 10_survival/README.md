# VPS72 Pan-Cancer Survival Analysis

R scripts supporting the evaluation of VPS72 expression, copy-number alterations, and chromatin-remodeling signatures in relation to overall survival across human cancers.

These analyses accompany *Multi-omic analyses reveal VPS72 associations with chromatin remodeling, proliferation, and stress-response programs across human cancers*.

## Analysis Scripts

| Script | Analysis |
|---|---|
| `00_config.R` | Input configuration, survival endpoints, and analysis parameters |
| `01_tcga_overall_survival.R` | Pan-cancer and cancer-specific Cox proportional hazards analyses of VPS72 expression |
| `02_prognostic_discrimination.R` | Concordance index and time-dependent ROC AUC at 1, 3, and 5 years |
| `03_copy_number_survival.R` | Associations between VPS72 copy-number alterations and overall survival |
| `04_signature_survival.R` | Survival associations involving VPS72 and SRCAP/TIP60/INO80 chromatin-remodeling signatures |
| `05_glioma_adjusted_survival.R` | Glioma survival analyses adjusted for IDH status, MGMT methylation, transcriptional subtype, and clinical covariates |
| `06_pcawg_optional_survival.R` | Independent PCAWG survival analysis |

## Data Requirements

- **TCGA RNA-seq and clinical data:** UCSC Xena VPS72 expression and overall survival annotations.
- **Copy-number data:** TCGA PanCancer Atlas continuous copy-number values from cBioPortal.
- **Chromatin signatures:** TCGA/GTEx TPM expression matrices and corresponding clinical annotations.
- **Glioma datasets:** TCGA glioma expression and clinical annotations, including available IDH, MGMT, and molecular subtype information.
- **PCAWG:** Independent expression and survival data available through cBioPortal.

Required filenames and input directories are specified in `00_config.R`.

## Execution

Configure the input and output paths in `00_config.R`, then execute:

```r
source("01_tcga_overall_survival.R")
source("02_prognostic_discrimination.R")
source("03_copy_number_survival.R")
source("04_signature_survival.R")
source("05_glioma_adjusted_survival.R")
source("06_pcawg_optional_survival.R")
```

Script `02` requires the output of script `01`. Scripts `03`–`06` analyze separate datasets and may be executed independently.

## Statistical Methods

Overall survival associations are assessed using Cox proportional hazards models, with VPS72 expression evaluated as a continuous variable and through expression-based stratification.

Prognostic discrimination is summarized using the concordance index and time-dependent AUC. Copy-number associations are evaluated using continuous values and predefined alteration thresholds.

Glioma analyses include multivariable adjustment for available molecular and clinical covariates. Proportional hazards assumptions are evaluated using Schoenfeld residuals.

The PCAWG analysis is conducted independently of TCGA.

## Outputs

Scripts generate tab-delimited results containing survival model estimates, hazard ratios, confidence intervals, statistical significance, C-index and AUC values, copy-number associations, and molecularly adjusted glioma survival results.

**Requirements:** R ≥ 4.2; `data.table`, `survival`, `timeROC`, `dplyr`, `tidyr`, `readr`, and `tibble`. Optional restricted mean survival time analysis requires `survRM2`.