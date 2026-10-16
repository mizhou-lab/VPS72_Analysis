# VPS72 Prognostic Independence and 1q21 Adjustment

R scripts supporting the prognostic analyses of VPS72, MCL1, SETDB1 and regional 1q21 copy number in *Multi-omic analyses reveal VPS72 associations with chromatin remodeling, proliferation, and stress-response programs across human cancers*.

## Scripts

| Script | Analysis |
|---|---|
| `00_config.R` | TCGA and cBioPortal input filenames, 1q21 neighbor gene panel and modeling thresholds |
| `01_prognostic_independence.R` | Patient-level data matching, multivariable Cox models, copy-number strata, within-cancer models, collinearity and numerical KM source data |

## Inputs

- cBioPortal TCGA PanCancer Atlas: `mRNA expression z-scores relative to all samples (log RNA Seq V2 RSEM).txt` and `Log2 copy-number values.txt`.
- TCGA clinical export: `combined_study_clinical_data (4).tsv`, containing patient identifiers, sample type, cancer type, overall survival, age, sex, stage, aneuploidy and fraction genome altered.
- Optional `Scores_160_Signatures.tsv` with `Module11_Prolif_score` for the proliferation-adjusted model.

Set `CBIO_DIR`, `TCGA_DIR` and `OUT_DIR` in `00_config.R` or use the documented environment variables there. The copy-number and proliferation files are optional, but models requiring them will not be fitted without valid data.

## Run

From this folder:

```r
source("01_prognostic_independence.R")
```

## Statistical methods

Only TCGA primary tumors are retained, with one expression and clinical record per patient. Overall survival uses months; death = 1, censoring = 0. VPS72, MCL1, SETDB1 and regional copy-number variables are standardized within cancer type. Pooled Cox models use cancer-type-specific baseline hazards, progressively adjusting for age, sex, AJCC stage, aneuploidy, fraction genome altered, MCL1/SETDB1, and a 20-gene local 1q21 copy-number proxy. Per-cancer models adjust for age, sex, MCL1 and SETDB1. Comparator models use the same complete-case patients; likelihood-ratio P values and VIF are provided.

Continuous VPS72 copy-number gain is defined as log2 CN > 0.2; neutral is -0.2 to 0.2. The local 20-gene copy-number proxy is not chromosome-arm 1q gain. Per-cancer Cox P values receive BH correction. Cohorts with fewer than 50 complete cases or 20 events are not fitted.

## Outputs

TSV files under `results/`: `matched_patient_survival_data.tsv`, `table3_cox_models.tsv`, `cox_same_sample_comparators.tsv`, `cox_per_cancer.tsv`, `cohort_summary.tsv`, `expression_cnv_correlations.tsv`, `kaplan_meier_numerical_curves.tsv`, and `model_skips.tsv`; plus a session record.

**Requirements:** R >= 4.2; `data.table`, `survival`.
