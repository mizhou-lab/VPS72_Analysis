# VPS72 Genomic Association Analysis

R scripts supporting the pan-cancer analysis of VPS72-associated somatic mutations, chromatin-remodeling complexes, and biological pathways.

These analyses accompany *Multi-omic analyses reveal VPS72 associations with chromatin remodeling, proliferation, and stress-response programs across human cancers*.

## Analysis Scripts

| Script | Analysis |
|---|---|
| `00_config.R` | Input configuration, gene sets, clinical covariates, and statistical parameters |
| `01_master_table.R` | Integration of TCGA expression, clinical information, and somatic mutation data |
| `02_complex_matrices.R` | Expression, copy-number, and amplification summaries for SRCAP/TIP60/INO80 complex genes |
| `03_disruption_score.R` | Composite epigenetic disruption scores and gene ranking |
| `04_vps72_residual.R` | VPS72 expression residuals adjusted for chromatin-complex signatures and cancer type |
| `05_mutation_association.R` | Associations between VPS72 expression and somatic mutation status |
| `06_discrimination_auc.R` | ROC and precision–recall analyses of mutation-associated genes |
| `07_tp53_variants.R` | Associations between recurrent TP53 protein variants and VPS72 expression |
| `08_comutation_fisher.R` | Driver-mutation enrichment in VPS72-high versus VPS72-low tumors |
| `09_pathway_enrichment.R` | GO Biological Process, KEGG, and Hallmark enrichment analyses |

## Data Requirements

- **TCGA RNA-seq:** EB++AdjustPANCAN Illumina HiSeq RNASeqV2 expression dataset.
- **Somatic mutations:** TCGA MC3 mutation calls (v0.2.8).
- **Copy-number alterations:** TCGA pan-cancer GISTIC2 dataset.
- **Clinical annotations:** TCGA patient-level demographic and cancer-type information.
- **Independent cohorts:** cBioPortal expression and mutation datasets for cohort-specific analyses.

Input directories, gene sets, and thresholds are specified in `00_config.R`.

## Execution

Configure the input paths in `00_config.R`, then execute:

```r
source("01_master_table.R")
source("02_complex_matrices.R")
source("03_disruption_score.R")
source("04_vps72_residual.R")
source("05_mutation_association.R")
source("06_discrimination_auc.R")
source("07_tp53_variants.R")
source("08_comutation_fisher.R")
source("09_pathway_enrichment.R")
```

Scripts are organized by analysis dependency. The co-mutation workflow (`08`) can also be run independently.

## Statistical Methods

Somatic mutation associations are assessed using logistic regression adjusted for age, sex, and cancer type. Recurrent TP53 variants are evaluated using multivariable linear regression.

Chromatin-complex associations are examined through expression and copy-number summaries, composite disruption scores, and VPS72 residual models adjusted for SRCAP/TIP60 signatures and cancer type.

Cohort-specific mutation enrichment is evaluated using Fisher's exact tests. Pathway over-representation is assessed using GO, KEGG, and Hallmark gene sets. Multiple-testing correction uses the Benjamini–Hochberg procedure.

## Outputs

Scripts generate tab-delimited tables containing integrated patient-level data, mutation associations, chromatin-complex scores, residual model estimates, discrimination metrics, TP53 variant associations, co-mutation statistics, and pathway enrichment results.

**Requirements:** R ≥ 4.2; `data.table`, `broom`, `pROC`, `PRROC`, `clusterProfiler`, `org.Hs.eg.db`, and `msigdbr`.