# Epigenetic regulator disruption ranking

Analysis code accompanying *Multi-omic analyses reveal VPS72 associations with chromatin remodeling, proliferation, and stress-response programs across human cancers*.

This workflow ranks a curated panel of epigenetic regulators using TCGA pan-cancer expression and copy-number measurements. It reports VPS72 and other SRCAP/TIP60/INO80-related genes within the ranking.

## Requirements

- R >= 4.2; packages: `data.table`, `matrixStats`, `AnnotationDbi`, `org.Hs.eg.db`.
- TCGA pan-cancer expression matrix: `EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv`.
- TCGA GISTIC2 matrix: `GDC-PANCAN.gistic.tsv`.
- TCGA barcode-to-cancer-type annotations: `merged_sample_quality_annotations.tsv`.

The TCGA source datasets are available through UCSC Xena and the NCI Genomic Data Commons. Input files are not redistributed in this folder. See `00_config.R` for the expected local input paths and gene panel.

## Run

1. Place the three input files at the locations configured in `00_config.R`.
2. Set `PROJ_DIR` in `00_config.R`.
3. From this analysis folder, run:

```r
source("01_epi_ranking.R")
```

## Analysis

For each gene and cancer type, the workflow summarizes expression (median expression z-score), copy number (mean GISTIC2 score), and high-level amplification frequency (fraction with GISTIC2 >= 2). Cancer-type summaries are averaged with equal weighting. Each component is standardized across the analyzed regulator panel, and the standardized components are summed to obtain the composite disruption score.

The analyzed gene set, eligibility thresholds, and handling of missing values are specified in `00_config.R` and `01_epi_ranking.R`.

## Output files

Results are written to `output/`:

| File | Contents |
|---|---|
| `epi_regulator_disruption_ranking.tsv` | Composite scores, component scores, and regulator ranks |
| `expression_by_cancer_type.tsv` | Expression summaries by gene and cancer type |
| `cnv_by_cancer_type.tsv` | GISTIC2 and amplification summaries by gene and cancer type |
| `genes_excluded.tsv` | Genes excluded by analysis eligibility criteria |

Gene-identifier mappings may be cached in `output/cache/` for repeat runs.
