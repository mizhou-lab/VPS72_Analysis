# VPS72 Pan-Cancer Single-Cell RNA-seq Analysis

Python scripts supporting the analysis of VPS72 expression, differential gene expression, Hallmark pathway enrichment, and immune-cell composition across pan-cancer single-cell RNA-seq datasets.

These analyses accompany *Multi-omic analyses reveal VPS72 associations with chromatin remodeling, proliferation, and stress-response programs across human cancers*.

## Analysis Scripts

| Script | Analysis |
|---|---|
| `config.py` | Input configuration, analysis parameters, and statistical thresholds |
| `01_epithelial_deg_hallmark.py` | Single-cell cohort QC, malignant versus non-malignant epithelial differential expression, VPS72-high/low comparisons, and Hallmark enrichment |
| `02_immune_composition.py` | Immune-cell composition associated with VPS72 expression, sample-level comparisons, and FDR correction |

## Data Requirements

Analyses use the merged pan-cancer single-cell RNA-seq AnnData object:

`full_32datasets_with_umap.h5ad`

The dataset must contain normalized gene expression, VPS72 expression, cell-type annotations, inferred copy-number status, sample and patient identifiers, and cancer-type metadata.

Hallmark enrichment uses the MSigDB Hallmark 2020 collection.

## Execution

Configure the input and output directories in `config.py`, then run:

```bash
python 01_epithelial_deg_hallmark.py
python 02_immune_composition.py
```

Both analyses can be executed independently.

## Statistical Methods

Malignant epithelial cells are identified using epithelial annotation and inferred tumor copy-number status. VPS72-high and VPS72-low malignant cells are defined using the upper and lower expression quantiles (80th and 20th percentiles).

Differential expression is evaluated using Scanpy's Wilcoxon rank-sum method, with FDR <escape><</escape> 0.05 and absolute log2 fold change ≥ 0.25. Hallmark pathway over-representation is evaluated using hypergeometric tests.

Immune-cell composition is compared between VPS72-high and VPS72-low tumor samples using two-sided Mann–Whitney tests. Samples with fewer than 50 classified immune cells are excluded. Multiple-testing correction uses the Benjamini–Hochberg procedure.

## Outputs

Scripts generate CSV files containing cohort characteristics, differential-expression results, VPS72-high/low group summaries, Hallmark pathway enrichment statistics, and sample-level immune-cell composition analyses.

**Requirements:** Python ≥ 3.9; `scanpy`, `anndata`, `numpy`, `pandas`, `scipy`, and `gseapy`.