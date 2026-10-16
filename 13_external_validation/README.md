# VPS72 External Validation Analysis

R scripts supporting the reanalysis of publicly available VPS72/YL1 perturbation datasets to evaluate transcriptional changes, pathway enrichment, and histone H4 acetylation in melanoma cells.

These analyses accompany *Multi-omic analyses reveal VPS72 associations with chromatin remodeling, proliferation, and stress-response programs across human cancers*.

## Analysis Scripts

| Script | Analysis |
|---|---|
| `00_config.R` | Input configuration, gene signatures, analysis parameters, and statistical thresholds |
| `01_download_geo.R` | Retrieval and preparation of public GEO datasets |
| `02_deg_tables.R` | Standardization and comparison of published differential-expression results |
| `03_overlap_partitions.R` | Gene overlap, perturbation-specific transcriptional programs, and Hallmark enrichment |
| `04_fgsea_hallmark.R` | Hallmark gene-set enrichment and normalized enrichment score concordance |
| `05_chip_h4ac.R` | Promoter-associated H4 acetylation changes following VPS72/YL1 depletion |
| `06_integrate.R` | Integration of transcriptional responses with promoter H4 acetylation changes |

## Data Requirements

Public datasets from Jostes et al. (*Genes & Development*, 2024):

- **GSE242227:** RNA-seq differential-expression results following shYL1, shSRCAP, and shH2AZ1 perturbation in SK-MEL-147 melanoma cells.
- **GSE246121:** H4 acetylation ChIP-seq coverage data from control and YL1-depleted SK-MEL-147 cells.
- **Gene sets:** MSigDB Hallmark and additional signatures specified in `00_config.R`.

The RNA-seq workflow uses the published day-6 DESeq2 results. ChIP-seq analysis uses hg19 genomic coordinates.

## Execution

Configure the input and output directories in `00_config.R`, then execute sequentially:

```r id="pgw9pt"
source("01_download_geo.R")
source("02_deg_tables.R")
source("03_overlap_partitions.R")
source("04_fgsea_hallmark.R")
source("05_chip_h4ac.R")
source("06_integrate.R")
```

The GEO download step requires an internet connection. Subsequent analyses use the downloaded expression tables and ChIP-seq coverage files.

## Statistical Methods

Differentially expressed genes are identified using the published DESeq2 results, applying adjusted P <escape><</escape> 0.05 and absolute log2 fold change > 0.5.

Gene overlap and perturbation-specific transcriptional programs are evaluated through gene-set partitioning and Hallmark over-representation analysis. Hallmark enrichment is additionally assessed using pre-ranked GSEA.

Promoter H4 acetylation is quantified from ChIP-seq coverage within transcription start site windows and compared between control and YL1-depleted cells. Expression and acetylation changes are subsequently integrated at the gene level.

Multiple-testing correction uses the Benjamini–Hochberg procedure.

## Outputs

Scripts generate tab-delimited tables containing standardized differential-expression results, gene overlap partitions, Hallmark enrichment statistics, H4 acetylation changes, and integrated transcriptional and chromatin-associated responses.

**Requirements:** R ≥ 4.2; `data.table`, `GEOquery`, `R.utils`, `fgsea`, `msigdbr`, `clusterProfiler`, `rtracklayer`, `GenomicRanges`, `GenomeInfoDb`, `TxDb.Hsapiens.UCSC.hg19.knownGene`, `org.Hs.eg.db`, and `AnnotationDbi`.