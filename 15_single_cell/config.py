"""Configuration for the two analysis-only VPS72 single-cell pipelines.

Change BASE_DIR to the directory containing full_32datasets_with_umap.h5ad.
Paths here are example placeholders; no personal directories are required.
"""
from pathlib import Path

BASE_DIR = Path("path/to/raw_data")
H5AD_FILE = BASE_DIR / "full_32datasets_with_umap.h5ad"
MASTER_OUT = BASE_DIR / "20260929_VPS72_MASTER_Epithelial"
IMMUNE_OUT = BASE_DIR / "20260908_VPS72_immune_HL"
DEG_FDR = 0.05
DEG_MIN_ABS_LOG2FC = 0.25
IMMUNE_MIN_CELLS_PER_SAMPLE = 50
# Optional: downloaded local version of MSigDB Hallmark 2020 gene sets in GMT format.
# If None, gseapy retrieves the original library (requires internet access).
HALLMARK_GMT = None
