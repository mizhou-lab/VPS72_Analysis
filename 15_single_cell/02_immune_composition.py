#!/usr/bin/env python3
"""VPS72 high/low immune composition in the merged pan-cancer single-cell atlas.

Analysis only: sample-level cell-type fractions, Mann-Whitney tests and BH FDR.
A sample is High/Low when its median tumor-epithelial VPS72 expression crosses
GLOBAL CELL-LEVEL q80/q20 (the original exploratory definition).
Percentages have a denominator of the eight listed immune-cell populations.
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.stats import mannwhitneyu
from config import H5AD_FILE, IMMUNE_OUT, IMMUNE_MIN_CELLS_PER_SAMPLE

INPUT = str(H5AD_FILE)
OUT = str(IMMUNE_OUT)
os.makedirs(OUT, exist_ok=True)
VPS72_FIELD = "VPS72_expr"
MIN_IMMUNE_CELLS = IMMUNE_MIN_CELLS_PER_SAMPLE
IMMUNE_ORDER = ["T cell", "Macrophage", "NK cell", "B cell",
                "Plasma cell", "Dendritic cell", "Mast", "pDC"]
GROUP_ORDER = ["VPS72-High", "VPS72-Low"]
if not H5AD_FILE.is_file():
    raise FileNotFoundError(f"Input H5AD not found: {H5AD_FILE}")

# 3. Read AnnData
############################################################
print("\n==============================================")
print("Reading:",INPUT)
adata=sc.read_h5ad(INPUT)

print("Cells:",adata.n_obs)
print("Genes:",adata.n_vars)
print("Samples:",adata.obs["Sample"].nunique())
print("Patients:",adata.obs["Patient"].nunique())
print("==============================================")

############################################################

# 4. Check required metadata
############################################################
required=[
"Dataset","Organ_origin","Sample","Patient",
"Tissue","Cancer type","cnv_status",
"Celltype","source_dataset"
]

for c in required:
    if c not in adata.obs.columns:
        raise ValueError(f"Missing metadata column: {c}")

if "VPS72" not in adata.var_names:
    raise ValueError("VPS72 not found in adata.var_names")
# Grouping uses the Sample field. Fail rather than silently combine different
# datasets/patients that reuse the same sample label.
links = adata.obs[["Sample", "Dataset", "Patient"]].drop_duplicates()
if (links.groupby("Sample", observed=True).size() > 1).any():
    raise ValueError("Sample IDs are not unique across Dataset/Patient. "
                     "Make globally unique sample IDs before this analysis.")


############################################################

# 5. Extract VPS72 expression
############################################################
x=adata[:,["VPS72"]].X
if sparse.issparse(x):
    x=x.toarray()

adata.obs[VPS72_FIELD]=np.asarray(x).reshape(-1)

print("\nVPS72 expression summary:")
print(adata.obs[VPS72_FIELD].describe())

############################################################

# 6. Define tumor epithelial cells
############################################################
tumor_epi_mask=(
(adata.obs["Celltype"].astype(str)=="Epithelial") &
(adata.obs["cnv_status"].astype(str).str.lower()=="tumor")
)

n_epi=int(tumor_epi_mask.sum())

if n_epi==0:
    raise ValueError("No tumor epithelial cells found")

print("\nTumor epithelial cells:",n_epi)

############################################################

# 7. Calculate q20/q80 exactly as previous VPS72 analysis
############################################################
vps72=adata.obs.loc[tumor_epi_mask,VPS72_FIELD]

q_low=vps72.quantile(0.20)
q_high=vps72.quantile(0.80)

if not (np.isfinite(q_low) and np.isfinite(q_high) and q_low < q_high):
    raise ValueError(f"Overlapping VPS72 quantiles: q20={q_low}, q80={q_high}. "
                     "Cannot define disjoint high/low groups.")

high_mask_epi=(
tumor_epi_mask &
(adata.obs[VPS72_FIELD]>=q_high)
)

low_mask_epi=(
tumor_epi_mask &
(adata.obs[VPS72_FIELD]<=q_low)
)

print("\nVPS72 thresholds:")
print("q20 =",q_low)
print("q80 =",q_high)

print("\nExtreme tumor epithelial cells:")
print("High =",int(high_mask_epi.sum()))
print("Low  =",int(low_mask_epi.sum()))
print("Middle excluded =",int(n_epi-high_mask_epi.sum()-low_mask_epi.sum()))

############################################################

# 8. Save extreme tumor epithelial cells
############################################################
epi_extreme=adata.obs.loc[
high_mask_epi | low_mask_epi,
[
"Dataset","Organ_origin","Sample","Patient",
"Tissue","Cancer type","cnv_status",
"Celltype","source_dataset",VPS72_FIELD
]
].copy()

epi_extreme["VPS72_group"]=np.where(
epi_extreme[VPS72_FIELD]>=q_high,
"VPS72-High",
"VPS72-Low"
)

epi_extreme.to_csv(
os.path.join(OUT,"01_tumor_epithelial_VPS72_HL_cells.csv")
)

############################################################

# 9. Tumor epithelial VPS72 per sample
############################################################
epi_sample=adata.obs.loc[
tumor_epi_mask,
[
"Sample","Patient","Cancer type",
"Organ_origin","Dataset",VPS72_FIELD
]
].copy()

sample_vps72=(
epi_sample
.groupby("Sample",observed=True)
.agg(
VPS72_median=(VPS72_FIELD,"median"),
VPS72_mean=(VPS72_FIELD,"mean"),
Tumor_Epithelial_Cells=(VPS72_FIELD,"size"),
Patient=("Patient","first"),
Cancer_type=("Cancer type","first"),
Organ_origin=("Organ_origin","first"),
Dataset=("Dataset","first")
)
.reset_index()
)

############################################################

# 10. Assign VPS72 High/Low samples
############################################################
sample_vps72["VPS72_group"]="Middle"

sample_vps72.loc[
sample_vps72["VPS72_median"]>=q_high,
"VPS72_group"
]="VPS72-High"

sample_vps72.loc[
sample_vps72["VPS72_median"]<=q_low,
"VPS72_group"
]="VPS72-Low"

sample_vps72.to_csv(
os.path.join(OUT,"02_all_samples_VPS72_status.csv"),
index=False
)

selected_samples=sample_vps72[
sample_vps72["VPS72_group"].isin(GROUP_ORDER)
].copy()

selected_samples.to_csv(
os.path.join(OUT,"03_selected_VPS72_HL_samples.csv"),
index=False
)

print("\nSample grouping:")
print(sample_vps72["VPS72_group"].value_counts().to_string())

############################################################

# 11. Sample -> VPS72 group mapping
############################################################
sample_group=dict(zip(
selected_samples["Sample"],
selected_samples["VPS72_group"]
))

############################################################

# 12. Select immune cells
############################################################
immune_mask=(
adata.obs["Celltype"].astype(str).isin(IMMUNE_ORDER) &
adata.obs["Sample"].isin(selected_samples["Sample"])
)

immune=adata.obs.loc[
immune_mask,
[
"Dataset","Organ_origin","Sample","Patient",
"Tissue","Cancer type","cnv_status",
"Celltype","source_dataset"
]
].copy()

immune["VPS72_group"]=immune["Sample"].map(sample_group)

immune=immune[
immune["VPS72_group"].isin(GROUP_ORDER)
].copy()

immune.to_csv(
os.path.join(OUT,"04_VPS72_HL_immune_cells.csv")
)

print("\nSelected immune cells:",len(immune))

############################################################

# 13. Immune-cell counts per sample
############################################################
sample_counts=(
immune
.groupby(["Sample","VPS72_group","Celltype"],observed=True)
.size()
.reset_index(name="Cell_Count")
)

############################################################

# 14. Create complete Sample x Celltype table
############################################################
comb=pd.MultiIndex.from_product(
[selected_samples["Sample"].unique(),IMMUNE_ORDER],
names=["Sample","Celltype"]
).to_frame(index=False)

comb=comb.merge(
selected_samples[["Sample","VPS72_group"]],
on="Sample",
how="left"
)

sample_counts_complete=comb.merge(
sample_counts,
on=["Sample","VPS72_group","Celltype"],
how="left"
)

sample_counts_complete["Cell_Count"]=(
sample_counts_complete["Cell_Count"]
.fillna(0)
.astype(int)
)

############################################################

# 15. Calculate immune-cell percentages
############################################################
sample_counts_complete["Total_Immune_Cells"]=(
sample_counts_complete
.groupby("Sample")["Cell_Count"]
.transform("sum")
)

sample_counts_complete["Percent"]=np.where(
sample_counts_complete["Total_Immune_Cells"]>0,
sample_counts_complete["Cell_Count"]/
sample_counts_complete["Total_Immune_Cells"]*100,
np.nan
)

############################################################

# 16. Add sample metadata
############################################################
sample_meta=selected_samples[
[
"Sample","Patient","Cancer_type","Organ_origin",
"Dataset","Tumor_Epithelial_Cells",
"VPS72_median","VPS72_mean"
]
]

sample_counts_complete=sample_counts_complete.merge(
sample_meta,
on="Sample",
how="left"
)

sample_counts_complete.to_csv(
os.path.join(OUT,"05_immune_population_per_sample_ALL.csv"),
index=False
)

############################################################

# 17. Sample QC
############################################################
sample_qc=(
sample_counts_complete[
["Sample","Total_Immune_Cells"]
]
.drop_duplicates()
)

valid_samples=sample_qc.loc[
sample_qc["Total_Immune_Cells"]>=MIN_IMMUNE_CELLS,
"Sample"
]

sample_qc_df=sample_counts_complete[
sample_counts_complete["Sample"].isin(valid_samples)
].copy()

sample_qc_df.to_csv(
os.path.join(OUT,"06_immune_population_per_sample_QC.csv"),
index=False
)

print("\nSample QC:")
print("Before:",sample_counts_complete["Sample"].nunique())
print("After :",sample_qc_df["Sample"].nunique())

print("\nSamples by VPS72 group:")
print(
sample_qc_df[["Sample","VPS72_group"]]
.drop_duplicates()["VPS72_group"]
.value_counts()
.to_string()
)

############################################################

# 18. Mann-Whitney statistics
############################################################
stats=[]

for ct in IMMUNE_ORDER:
    h=sample_qc_df.loc[
    (sample_qc_df["Celltype"]==ct) &
    (sample_qc_df["VPS72_group"]=="VPS72-High"),
    "Percent"
    ].dropna()

    l=sample_qc_df.loc[
    (sample_qc_df["Celltype"]==ct) &
    (sample_qc_df["VPS72_group"]=="VPS72-Low"),
    "Percent"
    ].dropna()

    if len(h)>0 and len(l)>0:
        u,p=mannwhitneyu(h,l,alternative="two-sided")
    else:
        u,p=np.nan,np.nan

    stats.append({
    "Celltype":ct,
    "High_N_samples":len(h),
    "Low_N_samples":len(l),
    "High_mean_percent":h.mean(),
    "Low_mean_percent":l.mean(),
    "High_median_percent":h.median(),
    "Low_median_percent":l.median(),
    "Difference_High_minus_Low":h.mean()-l.mean(),
    "MannWhitney_U":u,
    "P_value":p
    })

stats=pd.DataFrame(stats)

############################################################

# 19. BH-FDR correction
############################################################
stats["FDR"]=np.nan
valid=stats["P_value"].notna()

if valid.sum()>0:
    p=stats.loc[valid,"P_value"].values
    order=np.argsort(p)
    ranked=p[order]
    m=len(ranked)

    adj=ranked*m/np.arange(1,m+1)
    adj=np.minimum.accumulate(adj[::-1])[::-1]
    adj=np.minimum(adj,1)

    corrected=np.empty(m)
    corrected[order]=adj

    stats.loc[valid,"FDR"]=corrected

def significance(x):
    if pd.isna(x): return "NA"
    if x<0.0001: return "****"
    if x<0.001: return "***"
    if x<0.01: return "**"
    if x<0.05: return "*"
    return "ns"

stats["Significance"]=stats["FDR"].apply(significance)

stats.to_csv(
os.path.join(OUT,"07_immune_population_statistics_FDR.csv"),
index=False
)

############################################################

# 22. Pooled immune composition after QC
############################################################
immune_qc=immune[
immune["Sample"].isin(valid_samples)
].copy()

pooled=(
immune_qc
.groupby(["VPS72_group","Celltype"],observed=True)
.size()
.reset_index(name="Cell_Count")
)

pooled["Total"]=(
pooled.groupby("VPS72_group")["Cell_Count"]
.transform("sum")
)

pooled["Percent"]=(
pooled["Cell_Count"]/pooled["Total"]*100
)

pooled.to_csv(
os.path.join(OUT,"08_pooled_immune_composition_HL.csv"),
index=False
)

############################################################

# 29. Mean/median immune summary
############################################################
immune_summary=(
sample_qc_df
.groupby(["VPS72_group","Celltype"],observed=True)
.agg(
Mean_Percent=("Percent","mean"),
Median_Percent=("Percent","median"),
SD_Percent=("Percent","std"),
N_Samples=("Sample","nunique")
)
.reset_index()
)

immune_summary.to_csv(
os.path.join(OUT,"09_immune_population_mean_median_summary.csv"),
index=False
)

############################################################

# 30. Sample summary
############################################################
sample_summary=(
sample_qc_df[
[
"Sample","VPS72_group",
"Patient","Cancer_type",
"VPS72_median"
]
]
.drop_duplicates()
.groupby("VPS72_group",observed=True)
.agg(
Samples=("Sample","nunique"),
Patients=("Patient","nunique"),
Median_VPS72=("VPS72_median","median")
)
.reset_index()
)

sample_summary.to_csv(
os.path.join(OUT,"10_VPS72_HL_sample_summary.csv"),
index=False
)

############################################################

# 31. Cancer-type summary
############################################################
cancer_summary=(
sample_qc_df[
["Sample","Cancer_type","VPS72_group"]
]
.drop_duplicates()
.groupby(["Cancer_type","VPS72_group"],observed=True)
.size()
.reset_index(name="Samples")
)

cancer_summary.to_csv(
os.path.join(OUT,"11_VPS72_HL_cancer_type_summary.csv"),
index=False
)

############################################################

# 32. Sample-level immune-cell count summary
############################################################
sample_total_summary=(
sample_qc_df[
[
"Sample","VPS72_group",
"Patient","Cancer_type",
"Total_Immune_Cells"
]
]
.drop_duplicates()
)

sample_total_summary.to_csv(
os.path.join(OUT,"12_sample_total_immune_cells.csv"),
index=False
)

############################################################

# 33. Wide count table
############################################################
count_wide=(
sample_qc_df
.pivot_table(
index=["Sample","VPS72_group"],
columns="Celltype",
values="Cell_Count",
aggfunc="sum",
fill_value=0
)
.reset_index()
)

count_wide.to_csv(
os.path.join(OUT,"13_immune_counts_wide.csv"),
index=False
)

############################################################

# 34. Wide percentage table
############################################################
percent_wide=(
sample_qc_df
.pivot_table(
index=["Sample","VPS72_group"],
columns="Celltype",
values="Percent",
aggfunc="sum",
fill_value=0
)
.reset_index()
)

percent_wide.to_csv(
os.path.join(OUT,"14_immune_percent_wide.csv"),
index=False
)

############################################################

# 35. Save analysis parameters
############################################################
with open(
os.path.join(OUT,"00_analysis_parameters.txt"),
"w"
) as f:

    f.write(f"Input={INPUT}\n")
    f.write(f"Tumor_epithelial_cells={n_epi}\n")
    f.write(f"q20={q_low}\n")
    f.write(f"q80={q_high}\n")

    f.write(
    "Tumor_epithelial_definition="
    "Celltype Epithelial + cnv_status tumor\n"
    )

    f.write(
    "High_cell_definition="
    "tumor epithelial VPS72 >= q80\n"
    )

    f.write(
    "Low_cell_definition="
    "tumor epithelial VPS72 <= q20\n"
    )

    f.write(
    "High_sample_definition="
    "median tumor epithelial VPS72 >= global q80\n"
    )

    f.write(
    "Low_sample_definition="
    "median tumor epithelial VPS72 <= global q20\n"
    )

    f.write("Middle_samples=excluded\n")
    f.write(f"Minimum_immune_cells_per_sample={MIN_IMMUNE_CELLS}\n")
    f.write("Statistics=Mann-Whitney U + BH-FDR\n")


############################################################

print("\nSingle-cell immune-composition analysis complete.")
print("Output:", OUT)
print("VPS72 q20/q80:", q_low, q_high)
print("Sample groups after immune-cell QC:")
print(sample_summary.to_string(index=False))
print("Mann-Whitney + BH-FDR:")
print(stats[["Celltype", "P_value", "FDR"]].to_string(index=False))
