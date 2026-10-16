#!/usr/bin/env python3
"""VPS72 single-cell analysis: cohort QC, epithelial state, DEG and Hallmark ORA.

Analysis only: no UMAP, volcano, dotplot, heatmap, or Excel rendering.
The grouping/DEG/ORA calculations are retained from the submitted master script.
Unit of analysis for DEG: cells (exploratory; cells are NOT independent donors).
Input adata.X: normalized expression, not untransformed raw counts.
"""
import os
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.stats import hypergeom
from config import H5AD_FILE, MASTER_OUT, DEG_FDR, DEG_MIN_ABS_LOG2FC, HALLMARK_GMT

INPUT = str(H5AD_FILE)
OUT = str(MASTER_OUT)
FDR = DEG_FDR
LFC = DEG_MIN_ABS_LOG2FC
HIGH = "VPS72-High"
LOW = "VPS72-Low"
if not H5AD_FILE.is_file():
    raise FileNotFoundError(f"Input H5AD not found: {H5AD_FILE}")
for name in (
    "00_QC", "02_Epithelial_Status", "03_DEG_Malignant_vs_NonMalignant",
    "04_Hallmark_Malignant_vs_NonMalignant", "05_VPS72_HighLow",
    "06_DEG_VPS72_High_vs_Low", "07_Hallmark_VPS72_High_vs_Low",
    "09_Summary",
):
    (MASTER_OUT / name).mkdir(parents=True, exist_ok=True)

# 2. HELPERS
###############################################################################
def save(df,name,folder):
    p=f"{OUT}/{folder}/{name}"
    df.to_csv(p,index=False)
    print("Saved:",p)

def bh(p):
    p=np.asarray(p,float);out=np.full(len(p),np.nan);ok=np.isfinite(p)
    if ok.sum()==0:return out
    x=p[ok];o=np.argsort(x);r=x[o];n=len(r)
    a=np.minimum.accumulate((r*n/np.arange(1,n+1))[::-1])[::-1]
    z=np.empty(n);z[o]=np.minimum(a,1);out[np.where(ok)[0]]=z
    return out

def gene_vector(a,g):
    x=a[:,[g]].X
    if sparse.issparse(x):x=x.toarray()
    return np.asarray(x).ravel()


# 3. READ DATA
###############################################################################
print("\n===== READ DATA =====")
adata=sc.read_h5ad(INPUT)
print("Cells:",adata.n_obs)
print("Genes:",adata.n_vars)

req=["Dataset","Organ_origin","Sample","Patient","Tissue","Cancer type",
     "cnv_status","Celltype","source_dataset"]
missing=[x for x in req if x not in adata.obs.columns]
if missing:raise ValueError("Missing metadata: "+",".join(missing))
if "VPS72" not in adata.var_names:raise ValueError("VPS72 not found")


###############################################################################

# 4. GLOBAL QC
###############################################################################
qc=pd.DataFrame({"Metric":["Cells","Genes","Samples","Patients","Datasets",
"Cancer_types","Organ_origins","Celltypes"],"Value":[adata.n_obs,adata.n_vars,
adata.obs.Sample.nunique(),adata.obs.Patient.nunique(),adata.obs.Dataset.nunique(),
adata.obs["Cancer type"].nunique(),adata.obs.Organ_origin.nunique(),
adata.obs.Celltype.nunique()]})
save(qc,"01_input_summary.csv","00_QC")

save(adata.obs.Celltype.value_counts().rename_axis("Celltype").reset_index(name="Cells"),
"02_celltype_counts.csv","00_QC")
save(pd.crosstab(adata.obs.Celltype,adata.obs.cnv_status).reset_index(),
"03_celltype_x_cnv.csv","00_QC")
save(adata.obs["Cancer type"].value_counts().rename_axis("Cancer_type").reset_index(name="Cells"),
"04_cancer_counts.csv","00_QC")
save(adata.obs.Organ_origin.value_counts().rename_axis("Organ_origin").reset_index(name="Cells"),
"05_organ_counts.csv","00_QC")
save(adata.obs.Dataset.value_counts().rename_axis("Dataset").reset_index(name="Cells"),
"06_dataset_counts.csv","00_QC")

###############################################################################

# 8. EPITHELIAL CELLS
###############################################################################
epi=adata[adata.obs.Celltype.astype(str)=="Epithelial"].copy()
epi.obs["Epithelial_Status"]=epi.obs.cnv_status.astype(str).str.lower().map(
{"tumor":"Malignant Epithelial","normal":"Non-malignant Epithelial"})
epi=epi[epi.obs.Epithelial_Status.notna()].copy()

save(epi.obs.Epithelial_Status.value_counts().rename_axis("Status").reset_index(name="Cells"),
"01_epithelial_status_counts.csv","02_Epithelial_Status")

save(epi.obs.groupby(["Epithelial_Status","Sample"],observed=True).size().reset_index(name="Cells"),
"02_epithelial_status_by_sample.csv","02_Epithelial_Status")
save(epi.obs.groupby(["Epithelial_Status","Cancer type"],observed=True).size().reset_index(name="Cells"),
"03_epithelial_status_by_cancer.csv","02_Epithelial_Status")
save(epi.obs.groupby(["Epithelial_Status","Organ_origin"],observed=True).size().reset_index(name="Cells"),
"04_epithelial_status_by_organ.csv","02_Epithelial_Status")

###############################################################################

# 10. GENERIC DEG FUNCTION
###############################################################################
def run_deg(obj,groupby,group,reference,folder,prefix,pos_name,neg_name):
    print("\n===== DEG:",prefix,"=====")
    sc.tl.rank_genes_groups(obj,groupby=groupby,groups=[group],reference=reference,
    method="wilcoxon",use_raw=False,pts=True,tie_correct=True)

    d=sc.get.rank_genes_groups_df(obj,group=group)
    d=d.rename(columns={"names":"gene","scores":"score","logfoldchanges":"avg_log2FC",
    "pvals":"p_value","pvals_adj":"FDR","pct_nz_group":"pct_group",
    "pct_nz_reference":"pct_reference"})

    for c in ["score","avg_log2FC","p_value","FDR","pct_group","pct_reference"]:
        if c not in d.columns:d[c]=np.nan

    d["Significant"]=(d.FDR<FDR)&(d.avg_log2FC.abs()>=LFC)
    d["Direction"]=np.where(d.Significant&(d.avg_log2FC>0),pos_name,
                    np.where(d.Significant&(d.avg_log2FC<0),neg_name,"NS"))
    d["minusLog10FDR"]=-np.log10(np.clip(d.FDR.astype(float),1e-300,None))
    d=d.sort_values(["FDR","avg_log2FC"],ascending=[True,False])

    sig=d[d.Significant].copy()
    pos=sig[sig.avg_log2FC>0].copy()
    neg=sig[sig.avg_log2FC<0].copy()

    save(d,f"01_{prefix}_ALL_genes.csv",folder)
    save(sig,f"02_{prefix}_SIGNIFICANT.csv",folder)
    save(pos,f"03_{prefix}_{pos_name}.csv",folder)
    save(neg,f"04_{prefix}_{neg_name}.csv",folder)
    save(sig[["gene"]],f"05_{prefix}_significant_gene_names_ONLY.csv",folder)
    save(pos[["gene"]],f"06_{prefix}_{pos_name}_gene_names_ONLY.csv",folder)
    save(neg[["gene"]],f"07_{prefix}_{neg_name}_gene_names_ONLY.csv",folder)

    v=d[d.gene.astype(str).str.upper()=="VPS72"].copy()
    save(v,f"08_{prefix}_VPS72_result.csv",folder)

    idx1=np.where(obj.obs[groupby].astype(str).values==group)[0]
    idx2=np.where(obj.obs[groupby].astype(str).values==reference)[0]
    m1=np.asarray(obj.X[idx1].mean(axis=0)).ravel()
    m2=np.asarray(obj.X[idx2].mean(axis=0)).ravel()

    means=pd.DataFrame({"gene":obj.var_names.astype(str),
    f"{group}_mean":m1,f"{reference}_mean":m2,
    "Difference":m1-m2})
    save(means,f"09_{prefix}_ALL_gene_mean_expression.csv",folder)

    summ=pd.DataFrame({"Metric":["Group1","Group2","Group1_cells","Group2_cells",
    "Genes_tested","Significant_DEGs",pos_name,neg_name,"FDR","Abs_log2FC"],
    "Value":[group,reference,len(idx1),len(idx2),len(d),len(sig),len(pos),len(neg),FDR,LFC]})
    save(summ,f"10_{prefix}_summary.csv",folder)

    return d,sig,pos,neg,v,means

###############################################################################

# 11. COMPARISON A - MALIGNANT vs NON-MALIGNANT
###############################################################################
epi.obs["Epithelial_Status"]=pd.Categorical(epi.obs.Epithelial_Status,
categories=["Malignant Epithelial","Non-malignant Epithelial"])

degA,sigA,posA,negA,vpsA,meansA=run_deg(
epi,"Epithelial_Status","Malignant Epithelial","Non-malignant Epithelial",
"03_DEG_Malignant_vs_NonMalignant","Malignant_vs_NonMalignant",
"Higher_in_Malignant","Higher_in_NonMalignant")

###############################################################################

# 12. VPS72 EXPRESSION + MALIGNANT HIGH/LOW
###############################################################################
mal=epi[epi.obs.Epithelial_Status=="Malignant Epithelial"].copy()
mal.obs["VPS72_expr"]=gene_vector(mal,"VPS72")
q20=float(mal.obs.VPS72_expr.quantile(.20))
q80=float(mal.obs.VPS72_expr.quantile(.80))

if not (np.isfinite(q20) and np.isfinite(q80) and q20 < q80):
    raise ValueError(f"Overlapping VPS72 quantiles: q20={q20}, q80={q80}. Cannot form high/low groups safely.")
mal.obs["VPS72_group"]="Middle"
mal.obs.loc[mal.obs.VPS72_expr>=q80,"VPS72_group"]=HIGH
mal.obs.loc[mal.obs.VPS72_expr<=q20,"VPS72_group"]=LOW

hl=mal[mal.obs.VPS72_group.isin([HIGH,LOW])].copy()
hl.obs["VPS72_group"]=pd.Categorical(hl.obs.VPS72_group,categories=[HIGH,LOW])

save(pd.DataFrame({"Metric":["Malignant_cells","VPS72_positive",
"VPS72_positive_percent","q20","q80","High_cells","Low_cells"],
"Value":[mal.n_obs,(mal.obs.VPS72_expr>0).sum(),
100*(mal.obs.VPS72_expr>0).mean(),q20,q80,
sum(hl.obs.VPS72_group==HIGH),sum(hl.obs.VPS72_group==LOW)]}),
"01_VPS72_HL_summary.csv","05_VPS72_HighLow")

meta=hl.obs.copy();meta.insert(0,"Cell",meta.index)
save(meta,"02_ALL_VPS72_HL_cells.csv","05_VPS72_HighLow")
save(hl.obs.groupby(["Cancer type","VPS72_group"],observed=True).size().reset_index(name="Cells"),
"03_HL_by_cancer.csv","05_VPS72_HighLow")
save(hl.obs.groupby(["Organ_origin","VPS72_group"],observed=True).size().reset_index(name="Cells"),
"04_HL_by_organ.csv","05_VPS72_HighLow")
save(hl.obs.groupby(["Sample","VPS72_group"],observed=True).size().reset_index(name="Cells"),
"05_HL_by_sample.csv","05_VPS72_HighLow")

###############################################################################

# 15. COMPARISON B - VPS72 HIGH vs LOW
###############################################################################
degB,sigB,posB,negB,vpsB,meansB=run_deg(
hl,"VPS72_group",HIGH,LOW,
"06_DEG_VPS72_High_vs_Low","VPS72_High_vs_Low",
"Higher_in_VPS72_High","Higher_in_VPS72_Low")

###############################################################################

# 16. LOAD HALLMARK
###############################################################################
print("\n===== HALLMARK =====")
try:
    import gseapy as gp
    if HALLMARK_GMT is not None:
        hm={}
        with open(HALLMARK_GMT, encoding="utf-8") as fh:
            for line in fh:
                x=line.rstrip("\n").split("\t")
                if len(x) >= 3:
                    hm[x[0]]=x[2:]
        if not hm:
            raise ValueError("No gene sets found in HALLMARK_GMT")
    else:
        hm=gp.get_library(name="MSigDB_Hallmark_2020",organism="Human")
    hallmark={k:set(v) for k,v in hm.items()}
except Exception as e:
    raise RuntimeError("Install gseapy first: pip install gseapy\n"+str(e))

hall_long=pd.DataFrame([(p,g) for p,genes in hallmark.items() for g in genes],
columns=["Hallmark","gene"])
save(hall_long,"01_Hallmark_gene_sets.csv","09_Summary")

###############################################################################

# 17. HALLMARK FUNCTION
###############################################################################
def hallmark_ora(query,background,folder,prefix):
    bg=set(background);q=set(query)&bg;rows=[]
    for pathway,genes in hallmark.items():
        gs=set(genes)&bg;k=len(q&gs);M=len(bg);n=len(gs);N=len(q)
        if M==0 or n==0 or N==0:continue
        p=hypergeom.sf(k-1,M,n,N) if k>0 else 1
        rows.append([pathway,M,n,N,k,p,";".join(sorted(q&gs))])

    r=pd.DataFrame(rows,columns=["Hallmark","Background_N","Pathway_N",
    "Query_N","Overlap_N","P_value","Genes"])

    if len(r):
        r["FDR"]=bh(r.P_value)
        r["GeneRatio"]=r.Overlap_N/r.Query_N
        r["BackgroundRatio"]=r.Pathway_N/r.Background_N
        r["FoldEnrichment"]=r.GeneRatio/r.BackgroundRatio
        r=r.sort_values(["FDR","P_value","FoldEnrichment"],
                        ascending=[True,True,False])

    save(r,f"{prefix}_ALL.csv",folder)
    save(r[r.FDR<.05].copy() if len(r) else r,f"{prefix}_FDR05.csv",folder)
    return r

###############################################################################

# 18. HALLMARK - COMPARISON A
###############################################################################
hallA_all=hallmark_ora(sigA.gene,degA.gene,
"04_Hallmark_Malignant_vs_NonMalignant","01_ALL_significant_DEGs")
hallA_pos=hallmark_ora(posA.gene,degA.gene,
"04_Hallmark_Malignant_vs_NonMalignant","02_Malignant_HIGH_DEGs")
hallA_neg=hallmark_ora(negA.gene,degA.gene,
"04_Hallmark_Malignant_vs_NonMalignant","03_NonMalignant_HIGH_DEGs")

###############################################################################

# 19. HALLMARK - COMPARISON B
###############################################################################
hallB_all=hallmark_ora(sigB.gene,degB.gene,
"07_Hallmark_VPS72_High_vs_Low","01_ALL_significant_DEGs")
hallB_pos=hallmark_ora(posB.gene,degB.gene,
"07_Hallmark_VPS72_High_vs_Low","02_VPS72_HIGH_DEGs")
hallB_neg=hallmark_ora(negB.gene,degB.gene,
"07_Hallmark_VPS72_High_vs_Low","03_VPS72_LOW_DEGs")

###############################################################################

# 27. COMBINED DEG SUMMARY
###############################################################################
summary=pd.DataFrame({
"Comparison":[
"Malignant vs Non-malignant",
"VPS72-High vs VPS72-Low"
],
"Group1_cells":[
sum(epi.obs.Epithelial_Status=="Malignant Epithelial"),
sum(hl.obs.VPS72_group==HIGH)
],
"Group2_cells":[
sum(epi.obs.Epithelial_Status=="Non-malignant Epithelial"),
sum(hl.obs.VPS72_group==LOW)
],
"Genes_tested":[len(degA),len(degB)],
"Significant_DEGs":[len(sigA),len(sigB)],
"Positive_direction_DEGs":[len(posA),len(posB)],
"Negative_direction_DEGs":[len(negA),len(negB)]
})
save(summary,"01_MASTER_DEG_summary.csv","09_Summary")

###############################################################################

# 28. COMBINED VPS72 RESULTS
###############################################################################
vr=[]
if len(vpsA):
    x=vpsA.iloc[0].copy()
    vr.append({"Comparison":"Malignant vs Non-malignant",
    "gene":"VPS72","avg_log2FC":x.avg_log2FC,
    "p_value":x.p_value,"FDR":x.FDR,"Direction":x.Direction})

if len(vpsB):
    x=vpsB.iloc[0].copy()
    vr.append({"Comparison":"VPS72-High vs VPS72-Low",
    "gene":"VPS72","avg_log2FC":x.avg_log2FC,
    "p_value":x.p_value,"FDR":x.FDR,"Direction":x.Direction})

save(pd.DataFrame(vr),"02_VPS72_results_BOTH_comparisons.csv","09_Summary")

###############################################################################

# 29. ANALYSIS PARAMETERS
###############################################################################
params=pd.DataFrame({
"Parameter":[
"Input","DEG_method","Expression_matrix","FDR_cutoff",
"Abs_log2FC_cutoff","Malignant_definition",
"NonMalignant_definition","VPS72_High_definition",
"VPS72_Low_definition","VPS72_q20","VPS72_q80",
"Hallmark_library"
],
"Value":[
INPUT,"Scanpy Wilcoxon","normalized adata.X",FDR,LFC,
'Celltype=Epithelial + cnv_status=tumor',
'Celltype=Epithelial + cnv_status=normal',
"Malignant epithelial VPS72 >= q80",
"Malignant epithelial VPS72 <= q20",
q20,q80,"MSigDB_Hallmark_2020"
]})
save(params,"03_analysis_parameters.csv","09_Summary")

###############################################################################

# 30. OUTPUT MANIFEST
###############################################################################
manifest=[]
for root,_,files in os.walk(OUT):
    for f in files:
        manifest.append({
        "Folder":os.path.relpath(root,OUT),
        "File":f
        })
save(pd.DataFrame(manifest).sort_values(["Folder","File"]),
"04_OUTPUT_MANIFEST.csv","09_Summary")

###############################################################################


# Additional numerical source data: VPS72 expression by epithelial state.
from scipy.stats import mannwhitneyu
for_status = epi.obs["Epithelial_Status"].astype(str).to_numpy()
expr_vps = gene_vector(epi, "VPS72")
rows = []
for label in ("Malignant Epithelial", "Non-malignant Epithelial"):
    values = expr_vps[for_status == label]
    rows.append({"group": label, "cells": len(values),
                 "mean_expression": float(np.mean(values)),
                 "median_expression": float(np.median(values)),
                 "pct_detected": float(100*np.mean(values > 0))})
save(pd.DataFrame(rows), "03_VPS72_epithelial_expression_summary.csv", "09_Summary")

# 31. FINAL
###############################################################################
print("\n==========================================================")
print("VPS72 MASTER ANALYSIS COMPLETE")
print("==========================================================")
print("Epithelial cells:",epi.n_obs)
print("Malignant:",sum(epi.obs.Epithelial_Status=="Malignant Epithelial"))
print("Non-malignant:",sum(epi.obs.Epithelial_Status=="Non-malignant Epithelial"))
print("")
print("Malignant vs Non-malignant DEGs:",len(sigA))
print("Higher in Malignant:",len(posA))
print("Higher in Non-malignant:",len(negA))
print("")
print("VPS72 q20:",q20)
print("VPS72 q80:",q80)
print("VPS72-High cells:",sum(hl.obs.VPS72_group==HIGH))
print("VPS72-Low cells:",sum(hl.obs.VPS72_group==LOW))
print("VPS72 High vs Low DEGs:",len(sigB))
print("Higher in VPS72-High:",len(posB))
print("Higher in VPS72-Low:",len(negB))
print("")
print("Output:",OUT)
print("==========================================================")