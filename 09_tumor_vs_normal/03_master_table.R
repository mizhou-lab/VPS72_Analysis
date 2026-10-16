# One harmonized sample table: GEO normal references plus TCGA tumour and
# adjacent normal, with VPS72 standardized within each cancer.
source("00_config.R")

tcga <- fread(require_input(inputs$tcga))
meta <- fread(require_input(inputs$geo_meta))

# TCGA is standardized across tumour AND normal together within each cancer.
# Standardizing the two states separately would remove the tumour-normal
# difference the analysis is about.
tcga[,VPS72_z:=as.numeric(scale(VPS72)),by=cancer]
tcga[,VPS72_pct:=frank(VPS72,ties.method="average")/.N,by=cancer]
tcga <- merge(tcga,tissue_map,by="cancer",all.x=TRUE)
if (anyNA(tcga$tissue)) stop("Unmapped cancer code: ",
                             paste(unique(tcga[is.na(tissue)]$cancer),collapse=", "))

tcga_part <- tcga[,.(dataset=paste0("TCGA-",cancer),source="TCGA",
                     source_type=fifelse(state=="Tumor","tumor","adjacent_normal"),
                     sample,patient,cancer,tissue,state,
                     VPS72,VPS72_z,VPS72_pct,paired,primary_record)]

# The GEO arm carries no cancer label and no paired structure: its VPS72 values
# were standardized within series in 01 and are not comparable to the TCGA
# z-scores, so it is kept as a separate reference rather than pooled into the
# per-cancer contrasts.
geo_vps <- fread(require_input(inputs$geo_expr))[gene=="VPS72"]
if (!nrow(geo_vps)) stop("VPS72 not found in the harmonized GEO expression table")
geo_part <- merge(meta,geo_vps[,.(dataset,sample,VPS72=expr,VPS72_z=expr_z)],
                  by=c("dataset","sample"))
geo_part[,`:=`(source_type="normal_reference",patient=NA_character_,
               cancer=NA_character_,tissue=NA_character_,
               VPS72_pct=NA_real_,paired=FALSE,primary_record=TRUE)]

cols <- c("dataset","source","source_type","sample","patient","cancer","tissue",
          "state","VPS72","VPS72_z","VPS72_pct","paired","primary_record")
master <- rbindlist(list(tcga_part[,..cols],geo_part[,..cols]),fill=TRUE)
master[,label:=fifelse(source=="TCGA",paste0(cancer," ",state),dataset)]

write_tsv(master,"master_sample_table.tsv")
cat("Total samples:",nrow(master),
    "; TCGA:",sum(master$source=="TCGA"),
    "; GEO normal reference:",sum(master$source=="GEO"),"\n")
print(master[,.(samples=.N),by=.(source,source_type,state)])
