# Driver mutation frequency in VPS72-high vs VPS72-low tumours, per cohort.
# Marginal test: no covariates, so within one lineage this asks whether the
# groups differ, not why. VPS72-high is the top quartile within each cohort.
source("00_config.R")

if (!dir.exists(inputs$cohorts)) stop("Cannot see the cohort folder: ",inputs$cohorts,
                                      "\nDownload the per-study cBioPortal archives listed in the README.")

run_cohort <- function(label,study) {
  dir <- file.path(inputs$cohorts,study)
  f_expr <- file.path(dir,"data_mrna_seq_v2_rsem.txt")
  f_mut <- file.path(dir,"data_mutations.txt")
  if (!all(file.exists(f_expr,f_mut))) { cat("  skipped, files not found:",label,"\n"); return(NULL) }

  expr <- fread(f_expr)
  vps <- expr[toupper(Hugo_Symbol)=="VPS72"]
  if (!nrow(vps)) { cat("  VPS72 not in expression file:",label,"\n"); return(NULL) }
  cols <- setdiff(names(vps),c("Hugo_Symbol","Entrez_Gene_Id"))
  grp <- data.table(sample=cols,VPS72=as.numeric(vps[1L,..cols]))[is.finite(VPS72)]
  grp[,group:=fifelse(VPS72>=quantile(VPS72,high_quantile,na.rm=TRUE),"High","Low")]

  mut <- fread(f_mut,skip="Hugo_Symbol")
  mut <- unique(mut[Tumor_Sample_Barcode %chin% grp$sample & toupper(Hugo_Symbol) %chin% common_drivers,
                    .(sample=Tumor_Sample_Barcode,gene=toupper(Hugo_Symbol))])
  if (!nrow(mut)) { cat("  no driver mutations:",label,"\n"); return(NULL) }
  mut[,mutated:=1L]

  # Samples absent from the MAF are wild-type, not missing. Left as NA they
  # drop out of the contingency table and inflate mutation frequency.
  wide <- dcast(mut,sample~gene,value.var="mutated",fill=0L)
  dat <- merge(grp,wide,by="sample",all.x=TRUE)
  genes <- setdiff(names(wide),"sample")
  for (g in genes) dat[,(g):=fifelse(is.na(get(g)),0L,get(g))]

  out <- rbindlist(lapply(genes,function(g) {
    tab <- table(factor(dat[[g]],levels=c(0L,1L)),factor(dat$group,levels=c("Low","High")))
    if (sum(tab[2L,])<1L) return(NULL)
    ft <- fisher.test(tab)
    data.table(cohort=label,gene=g,
               mutated_high=tab[2L,"High"],total_high=sum(tab[,"High"]),
               mutated_low=tab[2L,"Low"],total_low=sum(tab[,"Low"]),
               odds_ratio=as.numeric(ft$estimate),
               ci_low=ft$conf.int[1L],ci_high=ft$conf.int[2L],p=ft$p.value)
  }))
  if (!nrow(out)) return(NULL)
  cat("  ",label,": ",nrow(out)," genes, n = ",nrow(dat),"\n",sep="")
  out[,FDR:=p.adjust(p,method="BH")]   # within cohort; each cohort is its own question
  out[]
}

res <- rbindlist(Map(run_cohort,names(cohorts),cohorts))
if (!length(res) || !nrow(res)) stop("No cohort produced results; check inputs$cohorts in 00_config.R")
setorder(res,cohort,p)

write_tsv(res,"comutation_fisher.tsv")
cat("Cohorts:",uniqueN(res$cohort),"; gene-cohort pairs at FDR <",fdr_cut,":",
    sum(res$FDR<fdr_cut,na.rm=TRUE),"\n")
