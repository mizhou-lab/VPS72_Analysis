# VPS72 dependency across cell lines, CRISPR and RNAi side by side, with the
# fraction of lines dependent at both thresholds and a common-essential versus
# selective classification.
source("00_config.R")

model <- read_model()
crispr <- read_gene_export(inputs$crispr,"crispr_chronos")
rnai   <- read_gene_export(inputs$rnai,"rnai_demeter2")

dep <- merge(model,crispr,by="ModelID",all.y=TRUE)
dep <- merge(dep,rnai,by="ModelID",all=TRUE)
dep <- dep[!is.na(crispr_chronos) | !is.na(rnai_demeter2)]
if (!nrow(dep)) stop("No cell line carried a VPS72 dependency score")
write_tsv(dep,"vps72_dependency.tsv")

# The headline numbers the dependency figure asserts. Both thresholds are
# reported: whether VPS72 reads as a near-universal dependency or a selective
# one depends entirely on which line is drawn, so neither is left implicit.
frac <- function(v,label) {
  v <- v[is.finite(v)]
  if (!length(v)) return(NULL)
  data.table(assay=label,n_lines=length(v),
             median=median(v),mean=mean(v),
             q05=quantile(v,0.05),q95=quantile(v,0.95),
             n_below_loose=sum(v<=dep_loose),pct_below_loose=100*mean(v<=dep_loose),
             n_below_strong=sum(v<=dep_strong),pct_below_strong=100*mean(v<=dep_strong))
}
summ <- rbindlist(list(frac(dep$crispr_chronos,"CRISPR_Chronos"),
                       frac(dep$rnai_demeter2,"RNAi_DEMETER2")))
summ[,`:=`(threshold_loose=dep_loose,threshold_strong=dep_strong)]
write_tsv(summ,"dependency_summary.tsv")
print(summ)

# Common-essential or selective. A gene depleted in nearly every line is a
# common essential; one depleted in a subset is selective, and only the second
# is a candidate target. The spread of scores separates the two, so it is
# reported rather than asserted.
cr <- dep$crispr_chronos[is.finite(dep$crispr_chronos)]
if (length(cr)>=20L) {
  classification <- data.table(
    assay="CRISPR_Chronos",
    pct_lines_below_strong=100*mean(cr<=dep_strong),
    pct_lines_below_loose=100*mean(cr<=dep_loose),
    iqr=IQR(cr),
    sd=sd(cr),
    skewness=mean((cr-mean(cr))^3)/sd(cr)^3,
    reading=fcase(
      mean(cr<=dep_strong)>0.90,"common_essential: depleted at the strong threshold in >90% of lines",
      mean(cr<=dep_loose)>0.90,"broadly_depleted: below the loose threshold in >90% of lines, but few reach the strong one",
      mean(cr<=dep_loose)>0.50,"frequently_depleted: over half of lines below the loose threshold",
      default="selective: fewer than half of lines below the loose threshold"))
  write_tsv(classification,"dependency_classification.tsv")
  cat("\n",classification$reading,"\n\n",sep="")
}

# Per-lineage medians, so a lineage-restricted dependency would be visible.
by_lineage <- dep[!is.na(lineage) & nzchar(lineage),
                  .(n=.N,
                    median_crispr=median(crispr_chronos,na.rm=TRUE),
                    median_rnai=median(rnai_demeter2,na.rm=TRUE),
                    pct_below_loose=100*mean(crispr_chronos<=dep_loose,na.rm=TRUE),
                    pct_below_strong=100*mean(crispr_chronos<=dep_strong,na.rm=TRUE)),
                  by=lineage][n>=min_lineage_n]
setorder(by_lineage,median_crispr)
write_tsv(by_lineage,"dependency_by_lineage.tsv")

# Agreement between the two screens. They use different reagents and scales, so
# a positive correlation is evidence the dependency is real rather than an
# artefact of either platform.
agree <- spearman(dep$crispr_chronos,dep$rnai_demeter2)
if (!is.null(agree)) {
  agree[,comparison:="CRISPR_vs_RNAi"]
  write_tsv(agree,"crispr_rnai_agreement.tsv")
  cat("CRISPR vs RNAi: rho =",round(agree$rho,3),"over",agree$n,"shared lines, p =",
      signif(agree$p,3),"\n")
}

# VPS72 against SRCAP in the same screen, from the full gene-effect matrix.
if (file.exists(inputs$crispr_matrix)) {
  ge <- fread(inputs$crispr_matrix)
  setnames(ge,1,"ModelID")
  syms <- strip_entrez(names(ge))
  pick <- function(g) { i <- which(syms==g); if (!length(i)) NULL else names(ge)[i[1L]] }
  cv <- pick("VPS72"); cs <- pick("SRCAP")
  if (is.null(cv) || is.null(cs)) {
    cat("VPS72 or SRCAP absent from the gene-effect matrix; pair skipped\n")
  } else {
    pair <- data.table(ModelID=ge$ModelID,
                       VPS72=suppressWarnings(as.numeric(ge[[cv]])),
                       SRCAP=suppressWarnings(as.numeric(ge[[cs]])))
    pair <- merge(model[,.(ModelID,lineage)],pair,by="ModelID")
    pair <- pair[is.finite(VPS72) & is.finite(SRCAP)]
    write_tsv(pair,"vps72_srcap_dependency_pair.tsv")
    r <- spearman(pair$VPS72,pair$SRCAP)
    if (!is.null(r)) {
      r[,comparison:="VPS72_vs_SRCAP_chronos"]
      write_tsv(r,"vps72_srcap_correlation.tsv")
      cat("VPS72 vs SRCAP dependency: rho =",round(r$rho,3),"over",r$n,"lines\n")
    }
  }
}
