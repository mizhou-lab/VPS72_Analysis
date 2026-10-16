# VPS72 co-dependencies from both screens, merged, with complex membership
# flagged and the two screens' agreement reported.
source("00_config.R")

read_codep <- function(p,value_name) {
  d <- clean_names(fread(require_input(p)))
  if (ncol(d)<2L) stop("Expected at least two columns in ",p)
  setnames(d,1:2,c("gene",value_name))
  d[,gene:=strip_entrez(gene)]
  d[[value_name]] <- suppressWarnings(as.numeric(d[[value_name]]))
  unique(d[nzchar(gene) & !is.na(get(value_name)),.SD,.SDcols=c("gene",value_name)])
}

cr <- read_codep(inputs$crispr_codep,"rho_crispr")
rn <- read_codep(inputs$rnai_codep,"rho_rnai")

# Outer join: a gene in the top 100 of one screen but not the other is
# informative, and dropping it would bias the overlap downward.
dat <- merge(cr,rn,by="gene",all=TRUE)
dat[,max_abs_rho:=pmax(abs(rho_crispr),abs(rho_rnai),na.rm=TRUE)]
dat[,in_both_lists:=!is.na(rho_crispr) & !is.na(rho_rnai)]
dat[,passes_cut:=(!is.na(rho_crispr) & abs(rho_crispr)>codep_cut) |
                 (!is.na(rho_rnai) & abs(rho_rnai)>codep_cut)]
dat[,complex_member:=gene %chin% complex_members]
# Only meaningful where both screens list the gene; sign agreement is otherwise NA.
dat[in_both_lists==TRUE,same_direction:=sign(rho_crispr)==sign(rho_rnai)]
setorder(dat,-max_abs_rho)
dat[,rank:=.I]

write_tsv(dat,"codependency_merged.tsv")

summ <- data.table(
  n_genes=nrow(dat),
  n_crispr_only=sum(!is.na(dat$rho_crispr) & is.na(dat$rho_rnai)),
  n_rnai_only=sum(is.na(dat$rho_crispr) & !is.na(dat$rho_rnai)),
  n_both=sum(dat$in_both_lists),
  n_above_cut=sum(dat$passes_cut),
  threshold=codep_cut,
  n_complex_members=sum(dat$complex_member),
  n_same_direction=sum(dat$same_direction,na.rm=TRUE),
  n_both_directional=sum(!is.na(dat$same_direction)))
write_tsv(summ,"codependency_summary.tsv")
print(summ)

# The complex members that appear, which is the claim the network figure makes.
members <- dat[complex_member==TRUE,.(rank,gene,rho_crispr,rho_rnai,in_both_lists,passes_cut)]
if (nrow(members)) {
  write_tsv(members,"codependency_complex_members.tsv")
  cat("\nComplex members among the co-dependencies:\n"); print(members)
} else {
  cat("\nNo gene from complex_members appears in either top-100 list.\n")
}

both <- dat[in_both_lists==TRUE]
if (nrow(both)>=10L) {
  r <- spearman(both$rho_crispr,both$rho_rnai)
  if (!is.null(r)) {
    r[,comparison:="codependency_CRISPR_vs_RNAi"]
    write_tsv(r,"codependency_screen_agreement.tsv")
    cat("Agreement over genes listed in both screens: rho =",round(r$rho,3),
        "over",r$n,"genes\n")
  }
} else {
  cat("Only",nrow(both),"genes appear in both lists; screen agreement not computed.\n")
}
