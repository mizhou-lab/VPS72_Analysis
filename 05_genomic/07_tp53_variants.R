# VPS72 expression by individual recurrent TP53 protein change:
#   VPS72_z ~ carries this variant + age + sex + cancer type
# The reference group is every other sample, including carriers of a different
# TP53 variant, so a coefficient is not a comparison against TP53 wild-type.
source("00_config.R")

master <- fread(require_input(inputs$master))

mc3 <- read_mc3()
prot <- find_col(mc3,c("HGVSp_Short","HGVSp","Protein Change","Protein.Change","amino_acid_change","AAChange"))
setnames(mc3,prot,"protein_change")
tp53 <- mc3[gene=="TP53"]
if (!nrow(tp53)) stop("No non-silent TP53 calls found")
tp53[,protein_change:=gsub("^p\\.","",trimws(as.character(protein_change)))]
tp53 <- unique(tp53[!is.na(protein_change) & nzchar(protein_change),.(patient,protein_change)])
rm(mc3); gc()

recurrent <- tp53[,.N,by=protein_change][N>=min_mutated]$protein_change
if (!length(recurrent)) stop("No TP53 variant reaches ",min_mutated," carriers")

carrier <- dcast(tp53[protein_change %chin% recurrent],patient~protein_change,
                 fun.aggregate=length,fill=0L)
carrier <- merge(master[,.(patient,VPS72_z,age,sex,cancer)],carrier,by="patient",all.x=TRUE)
# Patients absent from the TP53 table are non-carriers, not missing.
for (v in recurrent) carrier[,(v):=as.integer(fifelse(is.na(get(v)),0L,get(v))>0L)]

covars <- carrier[,.(VPS72_z,age,sex,cancer)]
fit_variant <- function(v) {
  d <- data.table(carries=carrier[[v]],covars)
  d <- d[complete.cases(d)]
  if (uniqueN(d$carries)<2L) return(NULL)
  m <- try(lm(VPS72_z~carries+age+sex+factor(cancer),data=d),silent=TRUE)
  if (inherits(m,"try-error")) return(NULL)
  tid <- as.data.table(tidy(m,conf.int=TRUE))[term=="carries"]
  if (!nrow(tid)) return(NULL)
  data.table(gene="TP53",variant=v,n_carriers=sum(d$carries==1L),
             beta=tid$estimate,se=tid$std.error,
             ci_low=tid$conf.low,ci_high=tid$conf.high,t=tid$statistic,p=tid$p.value)
}

res <- rbindlist(lapply(recurrent,fit_variant))
if (!nrow(res)) stop("No variant model returned a coefficient")
res[,FDR:=p.adjust(p,method="BH")]
setorder(res,p)

write_tsv(res,"tp53_variant_association.tsv")
cat("Recurrent variants tested:",length(recurrent),
    "; at FDR <",fdr_cut,":",sum(res$FDR<fdr_cut,na.rm=TRUE),"\n")
