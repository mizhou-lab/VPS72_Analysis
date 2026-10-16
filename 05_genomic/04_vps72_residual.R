# VPS72 expression not accounted for by its complex partners.
# Uses srcap_genes and tip60_genes from 00_config.R, which exclude VPS72 itself.
source("00_config.R")

annot <- read_annot()

expr <- fread(require_input(inputs$expr))
setnames(expr,1,"gene_id")
expr[,symbol:=toupper(sub("\\|.*$","",gene_id))]
sub <- expr[symbol %chin% unique(c(srcap_genes,tip60_genes,"VPS72"))]
if (!nrow(sub)) stop("No complex genes found in the expression matrix")

long <- melt(sub,id.vars=c("gene_id","symbol"),
             measure.vars=setdiff(names(sub),c("gene_id","symbol")),
             variable.name="sample",value.name="expr")
long[,sample16:=tcga16(sample)]
long[,expr:=suppressWarnings(as.numeric(expr))]
wide <- dcast(long,sample16~symbol,value.var="expr",fun.aggregate=mean)
wide <- merge(wide,annot,by="sample16")
rm(expr,sub,long); gc()

srcap_ok <- intersect(srcap_genes,names(wide))
tip60_ok <- intersect(tip60_genes,names(wide))
if (!"VPS72" %in% names(wide)) stop("VPS72 missing from the expression matrix")
if (length(srcap_ok)<2L || length(tip60_ok)<2L) stop("Too few complex members recovered to build scores")

wide[,SRCAP_score:=rowMeans(.SD,na.rm=TRUE),.SDcols=srcap_ok]
wide[,TIP60_score:=rowMeans(.SD,na.rm=TRUE),.SDcols=tip60_ok]
dat <- wide[is.finite(VPS72) & is.finite(SRCAP_score) & is.finite(TIP60_score) & nzchar(cancer)]
if (!nrow(dat)) stop("No samples with both complex scores and VPS72")

# Cancer type is in the model, so the residual is a deviation from what is
# expected for a tumour of that lineage with those complex scores.
fit <- lm(VPS72~SRCAP_score+TIP60_score+factor(cancer),data=dat)
dat[,VPS72_resid:=residuals(fit)]
write_tsv(dat[,.(sample16,cancer,VPS72,SRCAP_score,TIP60_score,VPS72_resid)],
          "vps72_residual_by_sample.tsv")

# Residuals centred on their grand mean, then each cancer type tested against
# zero. Cancer types with fewer than three samples are dropped.
dat[,rel_resid:=VPS72_resid-mean(VPS72_resid,na.rm=TRUE)]
by_cancer <- dat[,{
  if (.N<3L) NULL else {
    tt <- t.test(rel_resid,mu=0)
    .(n=.N,mean_resid=mean(rel_resid),sd_resid=sd(rel_resid),
      ci_low=tt$conf.int[1L],ci_high=tt$conf.int[2L],p=tt$p.value)
  }},by=cancer]
if (!nrow(by_cancer)) stop("No cancer type has three or more samples")
by_cancer[,FDR:=p.adjust(p,method="BH")]
by_cancer[,direction:=fifelse(mean_resid>0,"higher_than_expected","lower_than_expected")]
setorder(by_cancer,p)

write_tsv(by_cancer,"vps72_residual_by_cancer.tsv")
cat("Samples:",nrow(dat),"; adjusted R-squared:",round(summary(fit)$adj.r.squared,3),
    "; cancer types at FDR <",fdr_cut,":",sum(by_cancer$FDR<fdr_cut,na.rm=TRUE),
    "of",nrow(by_cancer),"\n")
