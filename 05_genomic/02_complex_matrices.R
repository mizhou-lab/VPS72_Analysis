# Gene x cancer-type matrices for the SRCAP/TIP60/INO80 genes:
# median expression z, mean GISTIC2 copy number, amplification frequency.
source("00_config.R")

annot <- read_annot()

expr <- fread(require_input(inputs$expr))
setnames(expr,1,"gene_id")
expr[,symbol:=toupper(sub("\\|.*$","",gene_id))]
panel <- expr[symbol %chin% panel_genes]
if (!nrow(panel)) stop("No panel genes found in the expression matrix")
absent <- setdiff(panel_genes,panel$symbol)
if (length(absent)) cat("not in expression file:",paste(absent,collapse=", "),"\n")

# z-scored per gene across all samples first, so the matrix is comparable
# between genes, then summarised to a median per cancer type.
long <- melt(panel,id.vars=c("gene_id","symbol"),
             measure.vars=setdiff(names(panel),c("gene_id","symbol")),
             variable.name="sample",value.name="expr")
long[,sample16:=tcga16(sample)]
long <- merge(long,annot,by="sample16")
long[,expr:=suppressWarnings(as.numeric(expr))]
long[,expr_z:=as.numeric(scale(expr)),by=symbol]

expr_mat <- dcast(long[,.(value=median(expr_z,na.rm=TRUE)),by=.(gene=symbol,cancer)],
                  gene~cancer,value.var="value")
write_tsv(expr_mat,"complex_expression_medianZ.tsv")
rm(expr,panel,long); gc()

cnv <- fread(require_input(inputs$gistic))
setnames(cnv,1,"symbol")
cnv[,symbol:=toupper(symbol)]
cnv_panel <- cnv[symbol %chin% panel_genes]
if (!nrow(cnv_panel)) stop("No panel genes found in the GISTIC2 table")

# Both copy-number summaries come from one melt.
cl <- melt(cnv_panel,id.vars="symbol",measure.vars=setdiff(names(cnv_panel),"symbol"),
           variable.name="sample",value.name="cnv")
cl[,sample16:=tcga16(sample)]
cl <- merge(cl,annot,by="sample16")
cl[,cnv:=suppressWarnings(as.numeric(cnv))]

write_tsv(dcast(cl[,.(value=mean(cnv,na.rm=TRUE)),by=.(gene=symbol,cancer)],
                gene~cancer,value.var="value"),"complex_cnv_mean.tsv")
write_tsv(dcast(cl[,.(value=100*mean(cnv>=amp_gistic,na.rm=TRUE)),by=.(gene=symbol,cancer)],
                gene~cancer,value.var="value"),"complex_ampfreq.tsv")
cat("Panel genes:",uniqueN(cl$symbol),"; cancer types:",uniqueN(cl$cancer),"\n")
