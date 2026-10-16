# Overlap of the three perturbation programmes, split into seven mutually
# exclusive partitions, then tested for Hallmark over-representation.
# This is the dependency logic: which VPS72 programmes are shared with SRCAP,
# with H2A.Z.1, with both, or with neither.
source("00_config.R")
suppressPackageStartupMessages({library(clusterProfiler);library(msigdbr)})

deg <- fread(require_input(inputs$deg))
needed <- c("VPS72","SRCAP","H2AZ1")
if (!all(needed %chin% deg$contrast)) stop("Missing contrast: ",
  paste(setdiff(needed,deg$contrast),collapse=", "))

gene_set <- function(ct,dir) deg[contrast==ct & direction==dir]$gene

sets <- rbindlist(lapply(c("Up","Down"),function(dir) {
  v <- gene_set("VPS72",dir); s <- gene_set("SRCAP",dir); h <- gene_set("H2AZ1",dir)
  parts <- make_partitions(v,s,h)
  rbindlist(lapply(names(parts),function(p)
    if (!length(parts[[p]])) NULL else data.table(direction=dir,partition=p,gene=parts[[p]])))
}))
if (!nrow(sets)) stop("No gene fell into any partition; check the significance thresholds")

write_tsv(sets,"regulated_gene_sets.tsv")

counts <- sets[,.(n=.N),by=.(direction,partition)]
# Set sizes per perturbation, alongside the partition counts.
per_contrast <- deg[direction!="NS",.(n=.N),by=.(contrast,direction)]
write_tsv(dcast(counts,partition~direction,value.var="n",fill=0L),"partition_counts.tsv")
write_tsv(dcast(per_contrast,contrast~direction,value.var="n",fill=0L),"contrast_counts.tsv")

# Which partitions imply dependence on which complex member. Stated as a table
# rather than left implicit in the partition names.
logic <- data.table(partition=c("VPS72_only","SRCAP_only","H2AZ1_only","VPS72_SRCAP",
                                "VPS72_H2AZ1","SRCAP_H2AZ1","CORE"))
logic[,`:=`(VPS72_regulated=partition %chin% c("VPS72_only","VPS72_SRCAP","VPS72_H2AZ1","CORE"),
            SRCAP_dependent=partition %chin% c("SRCAP_only","VPS72_SRCAP","SRCAP_H2AZ1","CORE"),
            H2AZ1_dependent=partition %chin% c("H2AZ1_only","VPS72_H2AZ1","SRCAP_H2AZ1","CORE"))]
write_tsv(logic,"dependency_logic.tsv")

# The background is every gene called significant in any of the three
# perturbations, not the whole transcriptome: the question is which partition a
# programme falls into, so the comparison is between partitions of that set.
universe <- unique(deg[direction!="NS"]$gene)
hall <- as.data.table(msigdbr(species="Homo sapiens",category="H"))
t2g <- unique(hall[,.(gs_name,gene_symbol=toupper(gene_symbol))])

enrich_one <- function(dir,part) {
  g <- sets[direction==dir & partition==part]$gene
  if (length(g)<min_set_n) return(NULL)
  e <- tryCatch(enricher(gene=g,universe=universe,TERM2GENE=t2g,
                         pAdjustMethod="BH",minGSSize=min_set_n),error=function(e) NULL)
  if (is.null(e)) return(NULL)
  d <- as.data.table(as.data.frame(e))
  if (!nrow(d)) return(NULL)
  ratio <- function(v) vapply(strsplit(as.character(v),"/"),
                              function(z) as.numeric(z[1L])/as.numeric(z[2L]),numeric(1))
  d[,`:=`(direction=dir,partition=part,n_partition=length(g),
          fold_enrichment=ratio(GeneRatio)/ratio(BgRatio))]
  d[,.(direction,partition,n_partition,term=ID,description=Description,
       gene_ratio=GeneRatio,bg_ratio=BgRatio,fold_enrichment,count=Count,
       p=pvalue,FDR=p.adjust,gene_ids=geneID)]
}

grid <- CJ(dir=c("Up","Down"),part=logic$partition,sorted=FALSE)
enrich <- rbindlist(Map(enrich_one,grid$dir,grid$part))
if (!nrow(enrich)) stop("No partition reached ",min_set_n," genes with an enriched term")
setorder(enrich,direction,partition,p)
write_tsv(enrich,"partition_enrichment.tsv")

cat("Partition sizes:\n"); print(dcast(counts,partition~direction,value.var="n",fill=0L))
cat("Partitions with at least one enriched Hallmark term:",
    uniqueN(enrich[,.(direction,partition)]),"\n")
