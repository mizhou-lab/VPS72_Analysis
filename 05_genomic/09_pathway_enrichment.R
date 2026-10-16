# GO BP, KEGG and Hallmark over-representation among the genes from 05.
# Background is every gene tested in 05, not the whole genome, so the
# comparison stays within genes that had enough carriers to be modelled.
source("00_config.R")
suppressPackageStartupMessages({library(clusterProfiler);library(org.Hs.eg.db);library(msigdbr)})

assoc <- fread(require_input(inputs$assoc))
hits <- assoc[FDR<=enrich_fdr]$gene
if (length(hits)<10L) stop("Only ",length(hits)," genes at FDR <= ",enrich_fdr,
                           "; too few for over-representation testing")
cat("genes at FDR <=",enrich_fdr,":",length(hits),"of",nrow(assoc),"tested\n")

to_entrez <- function(g) suppressMessages(
  bitr(g,fromType="SYMBOL",toType="ENTREZID",OrgDb=org.Hs.eg.db)$ENTREZID)
entrez_hits <- to_entrez(hits)
entrez_bg <- to_entrez(assoc$gene)
if (!length(entrez_hits)) stop("No hit gene symbol mapped to an Entrez ID")

ratio <- function(v) vapply(strsplit(as.character(v),"/"),
                            function(z) as.numeric(z[1L])/as.numeric(z[2L]),numeric(1))
tidy_enrich <- function(x,label) {
  if (is.null(x)) return(NULL)
  d <- as.data.table(as.data.frame(x))
  if (!nrow(d)) return(NULL)
  d[,fold_enrichment:=ratio(GeneRatio)/ratio(BgRatio)]
  d[,.(source=label,id=ID,description=Description,gene_ratio=GeneRatio,bg_ratio=BgRatio,
       fold_enrichment,count=Count,p=pvalue,FDR=p.adjust,gene_ids=geneID)]
}

go <- enrichGO(gene=entrez_hits,universe=entrez_bg,OrgDb=org.Hs.eg.db,keyType="ENTREZID",
               ont="BP",pAdjustMethod="BH",pvalueCutoff=0.05,qvalueCutoff=0.20,readable=TRUE)
# Redundant parent terms collapsed.
if (!is.null(go) && nrow(as.data.frame(go))) go <- simplify(go,cutoff=0.7,by="p.adjust",select_fun=min)

kegg <- enrichKEGG(gene=entrez_hits,universe=entrez_bg,organism="hsa",
                   pAdjustMethod="BH",pvalueCutoff=0.20,qvalueCutoff=0.20)
if (!is.null(kegg) && nrow(as.data.frame(kegg))) kegg <- setReadable(kegg,OrgDb=org.Hs.eg.db,keyType="ENTREZID")

msig <- as.data.frame(msigdbr(species="Homo sapiens",category="H"))
hall <- enricher(gene=entrez_hits,universe=entrez_bg,
                 TERM2GENE=unique(msig[,c("gs_name","entrez_gene")]),
                 pAdjustMethod="BH",pvalueCutoff=0.20,qvalueCutoff=0.20)

out <- rbindlist(list(tidy_enrich(go,"GO_BP"),tidy_enrich(kegg,"KEGG"),
                      tidy_enrich(hall,"Hallmark")),fill=TRUE)
if (!nrow(out)) stop("No term reached the cutoffs in any of the three collections")
setorder(out,source,p)

write_tsv(out,"pathway_enrichment.tsv")
print(out[,.N,by=source])
