if (!require("BiocManager", quietly = TRUE))
     install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")

library(clusterProfiler)
library(org.Hs.eg.db)

my_genes <- c("TP53", "EGFR", "TNF", "IL6", "VEGFA", "MYC", "GAPDH")

gene_entrez <- bitr(my_genes, 
                    fromType = "SYMBOL", 
                    toType   = "ENTREZID", 
                    OrgDb    = org.Hs.eg.db)


print(head(gene_entrez))

kegg_results <- enrichKEGG(
    gene         = gene_entrez$ENTREZID,
    organism     = 'hsa',
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
)

head(kegg_results@result)

dotplot(kegg_results, showCategory=10, title="Top 10 Enriched KEGG Pathways")

barplot(kegg_results, showCategory=10, title="Top 10 Enriched KEGG Pathways")
