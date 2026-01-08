if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("enrichplot")
BiocManager::install("org.Hs.eg.db")

library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(org.Hs.eg.db)

set.seed(123)

gene_symbols <- c("TP53", "BRCA1", "EGFR", "KRAS", "MYC",
                  "PTEN", "AKT1", "PIK3CA", "NRAS", "BRAF",
                  "CDK4", "CDKN2A", "RB1", "ATM", "MDM2",
                  "VEGFA", "HIF1A", "MTOR", "TSC1", "TSC2")

gene_entrez <- bitr(gene_symbols,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)

print(gene_entrez)

kegg_enrich <- enrichKEGG(gene = gene_entrez$ENTREZID,
                          organism = 'hsa',
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)

print(head(kegg_enrich))

kegg_results <- as.data.frame(kegg_enrich)
print(kegg_results)

if (nrow(kegg_results) > 0) {

  p1 <- dotplot(kegg_enrich, showCategory=10, title="KEGG Pathway Enrichment")
  print(p1)

  p2 <- barplot(kegg_enrich, showCategory=10, title="Top 10 Enriched KEGG Pathways")
  print(p2)

  if (nrow(kegg_results) >= 2) {
    p3 <- cnetplot(kegg_enrich, showCategory=5)
    print(p3)
  }

  if (nrow(kegg_results) >= 2) {
    kegg_enrich2 <- pairwise_termsim(kegg_enrich)
    p4 <- emapplot(kegg_enrich2, showCategory=10)
    print(p4)
  }

} else {
  print("No significant KEGG pathways found with current cutoffs")
}

gene_list <- rnorm(length(gene_entrez$ENTREZID), mean=0, sd=2)
names(gene_list) <- gene_entrez$ENTREZID
gene_list <- sort(gene_list, decreasing=TRUE)

#Ranked gene list (top 10):
print(head(gene_list, 10))

gsea_kegg <- gseKEGG(geneList = gene_list,
                     organism = 'hsa',
                     minGSSize = 10,
                     maxGSSize = 500,
                     pvalueCutoff = 0.05,
                     verbose = FALSE)

print(head(gsea_kegg))

if (nrow(as.data.frame(gsea_kegg)) > 0) {
  p5 <- ridgeplot(gsea_kegg, showCategory=10)
  print(p5)

  p6 <- gseaplot2(gsea_kegg, geneSetID=1:3, pvalue_table=TRUE)
  print(p6)
}

if (nrow(kegg_results) > 0) {
  write.csv(kegg_results, "KEGG_enrichment_results.csv", row.names=FALSE)
}
