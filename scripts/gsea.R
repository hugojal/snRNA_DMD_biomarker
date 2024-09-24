# data prep
gene_list <- de_results$avg_log2FC
names(gene_list) <- rownames(de_results)
gene_list <- sort(gene_list, decreasing = TRUE)

# GSEA
gsea_result <- gseGO(geneList = gene_list,
                     OrgDb = org.Hs.eg.db,
                     keyType = "SYMBOL",
                     ont = "BP",
                     pvalueCutoff = 0.05,
                     verbose = FALSE)

gsea_plot <- gseaplot2(gsea_result, geneSetID = 1:5, pvalue_table = TRUE)
save_plot(gsea_plot, "gsea_plot")
save_csv(as.data.frame(gsea_result), "gsea_results")
