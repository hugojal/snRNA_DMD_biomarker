# Correlation Matrix
expression_matrix <- GetAssayData(seurat_VL, layer = "data")
cor_matrix <- cor(as.matrix(expression_matrix))

all_genes <- rownames(seurat_VL)

# Matrix for AQP4
if ("AQP4" %in% all_genes) {
  all_expression <- GetAssayData(seurat_VL, layer = "data")
  
  # calculate correlation of AQP4 with all other genes
  aqp4_cor <- cor(as.vector(all_expression["AQP4",]), t(as.matrix(all_expression)))
  
  names(aqp4_cor) <- all_genes
  
  # sort and get top correlated genes
  top_correlated_4 <- sort(aqp4_cor, decreasing = TRUE)[1:100] 
  
  print("Top 100 genes correlated with AQP4:")
  print(top_correlated_4)
  
  # create a df
  top_correlated_4_df <- data.frame(
    Gene = names(top_correlated_4),
    Correlation = top_correlated_4
  )
  print(top_correlated_4_df)
  
  write.csv(top_correlated_4_df, "AQP4_top_correlated_genes.csv", row.names = FALSE)
  
} else {
  print("AQP4 gene not found in the dataset")
}

# Matrix for AQP1
if ("AQP1" %in% all_genes) {
  all_expression <- GetAssayData(seurat_VL, layer = "data")
  
  # calculate correlation of AQP4
  aqp1_cor <- cor(as.vector(all_expression["AQP1",]), t(as.matrix(all_expression)))
  
  names(aqp1_cor) <- all_genes
  
  # sort
  top_correlated_1 <- sort(aqp1_cor, decreasing = TRUE)[1:100]  
  
  print("Top 100 genes correlated with AQP1:")
  print(top_correlated_1)
  
  # create df
  top_correlated_1_df <- data.frame(
    Gene = names(top_correlated_1),
    Correlation = top_correlated_1
  )
  print(top_correlated_1_df)
  
  write.csv(top_correlated_1_df, "AQP1_top_correlated_genes.csv", row.names = FALSE)
  
} else {
  print("AQP1 gene not found in the dataset")
}

# Pathway enrichment of correlated genes
top_correlated_4 <- top_correlated_4_df$Gene[-1]

aqp4_ego <- enrichGO(gene = top_correlated_4,
                     OrgDb = org.Hs.eg.db,
                     keyType = "SYMBOL",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05)

aqp4_bar_plot <- barplot(aqp4_ego, showCategory = 20, main = "Enrichment Analysis of Top Correlated Genes to AQP4")
aqp4_bar_plot <- aqp4_bar_plot + ggtitle("Enrichment Analysis of Top Correlated Genes to AQP4") 
save_plot(aqp4_bar_plot, "Enrichment Analysis of Top Correlated Genes to AQP4")

top_correlated_1 <- top_correlated_1_df$Gene[-1]

aqp1_ego <- enrichGO(gene = top_correlated_1,
                     OrgDb = org.Hs.eg.db,
                     keyType = "SYMBOL",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05)

aqp1_bar_plot <- barplot(aqp1_ego, showCategory = 20, main = "Enrichment Analysis of Top Correlated Genes to AQP1")
aqp1_bar_plot <- aqp1_bar_plot + ggtitle("Enrichment Analysis of Top Correlated Genes to AQP1") 
save_plot(aqp1_bar_plot, "Enrichment Analysis of Top Correlated Genes to AQP1")
