# top markers for each cluster
markers <- FindAllMarkers(seurat_VL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
save_csv(top10, "top10_markers_per_cluster")

# de analysis
Idents(seurat_VL) <- "status"
de_results <- FindMarkers(seurat_VL, ident.1 = "Duchenne", ident.2 = "Healthy")

# fdr
p_values <- de_results$p_val
bh_adjusted <- p.adjust(p_values, method = "BH")
de_results$FDR_BH <- bh_adjusted

# set identity to celltype
Idents(seurat_VL) <- "celltype"

# calculate average expression for each cell type
avg_expr <- AverageExpression(seurat_VL, assays = "RNA")$RNA

# function to find the cell type with max expression for each gene
find_max_celltype <- function(gene) {
  celltype_expr <- avg_expr[gene,]
  return(names(celltype_expr)[which.max(celltype_expr)])
}

# add the cell type with highest expression for each gene
de_results$prominent_celltype <- sapply(rownames(de_results), find_max_celltype)

# save the results
save_csv(de_results, "differential_expression_analysis")
