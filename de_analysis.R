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

# volcano plot
volcano_plot <- EnhancedVolcano(de_results,
                                lab = rownames(de_results),
                                x = "avg_log2FC",
                                y = "p_val_adj",
                                pCutoff = 0.05,
                                pointSize = 3.0,
                                labSize = 6.0)
save_plot(volcano_plot, "de_volcano_plot")

# copy de_results
de_copy <- de_results

# visualization transformation
de_copy$log10_p_val <- -log10(de_results$p_val)
de_copy$log10_p_val_adj <- -log10(de_results$p_val_adj)

# selection
significant_genes <- rownames(de_copy[de_copy$p_val_adj < 0.05 & abs(de_copy$avg_log2FC) > 0.5, ])
print(paste("Number of significant genes:", length(significant_genes)))

# heatmap with ordered columns
top_genes <- head(significant_genes, 50)
heatmap_data <- as.matrix(GetAssayData(seurat_VL, layer = "data")[top_genes, ])

# check if any genes have all NA or zero values
genes_to_keep <- rowSums(!is.na(heatmap_data) & heatmap_data != 0) > 0
heatmap_data <- heatmap_data[genes_to_keep, ]

# scale the data
heatmap_data_scaled <- t(scale(t(heatmap_data)))

# create annotation data frame
annotation_col <- data.frame(
  status = seurat_VL$status,
  celltype = seurat_VL$celltype
)

rownames(annotation_col) <- colnames(heatmap_data)

# order the columns by status and then by celltype
order_cols <- order(annotation_col$status, annotation_col$celltype)
heatmap_data_scaled <- heatmap_data_scaled[, order_cols]
annotation_col <- annotation_col[order_cols, ]

annotation_col$status <- factor(annotation_col$status)

# check the levels of the status factor
print(levels(annotation_col$status))

# define colors for status and celltype
status_colors <- c("Healthy" = "#20B2AA", "Duchenne" = "#FF6347")
celltype_colors <- my_celltype_palette
names(celltype_colors) <- unique(annotation_col$celltype)

ann_colors <- list(
  status = status_colors,
  celltype = celltype_colors
)

# create heatmap
heatmap_plot <- pheatmap(heatmap_data_scaled,
                         annotation_col = annotation_col,
                         annotation_colors = ann_colors,
                         show_rownames = TRUE,
                         show_colnames = FALSE,
                         main = "Top Differentially Expressed Genes",
                         fontsize_row = 8,
                         cluster_cols = FALSE,
                         clustering_method = "ward.D2")
  
save_plot(heatmap_plot, "de_heatmap")

