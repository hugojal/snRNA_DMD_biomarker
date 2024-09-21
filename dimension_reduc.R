# pca
seurat_VL <- RunPCA(seurat_VL, features = VariableFeatures(object = seurat_VL))
print(seurat_VL[["pca"]], dims = 1:5, nfeatures = 5)
pca_plot <- DimPlot(seurat_VL, reduction = "pca", group.by = "status") +
  ggtitle("PCA Plot") +
  scale_color_manual(values = my_status_palette) +
  theme(plot.title = element_text(hjust = 0.5))
save_plot(pca_plot, "PCA_plot", width = 8, height = 8)

pca_elbow_plot <- ElbowPlot(seurat_VL, ndims = 50) +
  ggtitle("Elbow Plot for PCA") +
  theme(plot.title = element_text(hjust = 0.5))
save_plot(pca_elbow_plot, "PCA_elbow_plot")

# umap
seurat_VL <- RunUMAP(seurat_VL, dims = 1:10)
umap_plot.1 <- DimPlot(seurat_VL, reduction = "umap", group.by = "celltype", 
                       cols = my_celltype_palette, combine = TRUE)
umap_plot.1 <- umap_plot.1 + ggtitle("UMAP by Cell Type")
umap_plot.2 <- DimPlot(seurat_VL, reduction = "umap", group.by = c("status"), 
                       cols = my_status_palette, combine = TRUE)
umap_plot.2 <- umap_plot.2 + ggtitle("UMAP by Condition")
save_plot(umap_plot.1, "umap_visualization1", width = 8, height = 8)
save_plot(umap_plot.2, "umap_visualization2", width = 8, height = 8)
