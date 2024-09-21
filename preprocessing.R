# load and preprocess data
seurat <- readRDS("C:\\Users\\hugoj\\Desktop\\BioHacks\\3DMD_VL 3HealthyVL_TA.rds")
Idents(seurat) <- "ID"
seurat_VL <- subset(seurat, idents = grep("_VL$", Idents(seurat), value = TRUE))

# correct celltype names
seurat_VL$celltype <- recode(seurat_VL$celltype,
                             "myofibers" = "Myofibers",
                             "macrophage" = "Macrophage",
                             "Endothelial" = "Endothelial",
                             "Satellite cells" = "Satellite Cells",
                             "FAPS" = "FAPs",
                             "Smooth muscle" = "Smooth Muscle",
                             "t and b cells" = "T and B Cells",
                             "adipocytes" = "Adipocytes"
)

# correct status names
seurat_VL$status <- str_to_title(seurat_VL$status)

# qc metrics
qc_plot <- VlnPlot(seurat_VL, features = c("nFeature_RNA", "nCount_RNA"), cols = custom_palette, ncol = 2, combine = TRUE)
qc_plot <- qc_plot + ggtitle("Quality Metrics")
qc_plot <-  qc_plot + plot_annotation(tag_levels = 'a')
save_plot(qc_plot, "qc_metrics")

# clustering
seurat_VL <- FindNeighbors(seurat_VL, dims = 1:10)
seurat_VL <- FindClusters(seurat_VL, resolution = 0.5)
