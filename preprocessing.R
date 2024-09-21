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
