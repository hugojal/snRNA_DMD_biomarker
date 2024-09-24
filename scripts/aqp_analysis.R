# Aquaporin analysis
aquaporin_genes <- c("AQP1", "AQP2", "AQP3", "AQP4", "AQP5", "AQP6", 
                     "AQP7", "AQP8", "AQP9", "AQP10", "AQP11", "AQP12")
aquaporin_present <- aquaporin_genes[aquaporin_genes %in% rownames(seurat_VL)]

aqp_data <- FetchData(seurat_VL, vars = aquaporin_present)
aqp_data$celltype <- seurat_VL$celltype  # Add cell type column
aqp_data$status <- seurat_VL$status

# reshape data into long
aqp_data_long <- tidyr::pivot_longer(aqp_data, cols = all_of(aquaporin_present), names_to = "features", values_to = "expression")

aqp_data_long$status <- factor(aqp_data_long$status, levels = c("Healthy", "Duchenne"))

# filter out features where all values are the same
aqp_data_long <- aqp_data_long %>% dplyr::group_by(features) %>%
  dplyr::filter(var(expression) > 0)

# visualizations
aqp_vln_type_plot <- ggplot(aqp_data_long, aes(x = celltype, y = expression, fill = celltype)) +
  geom_violin(scale = "width") +
  theme_minimal() +
  facet_wrap(~ features, ncol = 1, strip.position = "right") +  # Align plots vertically
  theme(
    strip.text = element_text(size = 12, face = "bold"),  # Increase gene label size
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for readability
    axis.title.y = element_blank(),  # Remove y-axis title
    axis.title.x = element_text(size = 12),  # Increase size of x-axis title
    panel.grid = element_blank(),  # Remove grid lines
    panel.spacing.y = unit(0.5, "lines")  # Reduce vertical spacing
  ) +
  xlab("Cell Type") +
  ylab("Expression Level")

ggsave("aqp_vln_type_plot.png", plot = aqp_vln_type_plot, height = 14, width = 7)

aqp_vln_stat_plot <- ggplot(aqp_data_long, aes(x = status, y = expression, fill = celltype)) +
  geom_violin(scale = "width") +
  theme_minimal() +
  facet_wrap(~ features, ncol = 1, strip.position = "right") +  # Align plots vertically
  theme(
    strip.text = element_text(size = 12, face = "bold"),  # Increase gene label size
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for readability
    axis.title.y = element_blank(),  # Remove y-axis title
    axis.title.x = element_text(size = 12),  # Increase size of x-axis title
    panel.grid = element_blank(),  # Remove grid lines
    panel.spacing.y = unit(0.5, "lines")  # Reduce vertical spacing
  ) +
  xlab("Status") +
  ylab("Expression Level")
save_plot(aqp_vln_stat_plot, "aquaporing_expression_violin_by_status", height = 14, width = 7)

# AQP in UMAP
library(Seurat)
library(ggplot2)

# create list for plots
aquaporin_plots <- list()

for (gene in aquaporin_genes) {
  plot <- FeaturePlot(seurat_VL, features = gene, reduction = "umap") +
    ggtitle(gene) +
    theme(plot.title = element_text(hjust = 0.5))
  aquaporin_plots[[gene]] <- plot
}

# combine plots
combined_plot <- cowplot::plot_grid(plotlist = aquaporin_plots, ncol = 3)

ggsave("aquaporin_expression_umap.png", combined_plot, width = 15, height = 5 * ceiling(length(aquaporin_genes) / 3))
save_plot(combined_plot, "aquaporin_expression_umap.png")

# AQP Healthy vs DMD
for (gene in aquaporin_genes) {
  # Create separate plots for Duchenne and Healthy
  duchenne_plot <- FeaturePlot(seurat_VL[, seurat_VL$status == "duchenne"], 
                               features = gene, 
                               reduction = "umap") +
    ggtitle("Duchenne") +
    theme(plot.title = element_text(hjust = 0.5))
  
  healthy_plot <- FeaturePlot(seurat_VL[, seurat_VL$status == "healthy"], 
                              features = gene, 
                              reduction = "umap") +
    ggtitle("Healthy") +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Combine plots side by side
  combined_plot <- duchenne_plot + healthy_plot +
    plot_layout(ncol = 2) +
    plot_annotation(title = gene,
                    theme = theme(plot.title = element_text(hjust = 0.5)))
  
  # Save the combined plot
  ggsave(paste0(gene, "_split_umap.png"), combined_plot, width = 12, height = 5)
}
