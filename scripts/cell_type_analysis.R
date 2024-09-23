# cell type proportion analysis
cell_proportions <- prop.table(table(seurat_VL$celltype, seurat_VL$status), margin = 2)

cell_proportions_df <- as.data.frame(cell_proportions)
cell_proportions_df$Var2 <- factor(cell_proportions_df$Var2, levels = c("Healthy", "Duchenne"))

# create contingency table of cell counts (myofibers vs adipocytes) by condition
contingency_table <- table(seurat_VL$celltype, seurat_VL$status)

test_result <- chisq.test(contingency_table)

print(test_result)

cell_proportion_plot <- ggplot(cell_proportions_df, aes(x = Var2, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal(base_size = 14) +
  labs(x = "Condition", y = "Proportion", fill = "Cell Type") +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold")
  ) +
  geom_signif(
    comparisons = list(c("Healthy", "Duchenne")),
    annotation = "p < 2.2e-16",
    y_position = max(cell_proportions_df$Freq) + 0.05,
    tip_length = 0.02,
    textsize = 5
  )

save_plot(cell_proportion_plot, "cell_props_plot")

# cell type-specific de analysis
Idents(seurat_VL) <- "status"
cell_types <- unique(seurat_VL$celltype)
for (ct in cell_types) {
  ct_de <- FindMarkers(seurat_VL[, seurat_VL$celltype == ct],
                       ident.1 = "Healthy", ident.2 = "Duchenne")
  print(paste("Top DE genes for", ct))
  print(head(ct_de))
}
