# Pathway Enrichment Analysis
ego <- enrichGO(gene = significant_genes,
                OrgDb = org.Hs.eg.db,
                keyType = "SYMBOL",
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05)

ego_bar_plot <- barplot(ego, showCategory = 20, title = "Top Enriched Pathways")
save_plot(ego_bar_plot, "egobar")

# Get the top 10 GO terms
top_terms <- head(ego@result$Description, 10)

# Function to clean up and format GO terms: replace "." with space and capitalize words
clean_terms <- function(terms) {
  terms_clean <- str_replace_all(terms, "\\.", " ")  # Replace dots with spaces
  terms_clean <- str_to_title(terms_clean)  # Capitalize first letter of each word
  return(terms_clean)
}

# Clean up the top terms
top_terms_clean <- clean_terms(top_terms)

# Function to get top N genes for a term
get_top_n_genes <- function(term_genes, n = 10) {
  gene_list <- unlist(strsplit(term_genes, "/"))
  gene_lfc <- abs(lfc[gene_list])
  gene_lfc <- gene_lfc[!is.na(gene_lfc)]  # Remove NA values
  if(length(gene_lfc) == 0) return(character(0))  # Return empty if no genes have LFC
  top_genes <- names(sort(gene_lfc, decreasing = TRUE)[1:min(n, length(gene_lfc))])
  return(top_genes)
}

# Get top 10 genes for each term
top_genes_per_term <- lapply(ego@result$geneID[1:10], get_top_n_genes)

# Flatten the list of top genes
all_top_genes <- unique(unlist(top_genes_per_term))

# Create the matrix for the chord diagram
mat <- matrix(0, nrow = length(all_top_genes), ncol = length(top_terms_clean))
rownames(mat) <- all_top_genes
colnames(mat) <- top_terms_clean  # Use cleaned top terms as column names

for (i in 1:length(top_terms_clean)) {
  term_genes <- top_genes_per_term[[i]]
  mat[term_genes, i] <- 1
}

# Create chord_data
chord_data <- data.frame(mat)
chord_data$logFC <- lfc[all_top_genes]

# Remove any rows with NA in logFC (if any)
chord_data <- chord_data[!is.na(chord_data$logFC), ]

# Now create the chord diagram
ego_chord_plot <- GOChord(data = chord_data, 
                          gene.order = 'logFC', 
                          gene.space = 0.25, 
                          gene.size = 4,
                          ribbon.col = custom_palette[1:10],
                          process.label = 6,
                          limit = c(1, 1),
                          space = 0.01)

# Add title
title <- "Top 10 Genes per GO Term"
ego_chord_plot <- ego_chord_plot + ggtitle(title)

# Display the plot
print(ego_chord_plot)

# Save the plot
ggsave("GO_chord_plot_filtered.jpeg", ego_chord_plot, width = 16, height = 14, dpi = 300)
