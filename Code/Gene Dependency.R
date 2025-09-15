# Gene Dependency 
library(dplyr)
library(tibble)
library(pheatmap)
library(RColorBrewer)
library(tidyr)
library(ggplot2)

# --- Part 1: Analysis using local CSV files ---

# Load the DepMap data from local CSV files
# depmap_data: CRISPR gene effect scores. Rows are cell lines (ModelID), columns are genes.
depmap_data <- read.csv("./Data/depmap/CRISPRGeneEffect_25Q2.csv", row.names=1, check.names=FALSE, stringsAsFactors=FALSE)
# depmap_meta: Metadata for cell lines. Rows are cell lines (ModelID).
depmap_meta <- read.csv("./Data/depmap/Model.csv", row.names=1, check.names=FALSE, stringsAsFactors=FALSE)
essential_list <- read.csv("./Data/depmap/common_essentials.csv", stringsAsFactors=FALSE, header = TRUE)

message("--- Classifying DepMap data by OncotreeLineage---")

# 1. Find common cell lines that exist in both the data and metadata files.
common_cell_lines <- intersect(rownames(depmap_data), rownames(depmap_meta))
message(sprintf("Found %d common cell lines between dependency data and metadata.", length(common_cell_lines)))

# 2. Filter both dataframes to keep only these common cell lines for consistency.
depmap_data_common <- depmap_data[common_cell_lines, ]
depmap_meta_common <- depmap_meta[common_cell_lines, ]

# 3. Merge the disease type information into the gene dependency data.
depmap_data_annotated <- depmap_data_common %>%
  rownames_to_column("ModelID") %>%
  left_join(
    depmap_meta_common %>% 
      select(OncotreeLineage) %>%
      rownames_to_column("ModelID"),
    by = "ModelID"
  ) %>%
  # Move the disease column to the front for clarity
  select(ModelID, OncotreeLineage, everything())
depmap_data_annotated[1:5, 1:5]

# --- Part 1b: Subset for Essential Genes ---
message("\n--- Subsetting data for common essential genes ---")

# 1. Get the list of essential genes from the loaded file.
essential_genes <- essential_list$gene

# 2. Find which of these essential genes are present as columns in our data.
common_essential_genes <- intersect(essential_genes, colnames(depmap_data_annotated))
message(sprintf("Found %d common essential genes in the dependency data.", length(common_essential_genes)))

# 3. Create a new dataframe containing only the essential gene columns, plus identifiers.
depmap_essential_annotated <- depmap_data_annotated %>%
  select(ModelID, OncotreeLineage, all_of(common_essential_genes))

message("Created 'depmap_essential_annotated' with essential genes only.")
print(head(depmap_essential_annotated[, 1:6]))

# --- Part 1c: Calculate Average Essentiality by Lineage ---
message("\n--- Calculating average essentiality score for each gene by lineage ---")

# Group by lineage and calculate the mean score for each essential gene.
# A more negative score indicates higher essentiality.
lineage_essentiality_scores <- depmap_essential_annotated %>%
  group_by(OncotreeLineage) %>%
  summarise(across(all_of(common_essential_genes), ~mean(.x, na.rm = TRUE)))

message("Created 'lineage_essentiality_scores' with average gene scores per lineage.")
print(lineage_essentiality_scores[, 1:6])

# --- Part 1d: Visualize Essentiality Scores with a Heatmap ---
message("\n--- Generating heatmap of lineage essentiality scores ---")

# Create an 'Images' directory if it doesn't exist, similar to other scripts
if (!dir.exists("./Images")) {
  dir.create("./Images")
}

# Prepare the data for pheatmap. It needs to be a numeric matrix with lineages as rownames.
# It's also good practice to filter out lineages with very few cell lines (<5), as their
# averages can be noisy.
lineage_counts <- count(depmap_data_annotated, OncotreeLineage)
lineages_to_keep <- lineage_counts$OncotreeLineage[lineage_counts$n >= 5]

heatmap_matrix <- lineage_essentiality_scores %>%
  filter(OncotreeLineage %in% lineages_to_keep) %>%
  tibble::column_to_rownames("OncotreeLineage") %>%
  as.matrix()

# Generate the heatmap. A lower score (more negative) means more essential.
# We'll use a color scale where blue represents low scores (high essentiality).
pheatmap(
  heatmap_matrix,
  main = "Average Essentiality of Common Essential Genes Across Cancer Lineages",
  fontsize_row = 8,
  show_colnames = FALSE, # Too many genes to display cleanly
  cluster_cols = TRUE,   # Cluster genes by their essentiality profiles
  cluster_rows = TRUE,   # Cluster lineages by their dependency profiles
  scale = "none",        # The scores are already on a comparable scale
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
  filename = "./Images/lineage_essentiality_heatmap.png", # Save to file
  width = 10,
  height = 8
)

message("Heatmap saved to ./Images/lineage_essentiality_heatmap.png")

# --- Part 1e: Find Lineage-Specific Essential Genes (Example: Melanoma) ---
message("\n--- Finding top 5 genes specifically essential to Melanoma ---")

# To find specific dependencies, we calculate a "specificity score":
# specificity_score = (Gene Score in Melanoma) - (Mean Gene Score in Other Lineages)
# A more negative score indicates higher specific essentiality in Melanoma.

melanoma_specificity <- lineage_essentiality_scores %>%
  # Use only lineages with enough cell lines for a fair comparison (from heatmap)
  filter(OncotreeLineage %in% lineages_to_keep) %>%
  # Convert to long format for easier manipulation
  pivot_longer(cols = -OncotreeLineage, names_to = "gene", values_to = "score") %>%
  # Group by gene to compare Melanoma vs. Others
  group_by(gene) %>%
  # Calculate the score in Melanoma and the mean score in all other lineages
  summarise(
    melanoma_score = score[OncotreeLineage == "Skin"],
    other_lineages_mean_score = mean(score[OncotreeLineage != "Skin"], na.rm = TRUE)
  ) %>%
  mutate(specificity_score = melanoma_score - other_lineages_mean_score) %>%
  filter(!is.na(specificity_score)) %>%
  arrange(specificity_score)

top_melanoma_specific_genes <- head(melanoma_specificity, 50)
print(top_melanoma_specific_genes, n = 50)
top_melanoma_specific_genes$gene

# strip the gene names for easier reading
top_melanoma_specific_genes$gene <- gsub("^(.*) \\(.*\\)$", "\\1", top_melanoma_specific_genes$gene)
# export the gene names to a text file
write.table(top_melanoma_specific_genes$gene, file = "./Results/melanoma_specific_genes.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
top_melanoma_specific_genes$gene

# --- Part 1f: Visualize Specificity Scores for Top 5 Genes ---
message("\n--- Visualizing specificity of top 5 Skin-lineage genes ---")

# A grouped bar plot is ideal for comparing the dependency scores directly.

# 1. Take the top 5 genes from the previously calculated table.
# The gene names in `top_melanoma_specific_genes` have already been cleaned.
top_5_for_plot <- head(top_melanoma_specific_genes, 5)

# 2. Reshape the data from a wide to a long format, which is required for ggplot.
plot_data <- top_5_for_plot %>%
  select(gene, melanoma_score, other_lineages_mean_score) %>%
  pivot_longer(
    cols = c("melanoma_score", "other_lineages_mean_score"),
    names_to = "lineage_group",
    values_to = "mean_score"
  ) %>%
  # Make the gene names an ordered factor to preserve the ranking in the plot
  mutate(
    gene = factor(gene, levels = top_5_for_plot$gene),
    lineage_group = recode(lineage_group,
                           "melanoma_score" = "Skin Lineage",
                           "other_lineages_mean_score" = "Other Lineages (Mean)")
  )

# 3. Create the grouped bar plot.
specificity_plot <- ggplot(plot_data, aes(x = gene, y = mean_score, fill = lineage_group)) +
  geom_col(position = position_dodge(width = 0.9), alpha = 0.8) +
  labs(
    title = "Top 5 Most Specifically Essential Genes in Skin Lineage",
    subtitle = "Comparison of mean dependency score in Skin vs. other lineages",
    x = "Gene",
    y = "Mean Dependency Score (CERES)",
    fill = "Lineage Group"
  ) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top") +
  scale_fill_manual(values = c("Skin Lineage" = "#d95f02", "Other Lineages (Mean)" = "#7570b3"))

# 4. Save the plot and print it to the viewer.
ggsave("./Images/top5_skin_specific_genes_barplot.png", plot = specificity_plot, width = 8, height = 6)
message("Specificity bar plot saved to ./Images/top5_skin_specific_genes_barplot.png")
print(specificity_plot)


# --- Part 1g: Statistical Test and Distribution Visualization ---
message("\n--- Visualizing score distributions with statistical significance ---")

# 1. Get the original names of the top 5 specific genes (before cleaning).
top_5_original_names <- head(melanoma_specificity, 5)$gene

# 2. Prepare data for plotting and testing from the full cell line data.
distribution_plot_data <- depmap_essential_annotated %>%
  select(OncotreeLineage, all_of(top_5_original_names)) %>%
  mutate(lineage_group = ifelse(OncotreeLineage == "Skin", "Skin", "Non-Skin")) %>%
  pivot_longer(
    cols = all_of(top_5_original_names),
    names_to = "gene",
    values_to = "dependency_score"
  ) %>%
  mutate(gene = gsub(" \\(.*\\)$", "", gene))

# 3. Ensure the gene order in the plot matches the specificity ranking.
distribution_plot_data$gene <- factor(distribution_plot_data$gene, levels = top_5_for_plot$gene)

# 4. Perform Wilcoxon test for each gene and prepare annotations.
statistical_results <- distribution_plot_data %>%
  group_by(gene) %>%
  summarise(
    p_value = wilcox.test(
      dependency_score[lineage_group == "Skin"],
      dependency_score[lineage_group == "Non-Skin"]
    )$p.value,
    y_pos = max(dependency_score, na.rm = TRUE) + 0.1 # Position for label
  ) %>%
  ungroup() %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    p_label = scales::pvalue(p_adj, accuracy = 0.001, add_p = TRUE)
  )

message("Wilcoxon test results for top 5 skin-specific genes:")
print(statistical_results %>% select(gene, p_value, p_adj))

# 5. Create the violin plot with jittered points and p-value annotations.
distribution_plot <- ggplot(distribution_plot_data, aes(x = lineage_group, y = dependency_score, fill = lineage_group)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_jitter(height = 0, width = 0.2, size = 0.5, alpha = 0.3) +
  geom_text(
    data = statistical_results,
    aes(x = 1.5, y = y_pos, label = p_label),
    inherit.aes = FALSE,
    size = 3.5
  ) +
  facet_wrap(~gene, scales = "free_y", nrow = 1) +
  labs(
    title = "Dependency Score Distribution of Top 5 Skin-Specific Genes",
    subtitle = "Comparison between Skin and Non-Skin cell lineages (adjusted p-values from Wilcoxon test)",
    x = "Lineage Group",
    y = "Dependency Score (CERES)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", strip.text = element_text(face = "bold")) +
  scale_fill_manual(values = c("Skin" = "#d95f02", "Non-Skin" = "#7570b3"))

# 6. Save the plot and print it.
ggsave("./Images/top5_skin_specific_genes_distribution_with_pvals.png", plot = distribution_plot, width = 12, height = 5)
message("Distribution plot with p-values saved to ./Images/top5_skin_specific_genes_distribution_with_pvals.png")
print(distribution_plot)
