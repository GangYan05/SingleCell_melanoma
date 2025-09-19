# --- Script for Generating an OncoPrint from MSK-IMPACT Melanoma Data ---

# 1. Load Required Libraries
library(ComplexHeatmap)
library(grid)
library(dplyr)
library(tibble)
library(tidyr)
library(circlize)

# --- Step 1: Define Genes of Interest and Load Data ---

# Set this flag to TRUE to use the 'genes_of_interest' list.
# Set to FALSE to automatically find the top N most frequently altered genes.
use_specific_genes <- FALSE
genes_of_interest <- c("BRAF", "NRAS", "KIT", "CDKN2A", "PTEN", "TP53")

# Load genomic alteration data
mutation_data <- read.table('./Data/mel_mskimpact_2020/data_mutations.txt', header = TRUE, stringsAsFactors = FALSE, sep = '\t', quote = "")
cna_data_raw <- read.table('./Data/mel_mskimpact_2020/data_CNA.txt', header = TRUE, stringsAsFactors = FALSE, sep = '\t', check.names = FALSE)

cna_data <- cna_data_raw %>%
  # Remove non-sample ID columns before aggregation to prevent errors.
  select(-any_of(c("Entrez_Gene_Id", "Cytoband"))) %>%
  group_by(Hugo_Symbol) %>%
  # For each sample column, find the value from the row with the maximum absolute value.
  summarise(across(everything(), ~ .[which.max(abs(.))]), .groups = "drop") %>%
  # Now that Hugo_Symbol is unique, set it as rownames.
  column_to_rownames("Hugo_Symbol")

message("Loaded and processed CNA data, handling duplicate gene symbols.")

# --- Step 2: Process Raw Alteration Data ---

# A. Process Mutation Data
# Convert from long format to a wide matrix (genes x samples)
mut_matrix <- mutation_data %>%
  # A more refined approach: exclude silent mutations which are often non-functional.
  filter(Variant_Classification != "Silent") %>%
  # We only want one entry per gene-sample pair, with "MUT"
  select(Hugo_Symbol, Tumor_Sample_Barcode) %>% 
  distinct() %>%
  mutate(alteration = "MUT") %>%
  pivot_wider(names_from = Tumor_Sample_Barcode, values_from = alteration, values_fill = "")

# Convert to a standard matrix with genes as rownames
mut_matrix <- as.data.frame(mut_matrix)
rownames(mut_matrix) <- mut_matrix$Hugo_Symbol
mut_matrix$Hugo_Symbol <- NULL

# B. Process CNA Data
# Recode numeric values to alteration strings ("AMP", "HOMDEL")
cna_matrix <- cna_data %>%
  mutate(across(everything(), ~case_when(
    . == 2  ~ "AMP",    # High-level amplification
    . == 1  ~ "",       # Low-level gain (removed)
    . == -2 ~ "HOMDEL", # Deep deletion
    . == -1 ~ "",       # Shallow deletion (removed)
    TRUE    ~ ""
  )))

# --- Step 3: Select Genes for Plotting ---

if (use_specific_genes) {
  # Use the predefined list of genes.
  genes_to_plot <- genes_of_interest
  message(sprintf("Using %d user-defined genes of interest for inquiry.", length(genes_to_plot)))
  plot_title <- "Genomic Alterations in Genes of Interest"
  
} else {
  # OPTIMIZATION: Calculate frequencies first to find top altered genes
  # before creating the final combined matrix.
  message("Calculating alteration frequencies to find top altered genes...")
  
  # Find common genes and samples between the two matrices
  common_genes_freq <- intersect(rownames(mut_matrix), rownames(cna_matrix))
  common_samples_freq <- intersect(colnames(mut_matrix), colnames(cna_matrix))
  
  # Create boolean matrices indicating any alteration
  mut_altered <- mut_matrix[common_genes_freq, common_samples_freq] != ""
  cna_altered <- cna_matrix[common_genes_freq, common_samples_freq] != ""
  
  # Combine to see if a gene is altered by MUT or CNA in a sample
  any_alteration <- mut_altered | cna_altered
  
  # Calculate the frequency for each gene
  alteration_frequency <- rowSums(any_alteration) / ncol(any_alteration)
  
  # Sort genes by frequency in descending order
  sorted_genes <- names(sort(alteration_frequency, decreasing = TRUE))
  
  # Select the top N genes for the plot (e.g., 50)
  num_genes_to_plot <- 50
  genes_to_plot <- head(sorted_genes, num_genes_to_plot)
  
  message(sprintf("Selected the top %d most frequently altered genes for the OncoPrint.", num_genes_to_plot))
  plot_title <- sprintf("Genomic Landscape of Top %d Altered Genes in MSK-IMPACT Melanoma", num_genes_to_plot)
}

# --- Step 4: Create Final OncoPrint Matrix ---
main_plot_genes <- Reduce(intersect, list(genes_to_plot,
                                          rownames(mut_matrix),
                                          rownames(cna_matrix)))

final_common_samples <- Reduce(intersect, list(colnames(mut_matrix),
                                               colnames(cna_matrix)))

if(length(main_plot_genes) == 0){
  stop("None of the selected genes were found in the alteration data. Please check gene names.")
}

message(sprintf("Generating OncoPrint for %d genes and %d samples.",
                length(main_plot_genes), length(final_common_samples)))

# A. Create the final OncoPrint matrix
mut_matrix_subset <- mut_matrix[main_plot_genes, final_common_samples, drop = FALSE]
cna_matrix_subset <- cna_matrix[main_plot_genes, final_common_samples, drop = FALSE]

oncoprint_matrix_subset <- matrix("", nrow = nrow(mut_matrix_subset), ncol = ncol(mut_matrix_subset),
                                  dimnames = dimnames(mut_matrix_subset))

for (gene in rownames(oncoprint_matrix_subset)) {
  for (sample in colnames(oncoprint_matrix_subset)) {
    mut_val <- mut_matrix_subset[gene, sample]
    cna_val <- cna_matrix_subset[gene, sample]
    combined_val <- paste(c(mut_val, cna_val)[c(mut_val, cna_val) != ""], collapse = ";")
    oncoprint_matrix_subset[gene, sample] <- combined_val
  }
}

message("Final OncoPrint matrix created.")

# --- Step 5: Construct and Draw OncoPrint ---

# A. Define plotting parameters (colors and drawing functions)
col <- c("MUT" = "#008000", "AMP" = "red", "HOMDEL" = "blue")
alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
  },
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"), gp = gpar(fill = col["AMP"], col = NA))
  },
  HOMDEL = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"), gp = gpar(fill = col["HOMDEL"], col = NA))
  },
  MUT = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h * 0.33, gp = gpar(fill = col["MUT"], col = NA))
  }
)

# B. Create the OncoPrint object
# Note: We remove the column_title here as it will be added to the combined plot.
oncoprint_plot <- oncoPrint(oncoprint_matrix_subset,
          alter_fun = alter_fun, 
          # This is needed to resolve sample alignment issues with complex annotations
          alter_fun_is_vectorized = FALSE,
          col = col,
          # Control font sizes for row names and the frequency percentage text
          row_names_gp = gpar(fontsize = 8),
          pct_gp = gpar(fontsize = 8),
          # Keep all columns, even if empty, to match the annotation dimensions
          remove_empty_columns = FALSE,
          # This function tells oncoPrint how to parse strings with multiple alterations (e.g., "MUT;AMP")
          get_type = function(x) strsplit(x, ";")[[1]],
          heatmap_legend_param = list(title = "Alterations",
                                      at = c("AMP", "HOMDEL", "MUT"),
                                      labels = c("Amplification", "Deep Deletion", "Mutation")))

# C. Draw the plot to a file
png("./Images/msk_melanoma_oncoprint.png", width = 12, height = 8, units = "in", res = 300)
draw(oncoprint_plot, column_title = plot_title)
dev.off()

message("OncoPrint saved to ./Images/msk_melanoma_oncoprint.png")

# To display in an interactive R session, just run the oncoPrint command without png() and dev.off()
