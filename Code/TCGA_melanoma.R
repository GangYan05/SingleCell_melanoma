# --- Script for Generating an OncoPrint from TCGA Melanoma Data ---

# --- Step 1: Define Genes of Interest and Load Data ---

# Set this flag to TRUE to use the 'genes_of_interest' list.
# Set to FALSE to automatically find the top N most frequently altered genes.
use_specific_genes <- FALSE
genes_of_interest <- c("BRAF", "NRAS", "KIT", "CDKN2A", "PTEN", "TP53")

# Load genomic alteration data
mutation_data <- read.table('./Data/skcm_tcga_pan_can_atlas_2018/data_mutations.txt', header = TRUE, stringsAsFactors = FALSE, sep = '\t', quote = "")
cna_data_raw <- read.table('./Data/skcm_tcga_pan_can_atlas_2018/data_CNA.txt', header = TRUE, stringsAsFactors = FALSE, sep = '\t', check.names = FALSE)

cna_data <- cna_data_raw %>%
  # Remove non-sample ID columns before aggregation to prevent errors.
  select(-any_of(c("Entrez_Gene_Id", "Cytoband"))) %>%
  group_by(Hugo_Symbol) %>%
  # For each sample column, find the value from the row with the maximum absolute value.
  summarise(across(everything(), ~ .[which.max(abs(.))]), .groups = "drop") %>%
  # Now that Hugo_Symbol is unique, set it as rownames.
  column_to_rownames("Hugo_Symbol")

message("Loaded and processed CNA data, handling duplicate gene symbols.")

# --- Step 2: Load and Process RNA Expression Data ---

# Load RNA-seq data (log2-transformed RSEM values)
rna_data_raw <- read.table('./Data/skcm_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem.txt', header = TRUE, stringsAsFactors = FALSE, sep = '\t', check.names = FALSE)

rna_data <- rna_data_raw %>%
  select(-any_of("Entrez_Gene_Id")) %>%
  # For duplicate genes, take the mean expression value
  group_by(Hugo_Symbol) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  # Now that Hugo_Symbol is unique, set it as rownames.
  column_to_rownames("Hugo_Symbol")

message("Loaded and processed RNA-seq data.")

# --- Step 3: Load and Process Protein Expression Data ---

# Load RPPA data
protein_data_raw <- read.table('./Data/skcm_tcga_pan_can_atlas_2018/data_rppa.txt', header = TRUE, stringsAsFactors = FALSE, sep = '\t', check.names = FALSE)

# The first column contains protein/antibody IDs. We'll rename it for clarity
# and create a new 'Hugo_Symbol' column with corrections.
colnames(protein_data_raw)[1] <- "Protein_ID"

protein_data <- protein_data_raw %>%
  # Create a 'Hugo_Symbol' column by first extracting the gene symbol(s) before the '|'.
  mutate(Hugo_Symbol = gsub("\\|.*", "", Protein_ID)) %>%
  
  # For IDs with multiple genes (e.g., "AKT1 AKT2 AKT3"), create a separate row for each gene.
  # The expression values will be duplicated for each gene in the family.
  separate_rows(Hugo_Symbol, sep = " ") %>%
  
  # Then, on the result, strip phospho-site and other common suffixes.
  mutate(Hugo_Symbol = gsub("_p[sStTyY].*|_phospho|_.*", "", Hugo_Symbol)) %>%
  
  # Manually correct common inconsistencies where the base protein name
  # is not the official HUGO gene symbol.
  mutate(Hugo_Symbol = case_when(
    Hugo_Symbol == "4E-BP1" ~ "EIF4EBP1",
    Hugo_Symbol == "ACC" ~ "ACACA",
    Hugo_Symbol == "AMPKalpha" ~ "PRKAA1",
    Hugo_Symbol == "c-Kit" ~ "KIT",
    Hugo_Symbol == "c-Met" ~ "MET",
    Hugo_Symbol == "eIF4G" ~ "EIF4G1",
    Hugo_Symbol == "GSK-3alpha-beta" ~ "GSK3B", # Map to beta as it's more common
    Hugo_Symbol == "IGF-IR" ~ "IGF1R",
    Hugo_Symbol == "Lck" ~ "LCK",
    Hugo_Symbol == "MEK1" ~ "MAP2K1",
    Hugo_Symbol == "mTOR" ~ "MTOR",
    Hugo_Symbol == "p70S6K" ~ "RPS6KB1",
    Hugo_Symbol == "PDK1" ~ "PDPK1",
    Hugo_Symbol == "S6" ~ "RPS6",
    TRUE ~ Hugo_Symbol
  )) %>%
  
  # Now, group by the corrected Hugo_Symbol and average the values.
  # This handles multiple antibodies/phosphosites mapping to the same gene.
  group_by(Hugo_Symbol) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  
  # Now that Hugo_Symbol is unique, set it as rownames.
  column_to_rownames("Hugo_Symbol")

message("Loaded and processed RPPA protein data with improved gene symbol mapping.")

# --- Step 4: Process Raw Alteration Data ---

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

# --- Step 5: Select Genes for Plotting ---

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
  plot_title <- sprintf("Genomic Landscape of Top %d Altered Genes in TCGA Melanoma", num_genes_to_plot)
}

# --- Step 6: Create Final Matrices for Plotting ---

# The gene list for the OncoPrint and RNA heatmap will be determined by the
# intersection of alteration and RNA data, NOT filtered by protein data availability.
main_plot_genes <- Reduce(intersect, list(genes_to_plot,
                                          rownames(mut_matrix),
                                          rownames(cna_matrix),
                                          rownames(rna_data)))

# The samples for the main plot will be common across the alteration and RNA datasets.
final_common_samples <- Reduce(intersect, list(colnames(mut_matrix),
                                               colnames(cna_matrix),
                                               colnames(rna_data)))

if(length(main_plot_genes) == 0){
  stop("None of the selected genes were found across alteration and RNA datasets.")
}

message(sprintf("Generating OncoPrint/RNA plot for %d genes and %d samples.",
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

# B. Create the final Expression matrix
expression_subset <- rna_data[main_plot_genes, final_common_samples, drop = FALSE]

# Z-score the expression data by gene (row) for visualization.
# The t() is used because scale() works on columns.
expression_z_score <- t(scale(t(expression_subset)))
# Replace any potential NaN values (from genes with no variance) with 0.
expression_z_score[is.na(expression_z_score)] <- 0

message("Final OncoPrint and Expression matrices created.")

# --- Step 7: Construct and Combine Plots ---

# A. Define plotting parameters (colors and drawing functions)
col <- c("MUT" = "#008000", "AMP" = "red", "HOMDEL" = "blue")
expr_col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
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

# C. Create the Expression Heatmap object
expression_heatmap <- Heatmap(expression_z_score,
                              name = "Expression",
                              col = expr_col_fun,
                              show_row_names = FALSE, # Genes are already labeled in the oncoprint
                              cluster_rows = FALSE, # Keep the gene order from the oncoprint
                              show_column_names = FALSE,
                              cluster_columns = FALSE, # Use the column order from the oncoprint
                              height = unit(4, "cm"),
                              heatmap_legend_param = list(title = "Expr. (Z-score)"))

# D. Combine the plots vertically and draw to a file
png("./Images/melanoma_oncoprint_with_expression.png", width = 12, height = 8, units = "in", res = 300)
combined_plot_list <- oncoprint_plot %v% expression_heatmap
draw(combined_plot_list, 
     column_title = plot_title, 
     heatmap_legend_side = "right", 
     annotation_legend_side = "right")
dev.off()

message("OncoPrint with expression data saved to ./Images/melanoma_oncoprint_with_expression.png")

# --- Step 8: Independent Protein Data Analysis ---
message("\n--- Generating independent heatmap for protein data ---")

# Z-score the full protein data matrix for visualization
protein_full_z_score <- t(scale(t(protein_data)))
protein_full_z_score[is.na(protein_full_z_score)] <- 0

# Create the heatmap object
protein_full_heatmap <- Heatmap(protein_full_z_score,
                                name = "Protein (Z-score)",
                                col = expr_col_fun,
                                show_column_names = FALSE,
                                row_names_gp = gpar(fontsize = 4),
                                column_title = "Protein Expression Landscape (RPPA)",
                                cluster_rows = TRUE,
                                cluster_columns = TRUE)

# Draw the heatmap to a separate file
png("./Images/protein_expression_heatmap.png", width = 10, height = 8, units = "in", res = 300)
draw(protein_full_heatmap)
dev.off()

message("Independent protein expression heatmap saved to ./Images/protein_expression_heatmap.png")

# To display in an interactive R session, just run the oncoPrint command without png() and dev.off()
