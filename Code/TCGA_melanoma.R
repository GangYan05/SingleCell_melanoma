# --- Script for Generating an OncoPrint from TCGA Melanoma Data ---

# 1. Load Required Libraries
library(ComplexHeatmap)
library(grid)
library(dplyr)
library(tibble)
# --- Step 1: Load Gene List and Raw Alteration Data ---

# Load the list of top melanoma-specific genes you generated in 'Gene Dependency.R'
top_genes <- read.table("./Results/melanoma_specific_genes.txt", header = FALSE, col.names = "gene")$gene
message(sprintf("Loaded %d top melanoma-specific genes.", length(top_genes)))

# Load mutation data (long format)
mutation_data <- read.table('./Data/skcm_tcga_pan_can_atlas_2018/data_mutations.txt', header = TRUE, stringsAsFactors = FALSE, sep = '\t', quote = "")

# Load Copy Number Alteration (CNA) data (wide format)
# The CNA file can have duplicate gene symbols (Hugo_Symbol). We must read it in
# and aggregate the duplicates before processing. A common strategy is to keep
# the alteration with the largest absolute value (e.g., prefer Amp/Del over Gain/Loss).
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

# --- Step 2: Process and Combine Data into OncoPrint Matrix ---

# A. Process Mutation Data
# Convert from long format to a wide matrix (genes x samples)
mut_matrix <- mutation_data %>%
  filter(Hugo_Symbol %in% top_genes) %>%
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
  filter(rownames(.) %in% top_genes) %>%
  mutate(across(everything(), ~case_when(
    . == 2  ~ "AMP",
    . == -2 ~ "HOMDEL",
    TRUE    ~ ""
  )))

# C. Merge Mutation and CNA data
# Ensure both matrices have the same set of genes and samples
common_genes <- intersect(rownames(mut_matrix), rownames(cna_matrix))
common_samples <- intersect(colnames(mut_matrix), colnames(cna_matrix))

mut_matrix <- mut_matrix[common_genes, common_samples]
cna_matrix <- cna_matrix[common_genes, common_samples]

# Create the final matrix by combining the two.
# If a gene/sample has both a MUT and a CNA, they will be pasted together (e.g., "MUT;AMP")
oncoprint_matrix <- matrix("", nrow = length(common_genes), ncol = length(common_samples),
                           dimnames = list(common_genes, common_samples))

for (gene in common_genes) {
  for (sample in common_samples) {
    mut_val <- mut_matrix[gene, sample]
    cna_val <- cna_matrix[gene, sample]
    
    # Combine with a semicolon if both are present
    combined_val <- paste(c(mut_val, cna_val)[c(mut_val, cna_val) != ""], collapse = ";")
    oncoprint_matrix[gene, sample] <- combined_val
  }
}

message("Processed and combined mutation and CNA data into the final oncoprint matrix.")

# 3. Define Colors and Drawing Functions for Alterations
# This part is crucial for customizing the look of your OncoPrint.

# A named vector mapping alteration types to colors.
col <- c("MUT" = "#008000", "AMP" = "red", "HOMDEL" = "blue")

# A list of functions that tell oncoPrint how to draw each alteration type.
alter_fun <- list(
  # The background for cells with no alteration
  background = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # How to draw a HOMDEL (a solid blue rectangle)
  HOMDEL = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"), gp = gpar(fill = col["HOMDEL"], col = NA))
  },
  # How to draw an AMP (a solid red rectangle)
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"), gp = gpar(fill = col["AMP"], col = NA))
  },
  # How to draw a MUT (a small green rectangle, taking up 1/3 of the cell height)
  MUT = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h * 0.33, gp = gpar(fill = col["MUT"], col = NA))
  }
)

png("./Images/melanoma_oncoprint.png", width = 10, height = 6, units = "in", res = 300)
oncoPrint(oncoprint_matrix,
          alter_fun = alter_fun, 
          col = col,
          # This function tells oncoPrint how to parse strings with multiple alterations (e.g., "MUT;AMP")
          get_type = function(x) strsplit(x, ";")[[1]],
          column_title = "Genomic Alterations in Top Melanoma-Specific Dependency Genes",
          heatmap_legend_param = list(title = "Alterations", at = c("AMP", "HOMDEL", "MUT"),
                                      labels = c("Amplification", "Deep Deletion", "Mutation")))
dev.off()

message("OncoPrint saved to ./Images/melanoma_oncoprint.png")

# To display in an interactive R session, just run the oncoPrint command without png() and dev.off()
