# --- Script for Generating an OncoPrint from MSK-IMPACT Melanoma Data ---

# --- Step 1: Define Genes of Interest and Load Data ---

# Set this flag to TRUE to use the 'genes_of_interest' list.
# Set to FALSE to automatically find the top N most frequently altered genes.
use_specific_genes <- TRUE
genes_of_interest <- unique(
  c(
    "BRAF", "RET", "NTRK1", "NTRK2", "NTRK3",
    "RET", "KIT", "MTAP", "MAP2K1", "NRG1", "NRAS",
    "ERBB2", "TP53", "MDM2", "MTOR", "CCNE1", "CDKN2A",
    "MDM2", "KRAS", "NF1", "FGFR1", "FGFR2", "MET",
    "PIK3CA", "FBXW7", "CDK12", "ARID1A", "PPP2R1A", "FGFR3", "PTEN"
  )
)

# load clinical data
clinical_meta <- read.table("./Data/mel_mskimpact_2020/data_clinical_patient.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t", quote = "")
sample_meta <- read.table("./Data/mel_mskimpact_2020/data_clinical_sample.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t", quote = "")

# clinical data for annotations will be aligned in Step 5

# Load genomic alteration data
mutation_data <- read.table("./Data/mel_mskimpact_2020/data_mutations.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t", quote = "")
cna_data_raw <- read.table("./Data/mel_mskimpact_2020/data_CNA.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)

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
  mutate(across(everything(), ~ case_when(
    . == 2 ~ "AMP", # High-level amplification
    . == 1 ~ "", # Low-level gain (removed)
    . == -2 ~ "HOMDEL", # Deep deletion
    . == -1 ~ "", # Shallow deletion (removed)
    TRUE ~ ""
  )))

# --- Step 3: Select Genes for Plotting ---

if (use_specific_genes) {
  # Use the predefined list of genes.
  genes_to_plot <- genes_of_interest
  message(sprintf("Using %d user-defined genes of interest for inquiry.", length(genes_to_plot)))
  plot_title <- "MSK-IMPACT Melanoma"
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
main_plot_genes <- Reduce(intersect, list(
  genes_to_plot,
  rownames(mut_matrix),
  rownames(cna_matrix)
))

final_common_samples <- Reduce(intersect, list(
  colnames(mut_matrix),
  colnames(cna_matrix)
))

if (length(main_plot_genes) == 0) {
  stop("None of the selected genes were found in the alteration data. Please check gene names.")
}

message(sprintf(
  "Generating OncoPrint for %d genes and %d samples.",
  length(main_plot_genes), length(final_common_samples)
))

# A. Create the final OncoPrint matrix
mut_matrix_subset <- mut_matrix[main_plot_genes, final_common_samples, drop = FALSE]
cna_matrix_subset <- cna_matrix[main_plot_genes, final_common_samples, drop = FALSE]

oncoprint_matrix_subset <- matrix("",
  nrow = nrow(mut_matrix_subset), ncol = ncol(mut_matrix_subset),
  dimnames = dimnames(mut_matrix_subset)
)

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

# A. Prepare Column Annotations (TMB and OS)
# Align clinical data with the samples in the oncoprint
sample_ids <- final_common_samples
sample_idx <- match(sample_ids, sample_meta$SAMPLE_ID)
tmb_aligned <- sample_meta$TMB_NONSYNONYMOUS[sample_idx]

# Map samples to patients to get survival data
patient_ids <- sample_meta$PATIENT_ID[sample_idx]
patient_idx <- match(patient_ids, clinical_meta$PATIENT_ID)
os_aligned <- clinical_meta$OS_MONTHS[patient_idx]

# Create the annotation object
bottom_anno <- HeatmapAnnotation(
  TMB = anno_points(tmb_aligned, height = unit(1.5, "cm"), 
                     gp = gpar(fill = "grey"),
                     axis_param = list(side = "left", gp = gpar(fontsize = 8))),
  OS_Months = anno_points(os_aligned, height = unit(1.5, "cm"),
                        axis_param = list(side = "left", gp = gpar(fontsize = 8))),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 10)
)

# B. Define plotting parameters (colors and drawing functions)
col <- c("MUT" = "#008000", "AMP" = "red", "HOMDEL" = "blue")
alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
  },
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w - unit(1, "mm"), h - unit(1, "mm"), gp = gpar(fill = col["AMP"], col = NA))
  },
  HOMDEL = function(x, y, w, h) {
    grid.rect(x, y, w - unit(1, "mm"), h - unit(1, "mm"), gp = gpar(fill = col["HOMDEL"], col = NA))
  },
  MUT = function(x, y, w, h) {
    grid.rect(x, y, w - unit(1, "mm"), h * 0.66, gp = gpar(fill = col["MUT"], col = NA))
  }
)

# B. Create the OncoPrint object
# Note: We remove the column_title here as it will be added to the combined plot.
oncoprint_plot <- oncoPrint(oncoprint_matrix_subset,
  alter_fun = alter_fun,
  # Add the aligned clinical annotations to the bottom
  bottom_annotation = bottom_anno,
  # This is needed to resolve sample alignment issues with complex annotations
  alter_fun_is_vectorized = FALSE,
  col = col,
  # Control font sizes for row names and the frequency percentage text
  row_names_gp = gpar(fontsize = 10),
  pct_gp = gpar(fontsize = 10),
  # Keep all columns, even if empty, to match the annotation dimensions
  remove_empty_columns = FALSE,
  # This function tells oncoPrint how to parse strings with multiple alterations (e.g., "MUT;AMP")
  get_type = function(x) strsplit(x, ";")[[1]],
  heatmap_legend_param = list(
    title = "Alterations",
    at = c("AMP", "HOMDEL", "MUT"),
    labels = c("Amplification", "Deep Deletion", "Mutation")
  )
)

# C. Draw the plot to a file
png("./Images/MSK/driver_mutations_msk_melanoma_oncoprint.png", width = 12, height = 8, units = "in", res = 300)
draw(oncoprint_plot, column_title = plot_title)
dev.off()

message("OncoPrint saved to ./Images/MSK/msk_melanoma_oncoprint.png")

# --- Step 6: Calculate Mutation Co-occurrence ---

# Function ported from TCGA_melanoma.R
run_co_occurrence_module <- function(matrices, output_csv, output_plot) {
  message("Running Mutation Co-occurrence Analysis...")
  
  # 1. Binary Matrix (MUT only)
  genes <- rownames(matrices$mut)
  samples <- colnames(matrices$mut)
  mut_bin <- matrix(0, nrow = length(genes), ncol = length(samples), dimnames = list(genes, samples))
  
  for (g in genes) {
    for (s in samples) {
      if (matrices$mut[g, s] == "MUT") mut_bin[g, s] <- 1
    }
  }
  
  # 2. Fisher's Test
  res <- data.frame()
  n_genes <- length(genes)
  if (n_genes < 2) return(NULL)
  
  for (i in 1:(n_genes - 1)) {
    for (j in (i + 1):n_genes) {
      g1 <- genes[i]; g2 <- genes[j]
      tbl <- table(factor(mut_bin[g1, ], levels = c(1, 0)), factor(mut_bin[g2, ], levels = c(1, 0)))
      f <- fisher.test(tbl)
      res <- rbind(res, data.frame(Gene1 = g1, Gene2 = g2, Odds_Ratio = f$estimate, P_Value = f$p.value, CI_Lower = f$conf.int[1], CI_Upper = f$conf.int[2]))
    }
  }
  res$FDR <- p.adjust(res$P_Value, method = "fdr")
  write.csv(res, output_csv, row.names = FALSE)
  
  # 3. Filter for Plotting (Significant only)
  sig_res <- res %>% filter(FDR < 0.05)
  plot_genes <- sort(unique(c(sig_res$Gene1, sig_res$Gene2)))
  
  if (length(plot_genes) < 2) {
    message("Not enough significant interactions for heatmap.")
    return()
  }
  
  # 4. Heatmap Matrix
  mat_or <- matrix(0, length(plot_genes), length(plot_genes), dimnames = list(plot_genes, plot_genes))
  mat_sig <- matrix("", length(plot_genes), length(plot_genes), dimnames = list(plot_genes, plot_genes))
  
  for (i in 1:nrow(res)) {
    g1 <- res$Gene1[i]; g2 <- res$Gene2[i]
    if (g1 %in% plot_genes & g2 %in% plot_genes) {
      or <- log2(res$Odds_Ratio[i])
      if (is.infinite(or)) or <- ifelse(or > 0, 5, -5)
      mat_or[g1, g2] <- or; mat_or[g2, g1] <- or
      
      if (res$FDR[i] < 0.05) {
        # Use simple asterisk for wider compatibility
        mat_sig[g1, g2] <- "*"
        mat_sig[g2, g1] <- "*"
        
        # Add more stars for higher significance
        if (res$FDR[i] < 0.01) {
            mat_sig[g1, g2] <- "**"
            mat_sig[g2, g1] <- "**"
        }
        if (res$FDR[i] < 0.001) {
            mat_sig[g1, g2] <- "***"
            mat_sig[g2, g1] <- "***"
        }
      }
    }
  }
  
  if (!requireNamespace("circlize", quietly = TRUE)) library(circlize)
  col_fun <- colorRamp2(c(-3, 0, 3), c("purple", "white", "green"))
  ht <- Heatmap(mat_or, name = "Log2 OR", col = col_fun,
                cell_fun = function(j, i, x, y, width, height, fill) {
                  # Check if significance label exists for this cell
                  # Using numeric indices i, j on the matrix
                  # Note: ComplexHeatmap cell_fun receives i (row index) and j (col index) of the ORIGINAL matrix
                  sig_char <- mat_sig[i, j]
                  if (sig_char != "") {
                    grid.text(sig_char, x, y, gp = gpar(fontsize = 20, col = "black"))
                  }
                },
                cluster_rows = TRUE, cluster_columns = TRUE, 
                show_row_dend = FALSE, show_column_dend = FALSE,
                heatmap_legend_param = list(title = "Interaction", at = c(-3, 0, 3), labels = c("Mutual Exclusivity", "Random", "Co-occurrence"), legend_direction = "horizontal"))
  
  dir.create(dirname(output_plot), showWarnings = FALSE, recursive = TRUE)
  png(output_plot, width = 10, height = 11, units = "in", res = 300)
  draw(ht, heatmap_legend_side = "bottom", column_title = "MSK-IMPACT Mutation Co-occurrence")
  dev.off()
  message(paste("Co-occurrence heatmap saved to", output_plot))
}

# Execute the module using MSK data
# We construct a list 'matrices' to match the function's expected input format
matrices_msk <- list(
  mut = mut_matrix_subset
)

run_co_occurrence_module(
  matrices = matrices_msk,
  output_csv = "./Data/Results/MSK/mutation_co_occurrence_msk_melanoma.csv",
  output_plot = "./Images/MSK/mutation_co_occurrence_heatmap_msk_melanoma.png"
)

message("Analysis Complete.")
