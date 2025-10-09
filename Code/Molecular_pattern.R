# A dataset from the pubulication of <<Molecular patterns of resistance to immune checkpoint blockade in melanoma>>
# --- Part 1: Alteration Data and OncoPrint Generation ---

# 1.1 Load and Prepare Alteration Data
print("--- 1.1: Loading and preparing alteration data from Supplementary.xlsx ---")
alteration_data <- read_excel("./Data/molecular_pattern/Supplementary.xlsx", sheet = 2)


# Define column sets for different alteration types
columns_mutation <- c(
"BRAF_600", "BRAF", "NRAS_61", "NRAS", 
"HRAS", "KRAS", "NF1", "TP53", 
"CDKN2A", "PTEN", "DDX3X", "ARID2", 
"PPP6C", "RB1",
"RAC1", "IDH1", "SF3B1", "MAP2K1",
"CTNNB1", "KIT", "CCND1", "CDK4",
"MITF", "TERT", "HLA.A", "HLA.B",
"B2M", "TAP1", "TAP2", "IFNGR1",
"IFNGR2", "JAK1", "JAK2", "STAT1",
"IRF1", "APLNR", "SOCS1", "PTPN2",
"ADAR")

columns_copy <- c(
"BRAF_600_CN", "BRAF_CN", "NRAS_61_CN",
"NRAS_CN", "HRAS_CN", "KRAS_CN", "NF1_CN",
"TP53_CN", "CDKN2A_CN", "PTEN_CN", "DDX3X_CN",
"ARID2_CN", "PPP6C_CN", "RB1_CN", "RAC1_CN",
"IDH1_CN", "SF3B1_CN", "MAP2K1_CN", "CTNNB1_CN",
"KIT_CN", "CCND1_CN", "CDK4_CN", "MITF_CN",
"TERT_CN", "HLA.A_CN", "HLA.B_CN", "B2M_CN",
"TAP1_CN", "TAP2_CN", "IFNGR1_CN", "IFNGR2_CN",
"JAK1_CN", "JAK2_CN", "STAT1_CN", "IRF1_CN",
"APLNR_CN", "SOCS1_CN", "PTPN2_CN", "ADAR_CN")

columns_loh <- c("HLA.A_LOH", "HLA.B_LOH")

# Define the metadata columns we want to use for annotation
columns_meta <- c("treatment", "resistance", "type.of.primary", "prior.CTLA4i", "previous.BRAFi", "gender", "TMB")

# The first column of alteration_data is the sample ID, which we don't need in the data frames anymore
alteration_data <- alteration_data[, -1]

# Slice the main data frame into metadata and specific alteration types
meta_data <- alteration_data[, columns_meta]
mutation_data <- alteration_data[, columns_mutation]
copy_data <- alteration_data[, columns_copy]
loh_data <- alteration_data[, columns_loh]

# Clean the data: Convert string "NA" to actual NA values and ensure TMB is numeric
meta_data <- meta_data %>%
    mutate(across(where(is.character), ~na_if(., "NA"))) %>% # Convert "NA" strings to actual NA
    mutate(TMB = as.numeric(TMB)) # Ensure TMB is numeric

# 1.2 Prepare Data Matrices for OncoPrint
print("--- 1.2: Preparing data matrices for OncoPrint visualization ---")
copy_data_values <- table(unlist(copy_data), useNA = "ifany")
print("Categorical values and their frequencies in copy_data:")
print(copy_data_values)

# Remove the "_CN" suffix from the column names of copy_data
colnames(copy_data) <- sub("_CN$", "", colnames(copy_data))

# Verify that the gene lists are identical between mutation and copy number data
genes_mutation <- colnames(mutation_data)
genes_copy <- colnames(copy_data)

are_genes_identical <- setequal(genes_mutation, genes_copy)

print(paste("Are the gene lists in mutation and copy data identical?", are_genes_identical))

# Create the mutation matrix (MUT)
# A gene is considered mutated if the cell is not empty and not NA.
mut_matrix <- t(mutation_data) # Transpose to have genes as rows
is_mutated <- !is.na(mut_matrix) & mut_matrix != ""
mut_matrix[is_mutated] <- "MUT"
mut_matrix[!is_mutated] <- ""

# Create the copy number variation matrix (CNV)
# We map "Amplification" to "AMP" and "Deletion" to "DEL".
# Based on inspection, "hi amp" is amplification and "deep del" is deletion.
cnv_matrix <- t(copy_data) # Transpose
cnv_matrix[which(cnv_matrix == "hi amp", arr.ind = TRUE)] <- "AMP"
cnv_matrix[which(cnv_matrix == "deep del", arr.ind = TRUE)] <- "DEL"
# Set other non-alteration values to empty strings
cnv_matrix[!cnv_matrix %in% c("AMP", "DEL")] <- ""

# Combine matrices into a single matrix for oncoPrint
# The function expects a single matrix where alterations are semi-colon separated.
combined_matrix <- mut_matrix
# Add CNV data, separating with a semicolon if a mutation is already present
combined_matrix[cnv_matrix != ""] <- paste(combined_matrix[cnv_matrix != ""], cnv_matrix[cnv_matrix != ""], sep = ";")
# Clean up leading semicolons for cells that only have CNV data
combined_matrix <- gsub("^;", "", combined_matrix)

# 1.3 Define Aesthetics and Generate OncoPrint
print("--- 1.3: Defining aesthetics and generating OncoPrint ---")
alter_fun <- list(
    background = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
    },
    # Define colors and functions for each alteration type
    "AMP" = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "red", col = NA))
    },
    "DEL" = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "blue", col = NA))
    },
    "MUT" = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "#008000", col = NA))
    }
)

# Define the colors for the alterations (for barplots and legend)
col <- c(
    "MUT" = "#008000", # green
    "AMP" = "red",
    "DEL" = "blue"
)

# Define colors and create the Heatmap Annotation from metadata
anno_colors <- list(
    treatment = c("anti-CTLA4 resistant" = "#8DD3C7", "anti-PD1 resistant" = "#FFFFB3"),
    resistance = c("primary" = "#E41A1C", "acquired" = "#377EB8", "primary (ipi acq.)" = "#4DAF4A", "NA" = "grey"),
    type.of.primary = c("skin" = "#FDB462", "unknown" = "#B3B3B3", "mucosa" = "#FB8072"),
    prior.CTLA4i = c("yes" = "black", "no" = "white"),
    previous.BRAFi = c("no_" = "lightgrey", "no_During, day 7" = "#BEBADA", "yes_Before" = "#80B1D3"),
    gender = c("F" = "pink", "M" = "lightblue"),
    # For the numeric TMB column, we create a continuous color function
    TMB = colorRamp2(c(0, quantile(meta_data$TMB, 0.95, na.rm = TRUE)), c("white", "darkviolet"))
)

ha_top <- HeatmapAnnotation(
    df = meta_data,
    col = anno_colors,
    na_col = "grey90", # Set a color for NA values in categorical annotations
    annotation_name_gp = gpar(fontsize = 8),
    simple_anno_size = unit(0.5, "cm")
)

# To save the output to a file, we wrap the plotting function in a graphics device.
jpeg("oncoprint.jpg", width = 12, height = 8, units = "in", res = 300)
oncoPrint(combined_matrix, alter_fun = alter_fun, col = col, top_annotation = ha_top)
dev.off() # This closes the file and saves it.

print("--- OncoPrint successfully generated and saved as oncoprint.jpg ---")




# 9. Generate a PCA plot

vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup="groups", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=groups)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  ggtitle("PCA plot of samples") +
  theme_bw()

# Save the PCA plot to a file
ggsave("PCA_plot.jpg", width = 8, height = 6, units = "in", dpi = 300)





# --- Part 2: RNA Data Analysis ---
# This section is for parsing and analyzing the RNA-seq data from the series matrix file.
print("--- Part 2: Parsing RNA-seq metadata (placeholder for full analysis) ---")

# Load the raw series matrix file as lines
rna_data_matrix_lines <- readLines("./Data/molecular_pattern/GSE244982_series_matrix.txt")

# Find the line with sample IDs
sample_id_line_idx <- grep("^!Sample_title", rna_data_matrix_lines)
if (length(sample_id_line_idx) > 0) {
    sample_id_line <- rna_data_matrix_lines[sample_id_line_idx]
    sample_ids_from_matrix <- strsplit(sample_id_line, "\t")[[1]][-1] # Split by tab, remove first element
    sample_ids_from_matrix <- gsub('"', '', sample_ids_from_matrix)

    # Initialize a list to store metadata for each sample
    parsed_meta_list <- vector("list", length(sample_ids_from_matrix))
    names(parsed_meta_list) <- sample_ids_from_matrix

    # Find all lines containing sample characteristics
    char_lines_idx <- grep("^!Sample_characteristics_ch1", rna_data_matrix_lines)

    # Process each characteristics line
    for (line_idx in char_lines_idx) {
        line <- rna_data_matrix_lines[line_idx]
        elements <- strsplit(line, "\t")[[1]]
        char_values <- elements[-1] # Actual values start from the second element

        for (i in seq_along(char_values)) {
            sample_id <- sample_ids_from_matrix[i]
            char_str <- char_values[i]
            char_str_clean <- gsub('"', '', char_str) # Remove quotes
            parts <- strsplit(char_str_clean, ": ", fixed = TRUE)[[1]] # Split by ": "

            if (length(parts) == 2) {
                key <- trimws(parts[1])
                value <- trimws(parts[2])
                parsed_meta_list[[sample_id]][[key]] <- value
            } else {
                warning(paste("Skipping malformed or empty characteristic string:", char_str, "for sample:", sample_id))
            }
        }
    }

    # Convert the list of lists into a data frame
    all_keys <- unique(unlist(lapply(parsed_meta_list, names)))
    rna_meta_df <- data.frame(
        lapply(all_keys, function(key) {
            sapply(parsed_meta_list, function(sample_meta) {
                val <- sample_meta[[key]]
                if (is.null(val)) NA else val
            })
        }),
        stringsAsFactors = FALSE
    )
    colnames(rna_meta_df) <- all_keys
    rownames(rna_meta_df) <- sample_ids_from_matrix
    
    # Clean the rownames of rna_meta_df by removing quotes
    rownames(rna_meta_df) <- gsub('"', '', rownames(rna_meta_df))

    print("Successfully parsed RNA-seq metadata. First 5 rows:")
    print(head(rna_meta_df, 5))

} else {
    warning("Could not find '!Sample_geo_accession' line in GSE244982_series_matrix.txt. Skipping RNA metadata parsing.")
}

# Placeholder for loading the actual RNA expression matrix from the same file.
# The expression data typically starts after a line like '!series_matrix_table_begin'.
# Future analysis would involve parsing this table and aligning it with the metadata above.

print("--- Part 1 complete. Part 2 is ready for further RNA data analysis implementation. ---")

# 2.1 Load RNA-seq Raw Counts
RNAseq_raw <- read.table('./Data/molecular_pattern/GSE244982_ProcessedData_bulkRNAseq.txt', header = TRUE, stringsAsFactors = FALSE, 
sep = '\t', quote = "", check.names = FALSE, row.names = 1)

# 2.2 Filter and Align Metadata with Expression Data
print("--- 2.2: Filtering and aligning RNA-seq data ---")

# Filter out samples with "NA" group from the parsed metadata
rna_meta_filtered <- rna_meta_df[!is.na(rna_meta_df$groups) & rna_meta_df$groups != "NA", ]

# Find common samples between metadata and expression data
common_samples <- intersect(rownames(rna_meta_filtered), colnames(RNAseq_raw))

# Align both metadata and expression data to the common samples
rna_meta_aligned <- rna_meta_filtered[common_samples, ]
rnaseq_aligned <- RNAseq_raw[, common_samples]

# Ensure counts are integers for DESeq2
rnaseq_aligned <- round(rnaseq_aligned)

# Verify alignment
if (!all(rownames(rna_meta_aligned) == colnames(rnaseq_aligned))) {
    stop("Alignment failed between RNA metadata and expression data columns.")
} else {
    print(paste("Successfully aligned", length(common_samples), "samples."))
}

# 2.3 Prepare for Differential Expression Analysis (DEA)
print("--- 2.3: Preparing data for Differential Expression Analysis ---")

# Subset the data to the three main comparison groups
comparison_groups <- c("anti-CTLA4 resistant", "anti-PD1 resistant", "anti-PD1 resistant, BRAFi day7")
meta_for_dea <- rna_meta_aligned[rna_meta_aligned$groups %in% comparison_groups, ]

# Make sure the 'groups' column is a factor for DESeq2
meta_for_dea$groups <- factor(meta_for_dea$groups, levels = comparison_groups)

# Subset the expression matrix to match the DEA metadata
rnaseq_for_dea <- rnaseq_aligned[, rownames(meta_for_dea)]

# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = rnaseq_for_dea,
    colData = meta_for_dea,
    design = ~ groups)

print("DESeq2 object 'dds' created successfully. Ready for differential expression analysis.")
print("The 'mucosal' group is retained in 'rna_meta_aligned' for separate analysis (e.g., as an immune-cold reference).")

# 2.4 Perform Differential Expression Analysis (DEA)
print("--- 2.4: Performing Differential Expression Analysis ---")

# Run DESeq2 without estimating size factors or applying VST/rlog
dds <- estimateSizeFactors(dds, type = "poscounts")  # Use 'poscounts' to avoid errors with log-transformed data
dds <- DESeq(dds, fitType = "parametric", sfType = "poscounts",
    minReplicatesForReplace = Inf) # Specify fitType and sfType

print("DESeq2 analysis completed successfully.")


# Function to extract and print results for a given contrast
get_and_print_results <- function(dds, group1, group2) {
    res <- results(dds, contrast = c("groups", group1, group2))
    resOrdered <- res[order(res$pvalue), ]
    print(paste("Top genes for:", group1, "vs", group2))
    print(head(resOrdered, 10))
    return(res)
}

# Perform pairwise comparisons
group1 <- "anti-CTLA4 resistant"
group2 <- "anti-PD1 resistant"
group3 <- "anti-PD1 resistant, BRAFi day7"

# anti-CTLA4 resistant vs. anti-PD1 resistant
res_CTLA4_vs_PD1 <- get_and_print_results(dds, group1, group2)

# anti-CTLA4 resistant vs. anti-PD1 resistant, BRAFi day7
res_CTLA4_vs_BRAFi <- get_and_print_results(dds, group1, group3)

# anti-PD1 resistant vs. anti-PD1 resistant, BRAFi day7
res_PD1_vs_BRAFi <- get_and_print_results(dds, group2, group3)


# Summary of the results for each group
summary(res_CTLA4_vs_PD1)
summary(res_CTLA4_vs_BRAFi)
summary(res_PD1_vs_BRAFi)

# Count significant genes (p-value < 0.05) for each comparison
sig_CTLA4_vs_PD1 <- sum(res_CTLA4_vs_PD1$pvalue < 0.05, na.rm = TRUE)
sig_CTLA4_vs_BRAFi <- sum(res_CTLA4_vs_BRAFi$pvalue < 0.05, na.rm = TRUE)
sig_PD1_vs_BRAFi <- sum(res_PD1_vs_BRAFi$pvalue < 0.05, na.rm = TRUE)

print(paste("Number of significant genes (p < 0.05) for", group1, "vs", group2, ":", sig_CTLA4_vs_PD1))
print(paste("Number of significant genes (p < 0.05) for", group1, "vs", group3, ":", sig_CTLA4_vs_BRAFi))
print(paste("Number of significant genes (p < 0.05) for", group2, "vs", group3, ":", sig_PD1_vs_BRAFi))

# Order results by adjusted p-value




# Print the first 10 rows of the ordered results
print("First 10 rows of the differential expression results:")
print(head(resOrdered, 10))


# 10. Identify signature genes for each group
print("--- Identifying signature genes for each group ---")

num_signature_genes <- 50 # Number of signature genes to select per group

# Function to get top genes for a given comparison
get_top_genes <- function(res, n = num_signature_genes) {
    res_sig <- res[order(res$pvalue), ]
    head(rownames(res_sig), n)
}

# Get top genes for each comparison
top_CTLA4_vs_PD1 <- get_top_genes(res_CTLA4_vs_PD1)
top_CTLA4_vs_BRAFi <- get_top_genes(res_CTLA4_vs_BRAFi)
top_PD1_vs_BRAFi <- get_top_genes(res_PD1_vs_BRAFi)

# Combine the top genes into a single list and remove duplicates
signature_genes <- unique(c(top_CTLA4_vs_PD1, top_CTLA4_vs_BRAFi, top_PD1_vs_BRAFi))

# 11. Prepare data for heatmap
print("--- Preparing data for heatmap ---")

# Extract the expression data for the signature genes
heatmap_data <- rnaseq_aligned[signature_genes, ]

# Scale the data by row (gene)
heatmap_data <- t(scale(t(heatmap_data)))

# Remove genes/samples with NA values
heatmap_data <- heatmap_data[complete.cases(heatmap_data), complete.cases(t(heatmap_data))]

# 12. Generate Heatmap
print("--- Generating heatmap ---")

# Set color option
hmcols <- colorRampPalette(brewer.pal(11, "RdBu"))(256)

# Order columns by groups
column_order <- order(rna_meta_aligned$groups)



#Define colors for the heatmap annotation
anno_colors <- list(
    groups = c("anti-CTLA4 resistant" = "#8DD3C7",
               "anti-PD1 resistant" = "#FFFFB3",
               "anti-PD1 resistant" = "#FFFFB3", 
               "mucosal" = "#A6CEE3",
               "anti-PD1 resistant, BRAFi day7" = "#BEBADA")
)









# Define annotation for the heatmap
ha <- HeatmapAnnotation(

  groups = rna_meta_aligned$groups,

  col = list(groups = anno_colors$groups)
)

# Create the heatmap
Heatmap(heatmap_data, 
        name = "Z-score", # Title of the heatmap
        col = hmcols, # Colors for the heatmap
        top_annotation = ha, # Add the top annotation
        show_row_names = FALSE, # Hide row names
        show_column_names = FALSE, # Hide column names
        cluster_rows = TRUE, # Cluster rows
        cluster_columns = FALSE, # Do not cluster columns
        column_order = column_order) # Order columns by group

# Save the heatmap to a file
jpeg("heatmap_signature_genes.jpg", width = 10, height = 8, units = "in", res = 300)
Heatmap(heatmap_data, 
       name = "Z-score", # Title of the heatmap
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")), # Colors for the heatmap
        top_annotation = ha, # Add the top annotation
        show_row_names = FALSE, # Hide row names
        show_column_names = FALSE, # Hide column names
        cluster_rows = TRUE, # Cluster rows
        cluster_columns = FALSE, # Do not cluster columns
        column_order = column_order) # Order columns by group
dev.off()

# --- Additional diagnostic checks ---

# 1. Check the distribution of p-values
hist(res$pvalue, main = "Distribution of p-values", xlab = "P-value")

# 2. Check the size factors
print("Size factors:")
print(sizeFactors(dds))

# 3. Check the dispersion estimates

print("Dispersion estimates:")
print(dispersions(dds))

# 4. Check the contrast
print("The contrast used:")
print(resultsNames(dds))



# 13. Perform GO term enrichment analysis
print("--- Performing GO term enrichment analysis ---")

# Install and load necessary packages
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("clusterProfiler")
#BiocManager::install("org.Hs.eg.db")
library(clusterProfiler)
library(org.Hs.eg.db)

# Function to perform GO enrichment analysis
perform_go_enrichment <- function(gene_list, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05) {
  ego <- enrichGO(gene          = gene_list,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = "SYMBOL",
                  ont           = ont,
                  pAdjustMethod = pAdjustMethod,
                  pvalueCutoff  = pvalueCutoff,
                  qvalueCutoff  = qvalueCutoff,
                  readable      = TRUE)
  return(ego)
}

# Perform GO enrichment analysis for each signature
ego_CTLA4_vs_PD1 <- perform_go_enrichment(top_CTLA4_vs_PD1)
ego_CTLA4_vs_BRAFi <- perform_go_enrichment(top_CTLA4_vs_BRAFi)
ego_PD1_vs_BRAFi  <- perform_go_enrichment(top_PD1_vs_BRAFi)

# Print results for each signature
print("GO enrichment results for CTLA4_vs_PD1:")
print(head(summary(ego_CTLA4_vs_PD1), 10))

print("GO enrichment results for CTLA4_vs_BRAFi:")
print(head(summary(ego_CTLA4_vs_BRAFi), 10))

print("GO enrichment results for PD1_vs_BRAFi:")
print(head(summary(ego_PD1_vs_BRAFi), 10))

# You can visualize the results using dotplot
dotplot(ego_CTLA4_vs_PD1, showCategory=10)
dotplot(ego_CTLA4_vs_BRAFi, showCategory=10)
dotplot(ego_PD1_vs_BRAFi, showCategory=10)


# --- Part 3: Single Cell Data Analysis ---

print("--- Starting Part 3: Single Cell Analysis Workflow ---")

# 3.1 Load Single-Cell Data
print("--- 3.1: Loading single-cell data ---")
# Load the raw count matrix. The first column contains gene symbols, so we set it as row names.
sc_counts <- read.table("./Data/molecular_pattern/GSE244983_RawCounts_scRNAseq.txt",
                        header = TRUE,
                        row.names = 1,
                        sep = '\t',
                        check.names = FALSE)

sc_annotation <- read.table("./Data/molecular_pattern/GSE244983_SingleCellAnnotations.txt",
                        header = TRUE,
                        sep = '\t',
                        stringsAsFactors = FALSE)

# 3.1.1 Explore Raw Cell Annotations
print("--- 3.1.1: Exploring the raw sc_annotation file ---")

# Display the structure and first few rows of the annotation data
print("Structure of the sc_annotation data frame:")
str(sc_annotation)
colnames(sc_annotation)[1] <- "CellBarcode" # Rename the first column for clarity

print("First 6 rows of sc_annotation:")
print(head(sc_annotation))

# Summarize each column in the annotation data
print("Summary of each column in sc_annotation:")
for (col_name in colnames(sc_annotation)) {
    column_data <- sc_annotation[[col_name]]
    
    # Treat all columns as categorical unless they have too many unique values (e.g., cell names)
    if (length(unique(column_data)) > 50) {
        print(paste0("--- Skipping '", col_name, "' (high-cardinality with >50 unique values) ---"))
    } else {
        print(paste0("--- Value counts for '", col_name, "' (categorical) ---"))
        print(table(column_data, useNA = "ifany"))
    }
}

# 3.1.2 Prepare and Align Metadata
print("--- 3.1.2: Preparing and aligning cell metadata ---")

# The annotation data should contain the cell barcodes.
if (!"CellBarcode" %in% colnames(sc_annotation)) {
    stop("The 'sc_annotation' data frame must contain a 'NAME' column with cell barcodes.")
}

# Set the row names of the annotation data to be the cell barcodes.
cell_metadata <- sc_annotation %>%
    tibble::column_to_rownames("CellBarcode")

# Find the common cells that exist in both the count matrix and the metadata.
common_cells <- intersect(colnames(sc_counts), rownames(cell_metadata))

if (length(common_cells) == 0) {
    stop("No common cells found between count matrix and annotation file. Check barcode formats.")
}

# Align both the count matrix and the metadata to this common set of cells.
sc_counts_aligned <- sc_counts[, common_cells]
cell_metadata_aligned <- cell_metadata[common_cells, ]

sce <- SingleCellExperiment(
    assays = list(counts = as.matrix(sc_counts_aligned)),
    colData = cell_metadata_aligned
)

print(paste("Loaded single-cell data with", ncol(sce), "cells and", nrow(sce), "genes."))

# 3.2 Normalization
print("--- 3.2: Normalizing single-cell data ---")
# As the data is pre-QC'd, we proceed to normalization.
# We use deconvolution-based size factors for more accurate normalization across cell types.
set.seed(111) # for reproducibility
clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, clusters = clusters)
sce <- logNormCounts(sce)

print("Normalization complete.")

# 3.3 Feature Selection
print("--- 3.3: Selecting highly variable genes (HVGs) ---")
# Model gene variance to identify genes that drive biological heterogeneity.
dec <- modelGeneVar(sce)

# Select the top 2000 most variable genes for downstream analysis.
top_hvgs <- getTopHVGs(dec, n = 2000)

# Visualize the mean-variance trend
png("./Images/sc_mean_variance_plot.png", width = 6, height = 5, units = "in", res = 300)
plot(dec$mean, dec$total, xlab = "Mean log-expression", ylab = "Variance",
     main = "Mean-Variance Trend")
points(dec[top_hvgs, "mean"], dec[top_hvgs, "total"], col = "red")
curve(metadata(dec)$trend(x), col = "dodgerblue", add = TRUE)
dev.off()

print(paste("Identified", length(top_hvgs), "highly variable genes."))

# 3.4 Dimensionality Reduction
print("--- 3.4: Performing dimensionality reduction (PCA, UMAP) ---")
# Perform PCA on the HVGs to capture the main axes of variation.
set.seed(1234)
sce <- runPCA(sce, subset_row = top_hvgs)

# Run UMAP for visualization, using the first 20 principal components.
sce <- runUMAP(sce, dimred = "PCA", n_dimred = 20)

# 3.5 Clustering and Visualization
print("--- 3.5: Clustering cells and visualizing results ---")
# Perform graph-based clustering on the PCA results.
g <- buildSNNGraph(sce, k = 10, use.dimred = "PCA")
clusters_louvain <- igraph::cluster_louvain(g)
colLabels(sce) <- factor(clusters_louvain$membership)

# Visualize clusters on the UMAP plot, colored by the original study's cell types.
p_umap <- plotReducedDim(sce, "UMAP", colour_by = "Sample", text_by = "label",
                         text_size = 4, point_size = 0.5) +
    guides(text = "none") + # Hide legends for clarity
    ggtitle("UMAP of Single Cells by Original Annotation")

ggsave("./Images/melanoma_umap_by_celltype.png", plot = p_umap, width = 10, height = 8)

print("UMAP plot saved to ./Images/melanoma_umap_by_celltype.png")

# 3.6 Find Marker Genes
print("--- 3.6: Finding marker genes for each cell type ---")
# Use a Wilcoxon rank-sum test to find genes upregulated in each cluster vs. others.
markers <- findMarkers(sce, groups = colData(sce)$label, test.type = "wilcox", direction = "up")

# Extract top 10 markers for each cell type for visualization.
top_markers <- lapply(markers, function(x) {
    # Filter for significant markers and take the top 10
    x_df <- as.data.frame(x)
    sig_markers <- rownames(x_df[x_df$FDR < 0.05, ])
    head(sig_markers, 10)
})

# 3.7 Marker Gene Heatmap
print("--- 3.7: Generating marker gene heatmap ---")
# Generate a heatmap of the top marker genes.
p_heatmap <- plotHeatmap(sce,
            features = unique(unlist(top_markers)),
            columns = order(colData(sce)$label),
            colour_columns_by = c("label"),
            cluster_cols = FALSE,
            cluster_rows = TRUE,
            show_colnames = FALSE,
            main = "Top 10 Marker Genes per Cell Type")

# Ensure the output directory exists before saving the plot
if (!dir.exists("./Images")) {
    dir.create("./Images")
}

jpeg("./Images/melanoma_marker_heatmap.jpg", width = 10, height = 12, units = "in", res = 300)
p_heatmap
dev.off()

print("Marker gene heatmap saved to ./Images/melanoma_marker_heatmap.jpg")
print("--- Part 3: Single-cell analysis complete. ---")

# 3.8 Annotate Clusters Based on Marker Genes
print("--- 3.8: Annotating clusters based on marker genes ---")

# Display the top 5 markers for each cluster to aid in manual annotation
print("Top 5 markers for each cluster:")
for (cluster_name in names(top_markers)) {
    print(paste("Cluster", cluster_name, ":", paste(head(top_markers[[cluster_name]], 10), collapse = ", ")))
}

# Based on the markers above and known biology, create a mapping.
# This is a manual step. You will need to replace "Unknown" with the correct cell type.
# Example canonical markers:
# T-cells: CD3D, CD8A, CD4
# B-cells: MS4A1, CD19
# Myeloid: LYZ, CD68, CD14
# Melanoma: MLANA, PMEL
cluster_annotations <- c(
    "1" = "T-cells",
    "2" = "Melanoma",
    "3" = "Myeloid",
    "4" = "B-cells",
    "5" = "T-cells",
    "6" = "Unknown",
    "7" = "Unknown",
    "8" = "Unknown",
    "9" = "Unknown",
    "10" = "Unknown"
    # ... continue for all your clusters
)

# Add the new annotations to the colData of the sce object
sce$manual_annotation <- cluster_annotations[colLabels(sce)]

# Visualize the new manual annotations on the UMAP
p_umap_manual <- plotReducedDim(sce, "UMAP", colour_by = "manual_annotation", text_by = "manual_annotation",
                                text_size = 4, point_size = 0.5) +
    guides(text = "none") +
    ggtitle("UMAP of Single Cells by Manual Annotation")

ggsave("./Images/melanoma_umap_by_manual_annotation.png", plot = p_umap_manual, width = 10, height = 8)

print("UMAP plot with manual annotations saved to ./Images/melanoma_umap_by_manual_annotation.png")
print("Review the top markers and update the 'cluster_annotations' mapping in the script.")

# 3.9 Automated Cluster Annotation with SingleR
print("--- 3.9: Performing automated cluster annotation with SingleR ---")

# Load the Human Primary Cell Atlas data as a reference
# This might take a moment to download the first time you run it.
hpca.ref <- HumanPrimaryCellAtlasData()

# Perform the annotation. SingleR will compare each cell in your dataset
# to the reference and assign the most likely cell type.
predictions <- SingleR(test = sce, ref = hpca.ref, labels = hpca.ref$label.main)

# Add the predicted labels to your SingleCellExperiment object
sce$singler_annotation <- predictions$labels

# To see how the automated labels correspond to your unsupervised clusters,
# you can create a summary table.
print("Cross-tabulation of automated labels vs. unsupervised clusters:")
print(table(Predicted = sce$singler_annotation, Cluster = colLabels(sce)))

# Visualize the new automated annotations on the UMAP plot
p_umap_singler <- plotReducedDim(sce, "UMAP", colour_by = "singler_annotation", text_by = "singler_annotation",
                                 text_size = 4, point_size = 0.5) +
    guides(text = "none") +
    ggtitle("UMAP of Single Cells by SingleR Automated Annotation")

ggsave("./Images/melanoma_umap_by_singler_annotation.png", plot = p_umap_singler, width = 10, height = 8)

print("UMAP plot with SingleR annotations saved to ./Images/melanoma_umap_by_singler_annotation.png")

# You can also view the diagnostic plot to see the annotation scores
png("./Images/singler_diagnostic_plot.png", width = 10, height = 6, units = "in", res = 300)
plotScoreHeatmap(predictions)
dev.off()

print("SingleR diagnostic heatmap saved to ./Images/singler_diagnostic_plot.png")
