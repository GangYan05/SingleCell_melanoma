# A dataset from the pubulication of <<Molecular patterns of resistance to immune checkpoint blockade in melanoma>>
# --- Part 1: Alteration Data and OncoPrint Generation ---

# 1.1 Load and Prepare Alteration Data
print("--- 1.1: Loading and preparing alteration data from Alteration_data.xlsx ---")
alteration_data <- read_excel("./Data/molecular_pattern/Alteration_data.xlsx")

genes_of_interest <- unique(c(
  "BRAF", "RET", "NTRK1", "NTRK2", "NTRK3",
  "RET", "KIT", "MTAP", "MAP2K1", "NRG1", "NRAS",
  "ERBB2", "TP53", "MDM2", "MTOR", "CCNE1", "CDKN2A",
  "MDM2", "KRAS", "NF1", "FGFR1", "FGFR2", "MET",
  "PIK3CA", "FBXW7", "CDK12", "ARID1A", "PPP2R1A", "FGFR3", "PTEN"
))


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

# Add missing genes of interest to mutation and copy number column sets
# To ensure we have everything needed for the GOI plot and the original plot
all_mut_genes <- unique(c(columns_mutation, genes_of_interest))
all_copy_genes <- unique(c(columns_copy, paste0(genes_of_interest, "_CN")))

# Filter to keep only columns that actually exist in the data
columns_mutation <- intersect(colnames(alteration_data), all_mut_genes)
columns_copy <- intersect(colnames(alteration_data), all_copy_genes)

# Define the metadata columns we want to use for annotation
columns_meta <- c("treatment", "TMB")  # ("resistance", "type.of.primary", "prior.CTLA4i", "previous.BRAFi", "gender")

# Slice the main data frame into metadata and specific alteration types
meta_data <- alteration_data[, columns_meta]
mutation_data <- alteration_data[, columns_mutation]
copy_data <- alteration_data[, columns_copy]
loh_data <- alteration_data[, columns_loh]

# Clean the data: Convert string "NA" to actual NA values and ensure TMB is numeric
meta_data <- meta_data %>%
    mutate(across(where(is.character), ~na_if(., "NA"))) %>% # Convert "NA" strings to actual NA
    mutate(TMB = as.numeric(TMB)) # Ensure TMB is numeric

mutation_data <- mutation_data %>%
    mutate(across(everything(), ~na_if(., "NA")))

copy_data <- copy_data %>%
    mutate(across(everything(), ~na_if(., "NA")))

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

# Merge sub-variants: BRAF_600 into BRAF, NRAS_61 into NRAS
# We use a logical OR approach: if either is mutated, the merged gene is MUT.
if ("BRAF_600" %in% rownames(mut_matrix) && "BRAF" %in% rownames(mut_matrix)) {
  print("Merging BRAF_600 into BRAF mutation status...")
  mut_matrix["BRAF", ] <- ifelse(mut_matrix["BRAF", ] == "MUT" | mut_matrix["BRAF_600", ] == "MUT", "MUT", "")
}

if ("NRAS_61" %in% rownames(mut_matrix) && "NRAS" %in% rownames(mut_matrix)) {
  print("Merging NRAS_61 into NRAS mutation status...")
  mut_matrix["NRAS", ] <- ifelse(mut_matrix["NRAS", ] == "MUT" | mut_matrix["NRAS_61", ] == "MUT", "MUT", "")
}

# Create the copy number variation matrix (CNV)
# We map "Amplification" to "AMP" and "Deletion" to "DEL".
# Based on inspection, "hi amp" is amplification and "deep del" is deletion.
cnv_matrix <- t(copy_data) # Transpose
cnv_matrix[which(cnv_matrix == "hi amp", arr.ind = TRUE)] <- "AMP"
cnv_matrix[which(cnv_matrix == "deep del", arr.ind = TRUE)] <- "DEL"
# Set other non-alteration values to empty strings
cnv_matrix[!cnv_matrix %in% c("AMP", "DEL")] <- ""

# Merge sub-variants for CNV as well if they exist
if ("BRAF_600" %in% rownames(cnv_matrix) && "BRAF" %in% rownames(cnv_matrix)) {
  cnv_matrix["BRAF", ] <- ifelse(cnv_matrix["BRAF", ] != "" | cnv_matrix["BRAF_600", ] != "", 
                                 ifelse(cnv_matrix["BRAF", ] != "", cnv_matrix["BRAF", ], cnv_matrix["BRAF_600", ]), "")
}

if ("NRAS_61" %in% rownames(cnv_matrix) && "NRAS" %in% rownames(cnv_matrix)) {
  cnv_matrix["NRAS", ] <- ifelse(cnv_matrix["NRAS", ] != "" | cnv_matrix["NRAS_61", ] != "", 
                                 ifelse(cnv_matrix["NRAS", ] != "", cnv_matrix["NRAS", ], cnv_matrix["NRAS_61", ]), "")
}

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
    # resistance = c("primary" = "#E41A1C", "acquired" = "#377EB8", "primary (ipi acq.)" = "#4DAF4A", "NA" = "grey"),
    # type.of.primary = c("skin" = "#FDB462", "unknown" = "#B3B3B3", "mucosa" = "#FB8072"),
    # prior.CTLA4i = c("yes" = "black", "no" = "white"),
    # previous.BRAFi = c("no_" = "lightgrey", "no_During, day 7" = "#BEBADA", "yes_Before" = "#80B1D3"),
    # gender = c("F" = "pink", "M" = "lightblue"),
    # For the numeric TMB column, we create a continuous color function
    TMB = colorRamp2(c(0, quantile(meta_data$TMB, 0.95, na.rm = TRUE)), c("white", "darkviolet"))
)

ha_bottom <- HeatmapAnnotation(
    df = meta_data,
    col = anno_colors,
    na_col = "grey90", # Set a color for NA values in categorical annotations
    annotation_name_gp = gpar(fontsize = 8),
    simple_anno_size = unit(0.5, "cm")
)

# To save the output to a file, we wrap the plotting function in a graphics device.
jpeg("./Images/Molecular_pattern/oncoprint.jpg", width = 12, height = 8, units = "in", res = 300)
oncoPrint(combined_matrix, alter_fun = alter_fun, col = col, bottom_annotation = ha_bottom)
dev.off() # This closes the file and saves it.

print("--- OncoPrint successfully generated and saved as oncoprint.jpg ---")

# 1.4 Generate OncoPrint for Genes of Interest
print("--- 1.4: Generating OncoPrint for specific genes of interest ---")

# Slice the combined matrix for the genes of interest
# Only take genes that are actually in the matrix
genes_to_slice <- intersect(genes_of_interest, rownames(combined_matrix))
goi_matrix <- combined_matrix[genes_to_slice, , drop = FALSE]

# Generate and save the GOI OncoPrint
jpeg("./Images/Molecular_pattern/oncoprint_goi.jpg", width = 12, height = 8, units = "in", res = 300)
oncoPrint(goi_matrix, alter_fun = alter_fun, col = col, bottom_annotation = ha_bottom, 
          column_title = "Genomic Alterations in Genes of Interest")
dev.off()

print("--- Genes of Interest OncoPrint successfully generated and saved as oncoprint_goi.jpg ---")


# --- Part 2: Single Cell Data Analysis ---

print("--- Starting Part 2: Single Cell Analysis Workflow ---")

# 2.1 Load Single-Cell Data
print("--- 2.1: Loading single-cell data ---")
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

# 2.1.1 Explore Raw Cell Annotations
print("--- 2.1.1: Exploring the raw sc_annotation file ---")

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

# 2.1.2 Prepare and Align Metadata
print("--- 2.1.2: Preparing and aligning cell metadata ---")

# The annotation data should contain the cell barcodes.
if (!"CellBarcode" %in% colnames(sc_annotation)) {
    stop("The 'sc_annotation' data frame must contain a 'NAME' column with cell barcodes.")
}

# Set the row names of the annotation data to be the cell barcodes.
cell_metadata <- sc_annotation %>%
    tibble::column_to_rownames("CellBarcode")

summary(cell_metadata)


# Find the common cells that exist in both the count matrix and the metadata.
common_cells <- intersect(colnames(sc_counts), rownames(cell_metadata))

if (length(common_cells) == 0) {
    stop("No common cells found between count matrix and annotation file. Check barcode formats.")
}

# Align both the count matrix and the metadata to this common set of cells.
sc_counts_aligned <- sc_counts[, common_cells]
cell_metadata_aligned <- cell_metadata[common_cells, ]





# 2.2 Filter out B and T cells
print("--- 2.2: Filtering out B and T cells ---")

# Check if B/T cell columns exist
if (!all(c("Bcell", "Tcell") %in% colnames(cell_metadata))) {
   print("Warning: 'Bcell' or 'Tcell' columns not found in annotation. Checking available columns:")
   print(colnames(cell_metadata))
   # Fallback or stop if critical
} else {
    # Filter cells where Bcell == 0 AND Tcell == 0
    # The user specified: 0 (no) provided, so we keep 0.
    cells_to_keep <- rownames(cell_metadata[cell_metadata$Bcell == 0 & cell_metadata$Tcell == 0, ])
    print(paste("Original cell count:", nrow(cell_metadata)))
    print(paste("Cell count after removing B/T cells:", length(cells_to_keep)))
    
    # Subset aligned data
    sc_counts_aligned <- sc_counts_aligned[, cells_to_keep]
    cell_metadata_aligned <- cell_metadata_aligned[cells_to_keep, ]
}

sce <- SingleCellExperiment(
    assays = list(counts = as.matrix(sc_counts_aligned)),
    colData = cell_metadata_aligned
)

# Relabel the samples
print("--- Relabeling samples ---")
colData(sce)$Sample[colData(sce)$Sample == "Pat42"] <- "anti-PD1 res Pat42"
colData(sce)$Sample[colData(sce)$Sample == "Pat5"] <- "anti-CTLA4 res Pat5"
print("Relabeling complete. Checking unique samples:")
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
png("./Images/Molecular_pattern/sc_mean_variance_plot.png", width = 6, height = 5, units = "in", res = 300)
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

ggsave("./Images/Molecular_pattern/melanoma_umap_by_celltype.png", plot = p_umap, width = 10, height = 8)

print("UMAP plot saved to ./Images/Molecular_pattern/melanoma_umap_by_celltype.png")

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
if (!dir.exists("./Images/Molecular_pattern")) {
    dir.create("./Images/Molecular_pattern")
}

jpeg("./Images/Molecular_pattern/melanoma_marker_heatmap.jpg", width = 10, height = 12, units = "in", res = 300)
p_heatmap
dev.off()

print("Marker gene heatmap saved to ./Images/Molecular_pattern/melanoma_marker_heatmap.jpg")
print("--- Part 3: Single-cell analysis complete. ---")

# 2.8 Automated Cancer Cell Identification
print("--- 2.8: Automated Cancer Cell Identification ---")

# Define markers (Gene symbols)
marker_list <- list(
    Melanoma = c("MLANA", "PMEL", "MITF", "SOX10"),
    Endothelial = c("PECAM1", "VWF", "CLDN5"),
    CAF = c("COL1A1", "FAP", "ACTA2", "PDGFRB"),
    Myeloid = c("LYZ", "CD68", "CD14")
)

# Function to calculate average expression of markers per cluster
get_cluster_scores <- function(sce, markers, clusters) {
    unique_clusters <- unique(clusters)
    scores <- data.frame(Cluster = unique_clusters, row.names = unique_clusters)
    
    # Check valid markers
    valid_markers <- intersect(unlist(markers), rownames(sce))
    print(paste("Found", length(valid_markers), "valid markers out of", length(unlist(markers))))
    
    exprs_mat <- logcounts(sce)
    
    for (cell_type in names(markers)) {
        type_markers <- intersect(markers[[cell_type]], valid_markers)
        if (length(type_markers) > 0) {
            # For each cluster, calculate mean expression of these markers
            type_scores <- numeric(length(unique_clusters))
            for (i in seq_along(unique_clusters)) {
                clust <- unique_clusters[i]
                cells_in_clust <- names(clusters)[clusters == clust]
                
                # Mean of markers across all cells in cluster
                if (length(type_markers) == 1) {
                     # Handle single gene case
                     mean_expr <- mean(exprs_mat[type_markers, clusters == clust])
                } else {
                     # Handle multiple genes: mean of means or mean of all values
                     # Here take colMeans first (per cell) then mean of that
                     cell_means <- colMeans(exprs_mat[type_markers, clusters == clust, drop=FALSE])
                     mean_expr <- mean(cell_means)
                }
                type_scores[i] <- mean_expr
            }
            scores[[cell_type]] <- type_scores
        } else {
            scores[[cell_type]] <- 0
        }
    }
    return(scores)
}

# Run scoring
cluster_scores <- get_cluster_scores(sce, marker_list, colLabels(sce))
print("Mean expression scores per cluster:")
print(cluster_scores)

# Assign cell type to cluster based on max score
prob_cols <- names(marker_list)
cluster_scores$Assignment <- apply(cluster_scores[, prob_cols], 1, function(x) {
    if (max(x) < 0.1) return("Unknown/Other") # Threshold for very low expression
    prob_cols[which.max(x)]
})

print("Cluster Assignments:")
print(cluster_scores[, c("Cluster", "Assignment")])

# Map assignments back to cells
cluster_to_type <- setNames(cluster_scores$Assignment, cluster_scores$Cluster)
colData(sce)$CellType <- cluster_to_type[as.character(colLabels(sce))]

# Visualize automated annotation
p_umap_auto <- plotReducedDim(sce, "UMAP", colour_by = "CellType", text_by = "CellType") +
    ggtitle("UMAP by Automated Marker-based Annotation")
ggsave("./Images/Molecular_pattern/melanoma_umap_auto_annotation.png", plot = p_umap_auto, width = 10, height = 8)

print("Automated annotation plot saved to ./Images/Molecular_pattern/melanoma_umap_auto_annotation.png")

# 2.9 Subset Cancerous Cells
print("--- 2.9: Subsetting Cancerous Cells ---")
# Keep only Melanoma cells
sce_cancer <- sce[, colData(sce)$CellType == "Melanoma"]

print(paste("Original cell count (post-filter):", ncol(sce)))
print(paste("Cancerous cell count:", ncol(sce_cancer)))

if (ncol(sce_cancer) == 0) {
    warning("No cancerous cells identified! Checking thresholds or markers.")
} else {
    # Save the cancerous subset
    saveRDS(sce_cancer, "./Data/molecular_pattern/sce_cancer.rds")
    print("Cancerous cells subset saved to ./Data/molecular_pattern/sce_cancer.rds")
}

# --- Part 4: Differential Gene Expression Analysis (Melanoma Only) ---
print("--- Part 4: Differential Gene Expression Analysis on Melanoma Cells ---")

if (!exists("sce_cancer")) {
    stop("sce_cancer object not found. Please run the subsetting step first.")
}

# 4.1 Define Conditions
colData(sce_cancer)$Condition <- "Naive"
colData(sce_cancer)$Condition[grep("CTLA4", colData(sce_cancer)$Sample)] <- "CTLA4_Res"
colData(sce_cancer)$Condition[grep("PD1", colData(sce_cancer)$Sample)] <- "PD1_Res"

print("Condition counts in Melanoma subset:")
print(table(colData(sce_cancer)$Condition))

# 4.2 Perform DEA
print("--- Computing Differential Expression ---")

# Function to run DEA and save results
run_dea <- function(sce_obj, group1, group2, name) {
    print(paste("Computing DE:", group1, "vs", group2, "..."))
    subset_sce <- sce_obj[, sce_obj$Condition %in% c(group1, group2)]
    # Use findMarkers with wilcox test
    markers <- findMarkers(subset_sce, groups = subset_sce$Condition, test.type = "wilcox", direction = "any")
    stats <- markers[[group1]]
    
    # Convert to standard DF
    df <- as.data.frame(stats)
    
    # Check for logFC variations and standardize to 'logFC'
    # findMarkers often outputs 'summary.logFC' for multiple groups or 'logFC.<group>' for pairwise
    logfc_col <- grep("logFC", colnames(df), value = TRUE)
    if (length(logfc_col) > 0) {
        # Prefer summary.logFC if present
        if ("summary.logFC" %in% colnames(df)) {
             df$logFC <- df$summary.logFC
        } else {
             # Otherwise take the first logFC column available (e.g. logFC.Group2)
             df$logFC <- df[[logfc_col[1]]]
        }
    } else {
        # If no logFC column found at all (unlikely for findMarkers unless log.p=FALSE ?)
        # Manually calculate logFC if needed, or warn
        print("Warning: No logFC column found in findMarkers output. Calculating manual logFC...")
        # Simple manual calculation: mean(logcounts group1) - mean(logcounts group2)
        # Note: this is an approximation since we used wilcox test
        g1_cells <- which(subset_sce$Condition == group1)
        g2_cells <- which(subset_sce$Condition == group2)
        logcounts_mat <- logcounts(subset_sce)
        
        mean_g1 <- rowMeans(logcounts_mat[, g1_cells, drop=FALSE])
        mean_g2 <- rowMeans(logcounts_mat[, g2_cells, drop=FALSE])
        df$logFC <- mean_g1 - mean_g2
    }
    
    # Save CSV with standardized logFC
    filename <- paste0("./Data/molecular_pattern/DE_", name, "_Melanoma.csv")
    write.csv(df, filename)
    print(paste("Saved:", filename))
    return(df) # Return the standardized DF, not the raw stats object
}

stats_pd1 <- run_dea(sce_cancer, "PD1_Res", "Naive", "PD1_vs_Naive")
stats_ctla4 <- run_dea(sce_cancer, "CTLA4_Res", "Naive", "CTLA4_vs_Naive")

# 4.3 Volcano Plots
print("--- Generating Volcano Plots ---")

create_volcano_plot <- function(stats, title, filename, highlight_signatures = NULL) {
    df <- as.data.frame(stats)
    df$gene <- rownames(df)
    
    # Handle logFC column name variations
    if ("summary.logFC" %in% colnames(df)) {
        df$logFC <- df$summary.logFC
    }
    
    # Check if logFC is present
    if (!"logFC" %in% colnames(df)) {
         print(paste("Warning: 'logFC' (or 'summary.logFC') column not found for:", title))
         return(NULL)
    }
    
    # Ensure numeric columns
    if(is.character(df$logFC) | is.factor(df$logFC)) {
         df$logFC <- as.numeric(as.character(df$logFC))
    }

    df$p_val_adj <- df$FDR 
    
    # Determine Significance Direction (Red=Up, Blue=Down, Grey=NS)
    # BUT: If highlight_signatures is provided, we first set EVERYTHING to grey.
    df$ColorGroup <- "Not Significant"
    
    if (!is.null(highlight_signatures)) {
        # 1. Base State: Muted Background
        df$ColorGroup <- "Background"
        
        # 2. Identify Significant Treatment Signatures
        is_treatment <- df$gene %in% highlight_signatures
        is_sig <- df$FDR < 0.05
        is_up <- df$logFC > 0
        
        # Color Logic for Treatment Sigs
        df$ColorGroup[is_treatment & is_sig & is_up] <- "Upregulated"
        df$ColorGroup[is_treatment & is_sig & !is_up] <- "Downregulated"
        # Non-sig treatment or background remains "Background"
        
        # Define Colors
        color_map <- c("Upregulated" = "red", "Downregulated" = "blue", "Background" = "grey90")
        
        # Helper for plotting order (plot meaningful points LAST so they are on top)
        df$Order <- ifelse(df$ColorGroup == "Background", 1, 2)
        df <- df[order(df$Order), ]
        
    } else {
        # Default Behavior (No specific highlights)
        df$ColorGroup[df$FDR < 0.05 & df$logFC > 0.5] <- "Upregulated"
        df$ColorGroup[df$FDR < 0.05 & df$logFC < -0.5] <- "Downregulated"
        color_map <- c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")
    }
    
    # Base Plot
    p <- ggplot(df, aes(x = logFC, y = -log10(p_val_adj))) +
        geom_point(aes(color = ColorGroup), alpha = 1, size = 1.5) +
        scale_color_manual(values = color_map) +
        theme_bw() +
        labs(title = title, x = "Log2 Fold Change", y = "-Log10 FDR") +
        theme(legend.position = "none") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") # FDR 0.05 Line
        
    # Labelling Logic
    if (!is.null(highlight_signatures)) {
        # Label ONLY top 10 significant treatment signatures
        label_data <- df[df$ColorGroup %in% c("Upregulated", "Downregulated"), ]
        
        if (nrow(label_data) > 0) {
             # Take top 10 by significance (FDR)
             top_labels <- head(label_data[order(label_data$FDR), ], 10)
             
             p <- p + geom_text_repel(data = top_labels, 
                                      aes(label = gene), 
                                      max.overlaps = 50, box.padding = 0.5, size = 3)
        }
    } else {
        # Default Labeling (Top 10 overall)
        top_genes <- head(df[order(df$FDR), "gene"], 10)
        p <- p + geom_text_repel(data = df[df$gene %in% top_genes, ], aes(label = gene), max.overlaps = 20)
    }

    ggsave(filename, plot = p, width = 8, height = 6)
    print(paste("Saved volcano plot:", filename))
}

# Load Treatment Signatures for highlighting
treatment_sig_path <- "./Data/scData/post.treatment.sig.RData"
if (file.exists(treatment_sig_path)) {
  # Load the RData file (contains 'trt.sig')
  load(treatment_sig_path)
  
  # Combine UP and DOWN signatures from 'trt.sig'
  if (exists("trt.sig") && is.list(trt.sig)) {
      # Handle trt.up and trt.down
      up_sigs <- if("trt.up" %in% names(trt.sig)) trt.sig$trt.up else NULL
      down_sigs <- if("trt.down" %in% names(trt.sig)) trt.sig$trt.down else NULL
      
      treatment_sig_names <- unique(c(up_sigs, down_sigs))
      
      # Save signature names to text file
      sig_output_path <- "./Data/Results/Molecular_pattern/treatment_associated_signatures_list.txt"
      dir.create(dirname(sig_output_path), recursive = TRUE, showWarnings = FALSE)
      writeLines(treatment_sig_names, sig_output_path)
      
      message(paste("Loaded", length(treatment_sig_names), "unique treatment-associated signatures."))
      message(paste("Signature list saved to:", sig_output_path))
  } else {
      warning("Object 'trt.sig' not found or invalid format after loading.")
      treatment_sig_names <- NULL
  }
} else {
  warning("Treatment signature file not found (RData). Highlighting will be skipped.")
  treatment_sig_names <- NULL
}

# Generate Plots with Highlights
create_volcano_plot(stats_pd1, 
                    "Melanoma: Anti-PD1 Resistant vs Naive", 
                    "./Images/Molecular_pattern/Volcano_PD1_vs_Naive_Melanoma_Highlighted.jpg",
                    highlight_signatures = treatment_sig_names)

create_volcano_plot(stats_ctla4, 
                    "Melanoma: Anti-CTLA4 Resistant vs Naive", 
                    "./Images/Molecular_pattern/Volcano_CTLA4_vs_Naive_Melanoma_Highlighted.jpg",
                    highlight_signatures = treatment_sig_names)

# --- Save Checkpoint ---
# Save the current workspace 
save_path <- "./Data/Results/Molecular_pattern/checkpoint_pd1_ctla4_analysis.RData"
dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
message(paste("Saving workspace to:", save_path))
save.image(file = save_path)
message("Workspace saved. You can load this in future sessions using: load('", save_path, "')")


# 5. Load post treatment sig from public data of melanoma immune resistance
print("--- 5. Analysis of Post-Treatment Signature (Extraction) ---")

print("--- Loading post treatment sig ---")
# Ensure the file exists before loading
sig_file <- "./Data/scData/post.treatment.sig.RData"
if (file.exists(sig_file)) {
    load(sig_file)
    print("Signature loaded.")
    
    # Extract gene sets
    genes_up <- trt.sig$trt.up
    genes_down <- trt.sig$trt.down
    all_sig_genes <- unique(c(genes_up, genes_down))
    print(paste("Total signature genes:", length(all_sig_genes)))
    
    # Function to extract signature genes from stats and save
    extract_and_save_sig <- function(stats, sig_genes, up_genes, down_genes, name) {
        df <- as.data.frame(stats)
        
        # Match genes
        common_genes <- intersect(rownames(df), sig_genes)
        if (length(common_genes) == 0) {
            print(paste("Warning: No signature genes found in", name))
            return(NULL)
        }
        
        # Subset
        sig_df <- df[common_genes, ]
        sig_df$Gene <- rownames(sig_df)
        
        # Add Signature Annotation
        sig_df$Signature_Type <- ifelse(sig_df$Gene %in% up_genes, "Treatment_Up", 
                                  ifelse(sig_df$Gene %in% down_genes, "Treatment_Down", "Uncertain"))
        
        # Add Consistency Column
        # Check if Up matches positive LogFC or Down matches negative LogFC
        sig_df$Consistency <- ifelse(
            (sig_df$Signature_Type == "Treatment_Up" & sig_df$logFC > 0) | 
            (sig_df$Signature_Type == "Treatment_Down" & sig_df$logFC < 0), 
            "Yes", "No"
        )
        
        # Reorder columns (Gene, Signature, Consistency first)
        cols_to_front <- c("Gene", "Signature_Type", "Consistency", "logFC", "FDR")
        other_cols <- setdiff(colnames(sig_df), cols_to_front)
        sig_df <- sig_df[, c(cols_to_front, other_cols)]
        
        # Save
        filename <- paste0("./Data/molecular_pattern/DEA_Signature_", name, ".csv")
        write.csv(sig_df, filename, row.names = FALSE)
        print(paste("Saved signature subset:", filename))
        return(sig_df)
    }
    
    # Extract and Save for PD1 and CTLA4
    extract_and_save_sig(stats_pd1, all_sig_genes, genes_up, genes_down, "PD1_vs_Naive")
    extract_and_save_sig(stats_ctla4, all_sig_genes, genes_up, genes_down, "CTLA4_vs_Naive")

} else {
    print(paste("Error: Signature file not found at", sig_file))
}





















# 6. Signature Heatmap Visualization (Single Cell Expression)
print("--- 6. Generating Single Cell Signature Heatmap ---")

# Define conditions to include and set refined order
conditions_to_plot <- c("Naive", "PD1_Res", "CTLA4_Res")

# Subset sce_cancer for these conditions
sce_heatmap <- sce_cancer[, sce_cancer$Condition %in% conditions_to_plot]

# Set factor levels to ensure correct ordering (Naive -> PD1 -> CTLA4)
sce_heatmap$Condition <- factor(sce_heatmap$Condition, levels = conditions_to_plot)

# Explicitly order the SCE object by Condition to guarantee matrix order
sce_heatmap <- sce_heatmap[, order(sce_heatmap$Condition)]

# Extract signature genes
valid_genes <- intersect(all_sig_genes, rownames(sce_heatmap))

if (length(valid_genes) > 0) {
    print(paste("Plotting", length(valid_genes), "signature genes across", ncol(sce_heatmap), "cells."))
    
    # Get log-normalized counts
    expr_mat <- logcounts(sce_heatmap)[valid_genes, ]
    
    # Scale the matrix (Z-score by row/gene)
    expr_scaled <- t(scale(t(as.matrix(expr_mat))))
    expr_scaled[is.nan(expr_scaled)] <- 0
    expr_scaled[expr_scaled > 2] <- 2
    expr_scaled[expr_scaled < -2] <- -2

    # Column Annotation: Condition
    # Ensure factor levels are respected
    col_annot <- HeatmapAnnotation(
        Condition = factor(sce_heatmap$Condition, levels = conditions_to_plot),
        col = list(Condition = c("Naive" = "grey", "PD1_Res" = "red", "CTLA4_Res" = "blue")),
        show_annotation_name = FALSE
    )
    
    # Color function
    col_fun = circlize::colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red"))
    
    # Generate Heatmap
    print("Generating refined heatmap...")
    ht <- Heatmap(expr_scaled,
                  name = "Scaled Exp",
                  col = col_fun,
                  top_annotation = col_annot,
                  # No right_annotation (removed as requested)
                  column_split = sce_heatmap$Condition, # Split maintains the factor order if cluster_column_slices=FALSE
                  cluster_column_slices = FALSE, # This is key for ordering Naive -> PD1 -> CTLA4
                  cluster_rows = TRUE,
                  cluster_columns = FALSE, # Use split/factor order
                  show_row_names = FALSE,  
                  show_column_names = FALSE, 
                  column_title = c("Naive", "PD1_Res", "CTLA4_Res"), # Correctly label the splits
                  use_raster = TRUE 
    )
    
    # Save Plot
    pdf("./Images/Molecular_pattern/Signature_SingleCell_Heatmap.pdf", width = 8, height = 10)
    draw(ht)
    dev.off()
    print("Refined single cell heatmap saved to ./Images/Molecular_pattern/Signature_SingleCell_Heatmap.pdf")


# 8. Biological Function Analysis (Enrichment)
print("--- 8. Biological Function Analysis ---")

# Load required libraries for enrichment (if not already loaded)
if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    print("Loading clusterProfiler...")
    library(clusterProfiler)
}
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    print("Loading org.Hs.eg.db...")
    library(org.Hs.eg.db)
}
if (!requireNamespace("enrichplot", quietly = TRUE)) {
    print("Loading enrichplot...")
    library(enrichplot)
}

# Function to perform GO enrichment and save results
perform_enrichment <- function(genes, title_suffix, output_prefix) {
    print(paste("Running GO Enrichment for:", title_suffix))
    
    # 1. Convert Gene Symbols to Entrez IDs
    entrez_ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    
    if (nrow(entrez_ids) > 0) {
        # 2. Perform GO Enrichment (Biological Process)
        ego <- enrichGO(gene          = entrez_ids$ENTREZID,
                        OrgDb         = org.Hs.eg.db,
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.2,
                        readable      = TRUE) # Output gene symbols in result
        
        if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
            # 3. Save Results Table
            write.csv(as.data.frame(ego), 
                      file = paste0("./Data/Results/Molecular_pattern/GO_Enrichment_", output_prefix, ".csv"), 
                      row.names = FALSE)
            print(paste("Saved enrichment table to: ./Data/Results/Molecular_pattern/GO_Enrichment_", output_prefix, ".csv", sep=""))
            
            # 4. Generate Dotplot
            p_dot <- dotplot(ego, showCategory = 20) + 
                     ggtitle(paste("GO Enrichment:", title_suffix))
            
            ggsave(paste0("./Images/Molecular_pattern/GO_Dotplot_", output_prefix, ".pdf"), 
                   plot = p_dot, width = 8, height = 10)
            print(paste("Saved dotplot to: ./Images/Molecular_pattern/GO_Dotplot_", output_prefix, ".pdf", sep=""))
            
            # 5. Functional Classification (Simplify and Map)
            print("--- Simplifying and Classifying Functions ---")
            
            # Simplify GO terms to remove redundancy (cutoff=0.7)
            ego_sim <- tryCatch({
                # Explicitly call simplify from clusterProfiler
                clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
            }, error = function(e) {
                print(paste("Warning: Full simplify failed:", e$message))
                print("Trying simplified call...")
                tryCatch({
                    clusterProfiler::simplify(ego, cutoff=0.7)
                }, error = function(e2) {
                    print(paste("Warning: Simplify failed entirely:", e2$message))
                    return(ego) # Skip simplification
                })
            })
            
            # Ensure output directory exists before writing
            if (!dir.exists("./Data/Results/Molecular_pattern")) {
                dir.create("./Data/Results/Molecular_pattern", recursive = TRUE)
            }

            # Save Simplified Result
            if (!identical(ego_sim, ego)) {
                write.csv(as.data.frame(ego_sim), 
                          file = paste0("./Data/Results/Molecular_pattern/GO_Enrichment_Simplified_", output_prefix, ".csv"), 
                          row.names = FALSE)
                
                # Generate Dotplot for Simplified Terms
                p_dot_sim <- dotplot(ego_sim, showCategory = 15) + 
                             ggtitle(paste("Simplified GO Enrichment:", title_suffix))
                ggsave(paste0("./Images/Molecular_pattern/GO_Dotplot_Simplified_", output_prefix, ".pdf"), 
                       plot = p_dot_sim, width = 8, height = 8)
                
                # Use simplified result for classification
                class_source <- ego_sim
            } else {
                class_source <- ego
            }
            
            # Create Gene-Function Classification Table
            # Extract top 10 simplified terms
            top_terms <- as.data.frame(class_source)[1:min(10, nrow(class_source)), ]
            
            classification_list <- list()
            for (i in 1:nrow(top_terms)) {
                term_id <- top_terms$ID[i]
                term_desc <- top_terms$Description[i]
                genes_in_term <- strsplit(top_terms$geneID[i], "/")[[1]]
                
                if (length(genes_in_term) > 0) {
                    df_temp <- data.frame(
                        Gene = genes_in_term,
                        Biological_Function = term_desc,
                        GO_ID = term_id,
                        stringsAsFactors = FALSE
                    )
                    classification_list[[i]] <- df_temp
                }
            }
            
            if (length(classification_list) > 0) {
                classification_df <- do.call(rbind, classification_list)
                # Remove duplicates (a gene might be in multiple top terms, keep first/most significant occurrence or list all)
                # For "sectioning", usually unique assignment is preferred or we keep all to show multifunctionality.
                # Let's keep all for now so user can see all roles.
                
                write.csv(classification_df, 
                          file = paste0("./Data/Results/Molecular_pattern/Signature_Functional_Classification_", output_prefix, ".csv"), 
                          row.names = FALSE)
                print(paste("Saved functional classification to: ./Data/Results/Molecular_pattern/Signature_Functional_Classification_", output_prefix, ".csv", sep=""))
            }

            # 6. Functional Clustering (Enrichment Map & Treeplot)
            print("--- Generating Functional Clustering Plots ---")
            
            # Pairwise Term Similarity is required for both plots
            # Use simplified result if available, otherwise original
            target_ego <- if (exists("ego_sim") && !is.null(ego_sim)) ego_sim else ego
            
            # Calculate similarity
            target_ego <- tryCatch({
                enrichplot::pairwise_termsim(target_ego)
            }, error = function(e) {
                print(paste("Warning: pairwise_termsim failed:", e$message))
                return(NULL)
            })
            
            if (!is.null(target_ego)) {
                # Enrichment Map
                tryCatch({
                    p_emap <- emapplot(target_ego, showCategory = 30) + 
                              ggtitle(paste("Enrichment Map:", title_suffix))
                    ggsave(paste0("./Images/Molecular_pattern/Enrichment_Map_", output_prefix, ".pdf"), 
                           plot = p_emap, width = 12, height = 10)
                    print(paste("Saved Enrichment Map to: ./Images/Molecular_pattern/Enrichment_Map_", output_prefix, ".pdf", sep=""))
                }, error = function(e) {
                     print(paste("Warning: emapplot failed:", e$message))
                })
                
                # Treeplot (Hierarchical Clustering)
                tryCatch({
                    p_tree <- treeplot(target_ego, showCategory = 30) + 
                              ggtitle(paste("Functional Treeplot:", title_suffix))
                    ggsave(paste0("./Images/Molecular_pattern/Functional_Treeplot_", output_prefix, ".pdf"), 
                           plot = p_tree, width = 12, height = 8)
                    print(paste("Saved Treeplot to: ./Images/Molecular_pattern/Functional_Treeplot_", output_prefix, ".pdf", sep=""))
                }, error = function(e) {
                    print(paste("Warning: treeplot failed:", e$message))
                })
            }

        } else {
            print(paste("No significant enriched terms found for:", title_suffix))
        }
    } else {
        print(paste("No valid Entrez IDs found for:", title_suffix))
    }
}

# Ensure output directory exists
if (!dir.exists("./Data/Results/Molecular_pattern")) {
    dir.create("./Data/Results/Molecular_pattern", recursive = TRUE)
}

# Run Analysis for All Signature Genes
all_sig_genes_combined <- unique(c(genes_up, genes_down))

if (length(all_sig_genes_combined) > 0) {
    perform_enrichment(all_sig_genes_combined, "Post-Treatment Signature (All)", "Treatment_All")
} else {
    print("No signature genes to analyze.")
}

print("--- Analysis Complete ---")








# --- Moved Part 2: Bulk RNA Data Analysis (Non-executable) ---
# # 9. Generate a PCA plot
# 
# vsd <- vst(dds, blind=FALSE)
# pcaData <- plotPCA(vsd, intgroup="groups", returnData=TRUE)
# percentVar <- round(100 * attr(pcaData, "percentVar"))
# 
# ggplot(pcaData, aes(PC1, PC2, color=groups)) +
#   geom_point(size=3) +
#   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
#   ggtitle("PCA plot of samples") +
#   theme_bw()
# 
# # Save the PCA plot to a file
# ggsave("PCA_plot.jpg", width = 8, height = 6, units = "in", dpi = 300)
# 
# 
# 
# 
# 
# # --- Part 2: RNA Data Analysis ---
# # This section is for parsing and analyzing the RNA-seq data from the series matrix file.
# print("--- Part 2: Parsing RNA-seq metadata (placeholder for full analysis) ---")
# 
# # Load the raw series matrix file as lines
# rna_data_matrix_lines <- readLines("./Data/molecular_pattern/GSE244982_series_matrix.txt")
# 
# # Find the line with sample IDs
# sample_id_line_idx <- grep("^!Sample_title", rna_data_matrix_lines)
# if (length(sample_id_line_idx) > 0) {
#     sample_id_line <- rna_data_matrix_lines[sample_id_line_idx]
#     sample_ids_from_matrix <- strsplit(sample_id_line, "\t")[[1]][-1] # Split by tab, remove first element
#     sample_ids_from_matrix <- gsub('"', '', sample_ids_from_matrix)
# 
#     # Initialize a list to store metadata for each sample
#     parsed_meta_list <- vector("list", length(sample_ids_from_matrix))
#     names(parsed_meta_list) <- sample_ids_from_matrix
# 
#     # Find all lines containing sample characteristics
#     char_lines_idx <- grep("^!Sample_characteristics_ch1", rna_data_matrix_lines)
# 
#     # Process each characteristics line
#     for (line_idx in char_lines_idx) {
#         line <- rna_data_matrix_lines[line_idx]
#         elements <- strsplit(line, "\t")[[1]]
#         char_values <- elements[-1] # Actual values start from the second element
# 
#         for (i in seq_along(char_values)) {
#             sample_id <- sample_ids_from_matrix[i]
#             char_str <- char_values[i]
#             char_str_clean <- gsub('"', '', char_str) # Remove quotes
#             parts <- strsplit(char_str_clean, ": ", fixed = TRUE)[[1]] # Split by ": "
# 
#             if (length(parts) == 2) {
#                 key <- trimws(parts[1])
#                 value <- trimws(parts[2])
#                 parsed_meta_list[[sample_id]][[key]] <- value
#             } else {
#                 warning(paste("Skipping malformed or empty characteristic string:", char_str, "for sample:", sample_id))
#             }
#         }
#     }
# 
#     # Convert the list of lists into a data frame
#     all_keys <- unique(unlist(lapply(parsed_meta_list, names)))
#     rna_meta_df <- data.frame(
#         lapply(all_keys, function(key) {
#             sapply(parsed_meta_list, function(sample_meta) {
#                 val <- sample_meta[[key]]
#                 if (is.null(val)) NA else val
#             })
#         }),
#         stringsAsFactors = FALSE
#     )
#     colnames(rna_meta_df) <- all_keys
#     rownames(rna_meta_df) <- sample_ids_from_matrix
#     
#     # Clean the rownames of rna_meta_df by removing quotes
#     rownames(rna_meta_df) <- gsub('"', '', rownames(rna_meta_df))
# 
#     print("Successfully parsed RNA-seq metadata. First 5 rows:")
#     print(head(rna_meta_df, 5))
# 
//} else {
#    warning("Could not find '!Sample_geo_accession' line in GSE244982_series_matrix.txt. Skipping RNA metadata parsing.")
#}
# 
# # Placeholder for loading the actual RNA expression matrix from the same file.
# # The expression data typically starts after a line like '!series_matrix_table_begin'.
# # Future analysis would involve parsing this table and aligning it with the metadata above.
# 
# print("--- Part 1 complete. Part 2 is ready for further RNA data analysis implementation. ---")
# 
# # 2.1 Load RNA-seq Raw Counts
# RNAseq_raw <- read.table('./Data/molecular_pattern/GSE244982_ProcessedData_bulkRNAseq.txt', header = TRUE, stringsAsFactors = FALSE, 
# sep = '\t', quote = "", check.names = FALSE, row.names = 1)
# 
# # 2.2 Filter and Align Metadata with Expression Data
# print("--- 2.2: Filtering and aligning RNA-seq data ---")
# 
# # Filter out samples with "NA" group from the parsed metadata
# rna_meta_filtered <- rna_meta_df[!is.na(rna_meta_df$groups) & rna_meta_df$groups != "NA", ]
# 
# # Find common samples between metadata and expression data
# common_samples <- intersect(rownames(rna_meta_filtered), colnames(RNAseq_raw))
# 
# # Align both metadata and expression data to the common samples
# rna_meta_aligned <- rna_meta_filtered[common_samples, ]
# rnaseq_aligned <- RNAseq_raw[, common_samples]
# 
# # Ensure counts are integers for DESeq2
# rnaseq_aligned <- round(rnaseq_aligned)
# 
# # Verify alignment
# if (!all(rownames(rna_meta_aligned) == colnames(rnaseq_aligned))) {
#     stop("Alignment failed between RNA metadata and expression data columns.")
# } else {
#     print(paste("Successfully aligned", length(common_samples), "samples."))
# }
# 
# # 2.3 Prepare for Differential Expression Analysis (DEA)
# print("--- 2.3: Preparing data for Differential Expression Analysis ---")
# 
# # Subset the data to the three main comparison groups
# comparison_groups <- c("anti-CTLA4 resistant", "anti-PD1 resistant", "anti-PD1 resistant, BRAFi day7")
# meta_for_dea <- rna_meta_aligned[rna_meta_aligned$groups %in% comparison_groups, ]
# 
# # Make sure the 'groups' column is a factor for DESeq2
# meta_for_dea$groups <- factor(meta_for_dea$groups, levels = comparison_groups)
# 
# # Subset the expression matrix to match the DEA metadata
# rnaseq_for_dea <- rnaseq_aligned[, rownames(meta_for_dea)]
# 
# # Create a DESeqDataSet object
# dds <- DESeqDataSetFromMatrix(countData = rnaseq_for_dea,
#     colData = meta_for_dea,
#     design = ~ groups)
# 
# print("DESeq2 object 'dds' created successfully. Ready for differential expression analysis.")
# print("The 'mucosal' group is retained in 'rna_meta_aligned' for separate analysis (e.g., as an immune-cold reference).")
# 
# # 2.4 Perform Differential Expression Analysis (DEA)
# print("--- 2.4: Performing Differential Expression Analysis ---")
# 
# # Run DESeq2 without estimating size factors or applying VST/rlog
# dds <- estimateSizeFactors(dds, type = "poscounts")  # Use 'poscounts' to avoid errors with log-transformed data
# dds <- DESeq(dds, fitType = "parametric", sfType = "poscounts",
#     minReplicatesForReplace = Inf) # Specify fitType and sfType
# 
# print("DESeq2 analysis completed successfully.")
# 
# 
# # Function to extract and print results for a given contrast
# get_and_print_results <- function(dds, group1, group2) {
#     res <- results(dds, contrast = c("groups", group1, group2))
#     resOrdered <- res[order(res$pvalue), ]
#     print(paste("Top genes for:", group1, "vs", group2))
#     print(head(resOrdered, 10))
#     return(res)
# }
# 
# # Perform pairwise comparisons
# group1 <- "anti-CTLA4 resistant"
# group2 <- "anti-PD1 resistant"
# group3 <- "anti-PD1 resistant, BRAFi day7"
# 
# # anti-CTLA4 resistant vs. anti-PD1 resistant
# res_CTLA4_vs_PD1 <- get_and_print_results(dds, group1, group2)
# 
# # anti-CTLA4 resistant vs. anti-PD1 resistant, BRAFi day7
# res_CTLA4_vs_BRAFi <- get_and_print_results(dds, group1, group3)
# 
# # anti-PD1 resistant vs. anti-PD1 resistant, BRAFi day7
# res_PD1_vs_BRAFi <- get_and_print_results(dds, group2, group3)
# 
# 
# # Summary of the results for each group
# summary(res_CTLA4_vs_PD1)
# summary(res_CTLA4_vs_BRAFi)
# summary(res_PD1_vs_BRAFi)
# 
# # Count significant genes (p-value < 0.05) for each comparison
# sig_CTLA4_vs_PD1 <- sum(res_CTLA4_vs_PD1$pvalue < 0.05, na.rm = TRUE)
# sig_CTLA4_vs_BRAFi <- sum(res_CTLA4_vs_BRAFi$pvalue < 0.05, na.rm = TRUE)
# sig_PD1_vs_BRAFi <- sum(res_PD1_vs_BRAFi$pvalue < 0.05, na.rm = TRUE)
# 
# print(paste("Number of significant genes (p < 0.05) for", group1, "vs", group2, ":", sig_CTLA4_vs_PD1))
# print(paste("Number of significant genes (p < 0.05) for", group1, "vs", group3, ":", sig_CTLA4_vs_BRAFi))
# print(paste("Number of significant genes (p < 0.05) for", group2, "vs", group3, ":", sig_PD1_vs_BRAFi))
# 
# # Order results by adjusted p-value
# 
# 
# 
# 
# # Print the first 10 rows of the ordered results
# print("First 10 rows of the differential expression results:")
# print(head(resOrdered, 10))
# 
# 
# # 10. Identify signature genes for each group
# print("--- Identifying signature genes for each group ---")
# 
# num_signature_genes <- 50 # Number of signature genes to select per group
# 
# # Function to get top genes for a given comparison
# get_top_genes <- function(res, n = num_signature_genes) {
#     res_sig <- res[order(res$pvalue), ]
#     head(rownames(res_sig), n)
# }
# 
# # Get top genes for each comparison
# top_CTLA4_vs_PD1 <- get_top_genes(res_CTLA4_vs_PD1)
# top_CTLA4_vs_BRAFi <- get_top_genes(res_CTLA4_vs_BRAFi)
# top_PD1_vs_BRAFi <- get_top_genes(res_PD1_vs_BRAFi)
# 
# # Combine the top genes into a single list and remove duplicates
# signature_genes <- unique(c(top_CTLA4_vs_PD1, top_CTLA4_vs_BRAFi, top_PD1_vs_BRAFi))
# 
# # 11. Prepare data for heatmap
# print("--- Preparing data for heatmap ---")
# 
# # Extract the expression data for the signature genes
# heatmap_data <- rnaseq_aligned[signature_genes, ]
# 
# # Scale the data by row (gene)
# heatmap_data <- t(scale(t(heatmap_data)))
# 
# # Remove genes/samples with NA values
# heatmap_data <- heatmap_data[complete.cases(heatmap_data), complete.cases(t(heatmap_data))]
# 
# # 12. Generate Heatmap
# print("--- Generating heatmap ---")
# 
# # Set color option
# hmcols <- colorRampPalette(brewer.pal(11, "RdBu"))(256)
# 
# # Order columns by groups
# column_order <- order(rna_meta_aligned$groups)
# 
# 
# 
# #Define colors for the heatmap annotation
# anno_colors <- list(
#     groups = c("anti-CTLA4 resistant" = "#8DD3C7",
#                "anti-PD1 resistant" = "#FFFFB3",
#                "anti-PD1 resistant" = "#FFFFB3", 
#                "mucosal" = "#A6CEE3",
#                "anti-PD1 resistant, BRAFi day7" = "#BEBADA")
# )
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # Define annotation for the heatmap
# ha <- HeatmapAnnotation(
# 
#   groups = rna_meta_aligned$groups,
# 
#   col = list(groups = anno_colors$groups)
# )
# 
# # Create the heatmap
# Heatmap(heatmap_data, 
#         name = "Z-score", # Title of the heatmap
#         col = hmcols, # Colors for the heatmap
#         top_annotation = ha, # Add the top annotation
#         show_row_names = FALSE, # Hide row names
#         show_column_names = FALSE, # Hide column names
#         cluster_rows = TRUE, # Cluster rows
#         cluster_columns = FALSE, # Do not cluster columns
#         column_order = column_order) # Order columns by group
# 
# # Save the heatmap to a file
# jpeg("heatmap_signature_genes.jpg", width = 10, height = 8, units = "in", res = 300)
# Heatmap(heatmap_data, 
#        name = "Z-score", # Title of the heatmap
#         col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")), # Colors for the heatmap
#         top_annotation = ha, # Add the top annotation
#         show_row_names = FALSE, # Hide row names
#         show_column_names = FALSE, # Hide column names
#         cluster_rows = TRUE, # Cluster rows
#         cluster_columns = FALSE, # Do not cluster columns
#         column_order = column_order) # Order columns by group
# dev.off()
# 
# # --- Additional diagnostic checks ---
# 
# # 1. Check the distribution of p-values
# hist(res$pvalue, main = "Distribution of p-values", xlab = "P-value")
# 
# # 2. Check the size factors
# print("Size factors:")
# print(sizeFactors(dds))
# 
# # 3. Check the dispersion estimates
# 
# print("Dispersion estimates:")
# print(dispersions(dds))
# 
# # 4. Check the contrast
# print("The contrast used:")
# print(resultsNames(dds))
# 
# 
# 
# # 13. Perform GO term enrichment analysis
# print("--- Performing GO term enrichment analysis ---")
# 
# # Install and load necessary packages
# #if (!requireNamespace("BiocManager", quietly = TRUE))
# #    install.packages("BiocManager")
# #BiocManager::install("clusterProfiler")
# #BiocManager::install("org.Hs.eg.db")
# library(clusterProfiler)
# library(org.Hs.eg.db)
# 
# # Function to perform GO enrichment analysis
# perform_go_enrichment <- function(gene_list, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05) {
#   ego <- enrichGO(gene          = gene_list,
#                   OrgDb         = org.Hs.eg.db,
#                   keyType       = "SYMBOL",
#                   ont           = ont,
#                   pAdjustMethod = pAdjustMethod,
#                   pvalueCutoff  = pvalueCutoff,
#                   qvalueCutoff  = qvalueCutoff,
#                   readable      = TRUE)
#   return(ego)
# }
# 
# # Perform GO enrichment analysis for each signature
# ego_CTLA4_vs_PD1 <- perform_go_enrichment(top_CTLA4_vs_PD1)
# ego_CTLA4_vs_BRAFi <- perform_go_enrichment(top_CTLA4_vs_BRAFi)
# ego_PD1_vs_BRAFi  <- perform_go_enrichment(top_PD1_vs_BRAFi)
# 
# # Print results for each signature
# print("GO enrichment results for CTLA4_vs_PD1:")
# print(head(summary(ego_CTLA4_vs_PD1), 10))
# 
# print("GO enrichment results for CTLA4_vs_BRAFi:")
# print(head(summary(ego_CTLA4_vs_BRAFi), 10))
# 
# print("GO enrichment results for PD1_vs_BRAFi:")
# print(head(summary(ego_PD1_vs_BRAFi), 10))
# 
# # You can visualize the results using dotplot
# dotplot(ego_CTLA4_vs_PD1, showCategory=10)
# dotplot(ego_CTLA4_vs_BRAFi, showCategory=10)
# dotplot(ego_PD1_vs_BRAFi, showCategory=10)


