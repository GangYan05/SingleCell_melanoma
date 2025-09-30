# A dataset from the pubulication of <<Molecular patterns of resistance to immune checkpoint blockade in melanoma>>
library(ComplexHeatmap)
library(circlize) # For colorRamp2
library(dplyr) # For data manipulation
library(DESeq2) # For differential expression analysis

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




# --- Part 2: RNA Data Analysis ---
# This section is for parsing and analyzing the RNA-seq data from the series matrix file.
print("--- Part 2: Parsing RNA-seq metadata (placeholder for full analysis) ---")

# Load the raw series matrix file as lines
rna_data_matrix_lines <- readLines('./Data/molecular_pattern/GSE244982_series_matrix.txt')

# Find the line with sample IDs (GEO accession numbers)
sample_id_line_idx <- grep("^!Sample_geo_accession", rna_data_matrix_lines)
if (length(sample_id_line_idx) > 0) {
    sample_id_line <- rna_data_matrix_lines[sample_id_line_idx]
    sample_ids_from_matrix <- strsplit(sample_id_line, "\t")[[1]][-1] # Split by tab, remove first element
    # Clean the sample IDs by removing quotes
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
