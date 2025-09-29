# A dataset from the pubulication of <<Molecular patterns of resistance to immune checkpoint blockade in melanoma>>
library(ComplexHeatmap)
library(circlize) # For colorRamp2
library(dplyr) # For data manipulation

# --- Part 1: Alteration Data and OncoPrint Generation ---

# 1.1 Load and Prepare Alteration Data
print("--- 1.1: Loading and preparing alteration data from Supplementary.xlsx ---")
alteration_data <- read_excel("./Data/molecular_pattern/Supplementary.xlsx", sheet = 2)

# The first column is the sample ID. Let's set it as the row names for proper alignment.
if (!is.character(alteration_data[[1]])) {
    stop("The first column is not a character vector suitable for sample IDs.")
}
rownames(alteration_data) <- alteration_data[[1]]

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

# Slice the main data frame into specific alteration types
mutation_data <- alteration_data[, columns_mutation]
copy_data <- alteration_data[, columns_copy]
loh_data <- alteration_data[, columns_loh]

# Convert string "NA" to actual NA values for robust processing
mutation_data[mutation_data == "NA"] <- NA
copy_data[copy_data == "NA"] <- NA
loh_data[loh_data == "NA"] <- NA

# 1.2 Parse, Structure, and Align Metadata
print("--- 1.2: Parsing and aligning metadata from GSE244982_series_matrix.txt ---")

# Load the raw series matrix file as lines
rna_data_matrix_lines <- readLines('./Data/molecular_pattern/GSE244982_series_matrix.txt')

# Find the line with sample IDs (GEO accession numbers)
sample_id_line_idx <- grep("^!Sample_geo_accession", rna_data_matrix_lines)
if (length(sample_id_line_idx) == 0) {
    stop("Could not find '!Sample_geo_accession' line in GSE244982_series_matrix.txt")
}
sample_id_line <- rna_data_matrix_lines[sample_id_line_idx]
sample_ids_from_matrix <- strsplit(sample_id_line, "\t")[[1]][-1] # Split by tab, remove first element

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
            # Handle cases where a characteristic might be malformed or empty
            warning(paste("Skipping malformed or empty characteristic string (no clear key:value format):", char_str, "for sample:", sample_id))
        }
    }
}

# Convert the list of lists into a data frame
all_keys <- unique(unlist(lapply(parsed_meta_list, names)))
structured_rna_meta_df <- data.frame(
    lapply(all_keys, function(key) {
        sapply(parsed_meta_list, function(sample_meta) {
            val <- sample_meta[[key]]
            if (is.null(val)) NA else val
        })
    }),
    stringsAsFactors = FALSE
)
colnames(structured_rna_meta_df) <- all_keys
rownames(structured_rna_meta_df) <- sample_ids_from_matrix

# The oncoprint samples are defined by the rownames of `alteration_data`.
oncoprint_samples <- rownames(alteration_data)

# Subset and reorder the parsed metadata to match the oncoprint samples
# If a sample from oncoprint_samples is not found in structured_rna_meta_df, its metadata will be NA.
meta_data_temp <- structured_rna_meta_df[match(oncoprint_samples, rownames(structured_rna_meta_df)), , drop = FALSE]
rownames(meta_data_temp) <- oncoprint_samples # Ensure rownames are correct after match

# Define the metadata columns we want to use for annotation
columns_meta <- c("treatment", "resistance", "type.of.primary", "prior.CTLA4i", "previous.BRAFi", "gender", "TMB")
# Select only the columns specified in `columns_meta`
missing_meta_cols <- setdiff(columns_meta, colnames(meta_data_temp))
if (length(missing_meta_cols) > 0) {
    warning(paste("The following metadata columns were requested but not found in GSE244982_series_matrix.txt:", paste(missing_meta_cols, collapse = ", ")))
    for (col_name in missing_meta_cols) {
        meta_data_temp[[col_name]] <- NA # Add missing columns as NA
    }
}

# Ensure the order of columns in meta_data_temp matches columns_meta
meta_data <- meta_data_temp[, columns_meta, drop = FALSE]

# Convert string "NA" to actual NA values and ensure TMB is numeric for plotting
meta_data <- meta_data %>%
    mutate(across(where(is.character), ~na_if(., "NA"))) %>% # Convert "NA" strings to actual NA
    mutate(TMB = as.numeric(TMB)) # Ensure TMB is numeric

# 1.3 Prepare Data Matrices for OncoPrint
print("--- 1.3: Preparing data matrices for OncoPrint visualization ---")
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

# 1.4 Define Aesthetics and Generate OncoPrint
print("--- 1.4: Defining aesthetics and generating OncoPrint ---")
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

# This section is reserved for loading and analyzing the RNA-seq data from the series matrix file.
print("--- Part 1 complete. Part 2 (RNA Data Analysis) can be implemented below. ---")
