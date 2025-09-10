# ImmRes1_decomposed_walkthrough.R
#
# This script breaks down the analysis from `ImmRes1_denovoCellTypeSig.R` into
# sequential, runnable sections. Each section has defined inputs and outputs,
# allowing you to understand and adapt the workflow for your own analysis.
#
# INSTRUCTIONS:
# 1. Make sure this script is in the same directory as the other ImmRes scripts.
# 2. Run each section sequentially.
# 3. After each section, inspect the output objects (e.g., using `str()`, `head()`, `summary()`)
#    to understand what was generated.

################################################################################
### PREAMBLE: Load Libraries and Source Functions
################################################################################

# This step loads all necessary packages and custom functions from the project.

print("Loading libraries and sourcing project files...")
source('ImmRes_source.R')


################################################################################
### SECTION 1: Load Data and Perform Pairwise Differential Expression
################################################################################

# --- PURPOSE ---
# Load the main single-cell dataset and perform all pairwise differential
# expression (DE) tests between the major cell types. This is the most
# computationally intensive step.

# --- INPUT ---
# - A pre-processed single-cell RNA-seq object.
#   File: ../Data/scData/Mel.all.data.QC.rds

# --- OUTPUT ---
# - `r`: The loaded single-cell data object.
# - `de_results`: A list containing the Z-scores from all pairwise comparisons.

print("--- SECTION 1: Loading data and running pairwise DE ---")

# Load the main single-cell dataset
r <- readRDS("../Data/scData/Mel.all.data.QC.rds")

# --- INSPECT THE 'r' OBJECT ---
print("--- Structure of the loaded 'r' object: ---")
str(r, max.level = 1) # Use max.level=1 to see just the top-level elements

# Get all unique cell types and create all possible pairs for comparison
cell.types.u <- sort(unique(r$cell.types))
cell.type.pairs <- t(combn(unique(r$cell.types), 2))

# Perform pairwise DE analysis (Wilcoxon rank-sum test) for each pair
de_zscores <- apply(cell.type.pairs, 1, function(x) {
  print(paste("Comparing", x[1], "to", x[2]))
  b1 <- is.element(r$cell.types, x[1])
  b2 <- is.element(r$cell.types, x[2])
  de <- apply.ttest(r$tpm[, b1 | b2], b1[b1 | b2], ranksum.flag = TRUE)[, "zscores"]
  return(de)
})
rownames(de_zscores) <- r$genes
colnames(de_zscores) <- paste(cell.type.pairs[, 1], cell.type.pairs[, 2], sep = "_")

# Store results in a list
de_results <- list(
  cell.type.pairs = cell.type.pairs,
  cell.types = cell.types.u,
  de.zscores = de_zscores
)

print("--- Section 1 Complete. Inspect the 'de_results' object. ---")


################################################################################
### SECTION 2: Summarize DE and Generate Initial Cell Type Signatures
################################################################################

# --- PURPOSE ---
# For each cell type, summarize the DE results to find genes that are
# consistently upregulated against all other cell types.

# --- INPUT ---
# - `de_results`: The list containing pairwise Z-scores from Section 1.

# --- OUTPUT ---
# - `cell_type_de`: An updated list containing summarized DE scores (`de.sum`)
#   and the initial gene signatures (`sig`).

print("--- SECTION 2: Summarizing DE and creating initial signatures ---")

# This code is taken directly from the `get.cell.type.sig` function.
MIN_Z_SCORE_THRESHOLD <- 5

de_results$de.sum <- lapply(cell.types.u, function(x) {
  b1 <- is.element(de_results$cell.type.pairs[, 1], x)
  b2 <- is.element(de_results$cell.type.pairs[, 2], x)
  de <- cbind.data.frame(de_results$de.zscores[, b1], -de_results$de.zscores[, b2])
  colnames(de) <- c(de_results$cell.type.pairs[b1, 2], de_results$cell.type.pairs[b2, 1])
  de$Min <- rowMins(as.matrix(de))
  return(de)
})
names(de_results$de.sum) <- cell.types.u

de_results$sig <- lapply(de_results$de.sum, function(x) sort(rownames(x)[x$Min > MIN_Z_SCORE_THRESHOLD]))
de_results$sig <- remove.ribo(de_results$sig) # Remove ribosomal genes

cell_type_de <- de_results

print("--- Section 2 Complete. Inspect 'cell_type_de$sig' to see the first signatures. ---")


################################################################################
### SECTION 3: Refine T-Cell and Generate Supertype Signatures
################################################################################

# --- PURPOSE ---
# Run helper functions to create more specific T-cell signatures (e.g., a
# general T-cell signature) and signatures for broad "supertypes" like
# 'lymphocyte' and 'stroma'.

# --- INPUT ---
# - `cell_type_de`: The results object from Section 2.

# --- OUTPUT ---
# - `cell_type_de_full`: The final object containing all major and supertype
#   signatures. This object is also saved to a file.

print("--- SECTION 3: Refining T-cell and creating supertype signatures ---")

cell_type_de <- get.t.cell.sig(cell_type_de)
cell_type_de_full <- get.cell.supertype.sig(de = cell_type_de)

print("--- Section 3 Complete. Inspect 'cell_type_de_full$sig'. ---")
print("Intermediate results saved to '../Results/CellTypes/cell.type.sig.full.RData'")


################################################################################
### SECTION 4: Generate T-Cell Functional Subset Signatures
################################################################################

# --- PURPOSE ---
# Perform a deep dive into T-cells to identify signatures for functional states
# like exhaustion, memory, and regulatory T-cells. This function is a wrapper
# that loads T-cell specific data and runs a new round of DE analysis.

# --- INPUT ---
# - T-cell specific data files (loaded automatically by the function).

# --- OUTPUT ---
# - `final_cell_signatures`: The complete collection of all cell type, supertype,
#   and T-cell subset signatures. This is the main output of the entire script.

print("--- SECTION 4: Generating T-cell functional subset signatures ---")

# This function handles its own data loading.
final_cell_signatures <- get.all.t.subset.sig(rF = r, r4 = NULL, r8 = NULL)

print("--- Section 4 Complete. Inspect 'final_cell_signatures'. This is the final signature list. ---")
print("Final results saved to '../Results/Signatures/cell.type.sig.full.RData'")


################################################################################
### SECTION 5: (Optional) Compute Signature Scores and Generate Plots
################################################################################

# --- PURPOSE ---
# Score every cell in the original dataset against the final signatures and
# generate t-SNE plots to visualize the expression of these signatures.

# --- INPUT ---
# - `r`: The full single-cell object.
# - `final_cell_signatures`: The final signature list from Section 4.

# --- OUTPUT ---
# - `type_OE`: A matrix of signature scores for each cell.
# - PDF files with t-SNE plots in `../Output/Figures/`.

print("--- SECTION 5: Computing signature scores and generating plots ---")

type_OE <- compute.cell.type.OE(r = r, cell.sig = final_cell_signatures)

print("--- Section 5 Complete. Analysis finished. ---")
print("Signature scores are in 'type_OE' and plots are in the Output/Figures directory.")


################################################################################
### SECTION 6: (Optional) Convert to Standard Single-Cell Object Classes
################################################################################

# --- PURPOSE ---
# Convert the custom 'r' list object into standard, widely-used single-cell
# object classes like Seurat or SingleCellExperiment for better interoperability
# with other bioinformatics tools.

# --- INPUT ---
# - `r`: The full single-cell object from Section 1.

# --- OUTPUT ---
# - `seurat_obj`: A Seurat object containing the data.
# - `sce_obj`: A SingleCellExperiment object containing the data.

print("--- SECTION 6: Converting to standard object classes ---")

# --- Option A: Convert to a Seurat Object ---

# First, ensure the Seurat package is installed.
# if (!requireNamespace("Seurat", quietly = TRUE)) {
#   install.packages("Seurat")
# }
library(Seurat)

print("Creating Seurat object...")

# Note: Seurat prefers raw counts, but we use the provided TPM matrix here.
seurat_obj <- CreateSeuratObject(counts = r$tpm, project = "MelanomaImmuno")

# Add cell-level metadata
seurat_obj$cell.types <- r$cell.types
seurat_obj$samples <- r$samples

# Add the t-SNE dimensional reduction
tsne_reduc <- CreateDimReducObject(embeddings = r$tsne, key = "tSNE_", assay = DefaultAssay(seurat_obj))
seurat_obj[["tsne"]] <- tsne_reduc

print("Seurat object created. Inspect 'seurat_obj'.")


# --- Option B: Convert to a SingleCellExperiment Object ---

# First, ensure the SingleCellExperiment package is installed.
# if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
#   BiocManager::install("SingleCellExperiment")
# }
library(SingleCellExperiment)

print("Creating SingleCellExperiment object...")
sce_obj <- SingleCellExperiment(assays = list(tpm = r$tpm),
                                colData = data.frame(cell.types = r$cell.types, samples = r$samples),
                                reducedDims = list(TSNE = r$tsne))

print("SingleCellExperiment object created. Inspect 'sce_obj'.")
print("--- Section 6 Complete. ---")
