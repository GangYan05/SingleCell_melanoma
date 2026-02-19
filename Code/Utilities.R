# Function to install packages (Bioconductor and CRAN)
install_and_load_packages <- function(bioc_pkgs = NULL, cran_pkgs = NULL) {
    # Install and load CRAN packages
    if (!is.null(cran_pkgs) && length(cran_pkgs) > 0) {
        new_pkgs <- cran_pkgs[!(cran_pkgs %in% installed.packages()[, "Package"])]
        if (length(new_pkgs) > 0) {
            message("Installing missing CRAN packages: ", paste(new_pkgs, collapse = ", "))
            install.packages(new_pkgs)
        } else {
            message("All required CRAN packages are already installed.")
        }
        sapply(cran_pkgs, require, character.only = TRUE)
    }

    # Install and load Bioconductor packages
    if (!is.null(bioc_pkgs) && length(bioc_pkgs) > 0) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
            install.packages("BiocManager")
        }
        new_pkgs <- bioc_pkgs[!(bioc_pkgs %in% installed.packages()[, "Package"])]
        if (length(new_pkgs) > 0) {
            message("Installing missing Bioconductor packages: ", paste(new_pkgs, collapse = ", "))
            # Use ask = FALSE and update = FALSE for non-interactive execution
            BiocManager::install(new_pkgs, ask = FALSE, update = FALSE)
        } else {
            message("All required Bioconductor packages are already installed.")
        }
        sapply(bioc_pkgs, require, character.only = TRUE)
    }
}



bioc_packages <- c("SingleCellExperiment", "scuttle", "scran", "scater", "rtracklayer", 
                   "DropletUtils", "batchelor", "bluster", "ensembldb", "DESeq2", "clusterProfiler",
                   "org.Mm.eg.db", "org.Hs.eg.db", "DropletTestFiles", "scRNAseq", "AnnotationHub",
                   "PCAtools", "celldex", "SingleR", "TENxPBMCData", "depmap", "ComplexHeatmap")
cran_packages <- c("uwot", "dynamicTreeCut", "dplyr", "pheatmap", "Seurat", "ggplot2", 
                   "gridExtra", "RColorBrewer", "tidyr", "tibble", "magrittr", "ggrepel",
                   "circlize", "grid", "readxl", "corrplot")

install_and_load_packages(bioc_pkgs = bioc_packages, cran_pkgs = cran_packages)
