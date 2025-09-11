# Function to install packages (Bioconductor and CRAN)
install_and_load_packages <- function(bioc_pkgs = NULL, cran_pkgs = NULL) {
    if (!is.null(bioc_pkgs)) {
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
        for (pkg in bioc_pkgs) {
          if (!requireNamespace(pkg, quietly = TRUE)) {
            BiocManager::install(pkg)
          } else {
            print(paste("package", pkg, "already installed"))
          }
          library(pkg, character.only = TRUE)
        }
    }

    if (!is.null(cran_pkgs)) {
      for (pkg in cran_pkgs) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            install.packages(pkg)
          } else {
            print(paste("package", pkg, "already installed"))
          }
         library(pkg, character.only = TRUE)
      }
    }
}



bioc_packages <- c("SingleCellExperiment", "scuttle", "scran", "scater", "uwot", 
                   "rtracklayer", "DropletUtils", "batchelor", "bluster", "ensembldb", 
                   "org.Mm.eg.db", "org.Hs.eg.db", "DropletTestFiles", "scRNAseq", "AnnotationHub",
                   "PCAtools", "celldex", "SingleR", "TENxPBMCData", "depmap")
cran_packages <- c("uwot", "dynamicTreeCut", "dplyr", "pheatmap", "Seurat")

install_and_load_packages(bioc_pkgs = bioc_packages, cran_pkgs = cran_packages)
