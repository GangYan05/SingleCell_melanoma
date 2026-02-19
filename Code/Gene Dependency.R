# Gene Dependency Analysis Script
library(dplyr)
library(tibble)
library(pheatmap)
library(RColorBrewer)
library(tidyr)
library(ggplot2)
library(circlize)

# =============================================================================
# Module 1: Data Loading and Preprocessing
# =============================================================================

load_and_annotate_depmap <- function(crispr_path, model_path, essential_path) {
  message("Loading and annotating DepMap data...")
  
  # 1. Load Files
  crispr_data <- read.csv(crispr_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  depmap_meta <- read.csv(model_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  essential_list <- read.csv(essential_path, stringsAsFactors = FALSE, header = TRUE)
  
  # 2. Gene Name Cleanup (Remove " (ID)")
  colnames(crispr_data) <- gsub(" \\(.*\\)", "", colnames(crispr_data))
  essential_clean <- gsub(" \\(.*\\)", "", essential_list$gene)
  
  # 3. Align and Merge metadata
  common_cell_lines <- intersect(rownames(crispr_data), rownames(depmap_meta))
  message(sprintf("Found %d common cell lines.", length(common_cell_lines)))
  
  depmap_annotated <- crispr_data[common_cell_lines, ] %>%
    rownames_to_column("ModelID") %>%
    left_join(
      depmap_meta[common_cell_lines, ] %>%
        select(OncotreeLineage, OncotreeSubtype) %>%
        rownames_to_column("ModelID"),
      by = "ModelID"
    ) %>%
    select(ModelID, OncotreeLineage, OncotreeSubtype, everything())
  
  list(
    data = depmap_annotated,
    essential_genes = essential_clean
  )
}

# =============================================================================
# Module 2: Essential Genes Analysis
# =============================================================================

run_essential_genes_module <- function(depmap_obj, out_img_dir) {
  message("\n--- Running Essential Genes Module ---")
  
  # 1. Filter for common essential genes
  shared_essential <- intersect(depmap_obj$essential_genes, colnames(depmap_obj$data))
  message(sprintf("Found %d shared essential genes.", length(shared_essential)))
  
  if(length(shared_essential) == 0) {
    stop("No shared essential genes found. Check column name format.")
  }
  
  essential_data <- depmap_obj$data %>% select(ModelID, OncotreeLineage, all_of(shared_essential))
  
  # 2. Lineage Average
  lineage_scores <- essential_data %>%
    group_by(OncotreeLineage) %>%
    summarise(across(all_of(shared_essential), ~ mean(.x, na.rm = TRUE)))
  
  # 3. Specificity Rank (Melanoma/Skin Specificity)
  lineage_counts <- count(depmap_obj$data, OncotreeLineage)
  valid_lineages <- lineage_counts$OncotreeLineage[lineage_counts$n >= 5]
  
  melanoma_spec <- lineage_scores %>%
    filter(OncotreeLineage %in% valid_lineages) %>%
    pivot_longer(cols = -OncotreeLineage, names_to = "gene", values_to = "score") %>%
    group_by(gene) %>%
    summarise(
      skin_score = score[OncotreeLineage == "Skin"],
      others_mean = mean(score[OncotreeLineage != "Skin"], na.rm = TRUE)
    ) %>%
    mutate(specificity = skin_score - others_mean) %>%
    arrange(specificity)
  
  # Save specificity results
  write.csv(melanoma_spec, file.path("./Data/Results/DepMap", "melanoma_essential_specificity.csv"), row.names = FALSE)
  
  # 4. Visualization: Heatmap
  heatmap_mat <- lineage_scores %>%
    filter(OncotreeLineage %in% valid_lineages) %>%
    column_to_rownames("OncotreeLineage") %>%
    as.matrix()
    
  pheatmap(
    heatmap_mat,
    main = "Average Essentiality of Common Essential Genes Across Lineages",
    show_colnames = FALSE,
    color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
    filename = file.path(out_img_dir, "lineage_essentiality_heatmap.png"),
    width = 10, height = 8
  )
  
  message("Essential Genes Heatmap saved.")
  return(melanoma_spec)
}

# =============================================================================
# Module 3: Genes of Interest (Inquiry) Module
# =============================================================================

run_genes_of_interest_module <- function(depmap_obj, genes_of_interest, out_img_dir, out_res_dir) {
  message("\n--- Running Genes of Interest Inquiry Module ---")
  
  # 1. Subset for Genes of Interest
  valid_genes <- intersect(genes_of_interest, colnames(depmap_obj$data))
  if(length(valid_genes) == 0) {
    message("None of the genes of interest found in DepMap data.")
    return(NULL)
  }
  
  goi_data <- depmap_obj$data %>% 
    select(ModelID, OncotreeLineage, all_of(valid_genes)) %>%
    mutate(lineage_group = ifelse(OncotreeLineage == "Skin", "Skin", "Other"))
  
  # 2. Stats: Skin vs Others
  stats_results <- data.frame()
  for (g in valid_genes) {
    skin_vals <- goi_data %>% filter(lineage_group == "Skin") %>% pull(!!sym(g))
    other_vals <- goi_data %>% filter(lineage_group == "Other") %>% pull(!!sym(g))
    
    skin_vals <- skin_vals[!is.na(skin_vals)]
    other_vals <- other_vals[!is.na(other_vals)]
    
    if (length(skin_vals) > 1 & length(other_vals) > 1) {
      wt <- wilcox.test(skin_vals, other_vals)
      stats_results <- rbind(stats_results, data.frame(
        Gene = g,
        Skin_Mean = mean(skin_vals),
        Others_Mean = mean(other_vals),
        Diff = mean(skin_vals) - mean(other_vals),
        P_Value = wt$p.value
      ))
    }
  }
  stats_results$FDR <- p.adjust(stats_results$P_Value, method = "BH")
  stats_results <- stats_results[order(stats_results$P_Value), ]
  write.csv(stats_results, file.path(out_res_dir, "goi_dependency_stats_latest.csv"), row.names = FALSE)
  print(head(stats_results))
  
  # 3. Visualization: Heatmap for Skin-lineage dependency
  skin_only <- goi_data %>% 
    filter(lineage_group == "Skin") %>%
    column_to_rownames("ModelID") %>%
    select(all_of(valid_genes)) %>%
    as.matrix()
  
  pheatmap(
    skin_only,
    main = "Dependency of Genes of Interest in Skin Lineage Cell Lines",
    show_colnames = TRUE,
    show_rownames = FALSE,
    color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
    filename = file.path(out_img_dir, "genes_of_interest_skin_dependency.png"),
    width = 10, height = 7
  )
  
  # 4. Visualization: Faceted Violin Plot for ALL Genes of Interest
  # Order genes by P-Value for the plot
  ordered_goi <- stats_results$Gene
  goi_data_long <- goi_data %>%
    pivot_longer(cols = all_of(valid_genes), names_to = "Gene", values_to = "Dependency") %>%
    mutate(Gene = factor(Gene, levels = ordered_goi))
  
  # Add p-value labels to the plot data
  p_labels <- stats_results %>%
    mutate(p_label = ifelse(FDR < 0.001, "FDR < 0.001", paste0("FDR = ", round(FDR, 3)))) %>%
    select(Gene, p_label)
  
  # Create the faceted plot
  p_violin <- ggplot(goi_data_long, aes(x = lineage_group, y = Dependency, fill = lineage_group)) +
    geom_violin(trim = FALSE, alpha = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.2, size = 0.5) +
    facet_wrap(~Gene, scales = "free_y") +
    geom_text(data = p_labels, aes(x = 1.5, y = 1.0, label = p_label), 
              inherit.aes = FALSE, size = 3, fontface = "italic") +
    labs(title = "Dependency Comparison: Skin vs Others",
         subtitle = "Genes ordered by significance (Wilcoxon Test)",
         x = "Lineage Group",
         y = "Dependency Score (CERES)") +
    theme_bw() +
    theme(strip.background = element_rect(fill = "gray90"),
          strip.text = element_text(face = "bold"),
          legend.position = "none") +
    scale_fill_manual(values = c("Skin" = "#d95f02", "Other" = "#7570b3"))

  ggsave(file.path(out_img_dir, "all_goi_dependency_violin.png"), plot = p_violin, width = 14, height = 10)
  
  message("Genes of Interest Inquiry completed.")
}

# =============================================================================
# Main Execution
# =============================================================================

main <- function() {
  # Paths
  crispr_path <- "./Data/depmap/CRISPRGeneEffect_25Q2.csv"
  model_path <- "./Data/depmap/Model.csv"
  essential_path <- "./Data/depmap/common_essentials.csv"
  out_img_dir <- "./Images/DepMap"
  out_res_dir <- "./Data/Results/DepMap"
  dir.create(out_img_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_res_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Genes of Interest
  genes_of_interest <- unique(c(
    "BRAF", "RET", "NTRK1", "NTRK2", "NTRK3",
    "KIT", "MTAP", "MAP2K1", "NRG1", "NRAS",
    "ERBB2", "TP53", "MDM2", "MTOR", "CCNE1", "CDKN2A",
    "KRAS", "NF1", "FGFR1", "FGFR2", "MET",
    "PIK3CA", "FBXW7", "CDK12", "ARID1A", "PPP2R1A", "FGFR3", "PTEN"
  ))
  
  # 1. Load and Clean
  depmap_obj <- load_and_annotate_depmap(crispr_path, model_path, essential_path)
  
  # 2. Run Module: Common Essentials
  spec_results <- run_essential_genes_module(depmap_obj, out_img_dir)
  
  # 3. Run Module: Genes of Interest Inquiry
  run_genes_of_interest_module(depmap_obj, genes_of_interest, out_img_dir, out_res_dir)
  
  message("\nPipeline completed successfully.")
}

# Run execution
main()
