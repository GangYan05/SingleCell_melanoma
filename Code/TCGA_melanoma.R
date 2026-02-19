# --- Script for Generating an OncoPrint from TCGA Melanoma Data ---

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(grid)

# =============================================================================
# Module 1: Data Loading and Processing
# =============================================================================

load_and_process_data <- function(clinical_pat_path, clinical_samp_path, cna_path, mut_path, rna_path, rppa_path) {
  message("Loading and processing data...")
  
  # 1. Clinical Data
  clinical_meta <- read.table(clinical_pat_path, header = TRUE, stringsAsFactors = FALSE, sep = "\t", quote = "")
  sample_meta <- read.table(clinical_samp_path, header = TRUE, stringsAsFactors = FALSE, sep = "\t", quote = "")
  
  os_data <- clinical_meta %>% select(Patient_ID = PATIENT_ID, OS_Months = OS_MONTHS) %>% distinct()
  tmb_data <- sample_meta %>% select(Sample_ID = SAMPLE_ID, TMB = TMB_NONSYNONYMOUS) %>% distinct()
  sample_to_patient <- sample_meta %>% select(Sample_ID = SAMPLE_ID, Patient_ID = PATIENT_ID) %>% distinct()
  
  # 2. CNA Data
  cna_raw <- read.table(cna_path, header = TRUE, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)
  cna_data <- cna_raw %>%
    select(-any_of(c("Entrez_Gene_Id", "Cytoband"))) %>%
    group_by(Hugo_Symbol) %>%
    summarise(across(everything(), ~ .[which.max(abs(.))]), .groups = "drop") %>%
    column_to_rownames("Hugo_Symbol")
  
  # 3. Mutation Data
  mut_raw <- read.table(mut_path, header = TRUE, stringsAsFactors = FALSE, sep = "\t", quote = "", fill = TRUE)
  mut_matrix <- mut_raw %>%
    filter(Variant_Classification != "Silent") %>%
    select(Hugo_Symbol, Tumor_Sample_Barcode) %>%
    distinct() %>%
    mutate(alteration = "MUT") %>%
    pivot_wider(names_from = Tumor_Sample_Barcode, values_from = alteration, values_fill = "") %>%
    as.data.frame()
  rownames(mut_matrix) <- mut_matrix$Hugo_Symbol
  mut_matrix$Hugo_Symbol <- NULL
  
  # 4. RNA Data
  rna_raw <- read.table(rna_path, header = TRUE, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)
  rna_data <- rna_raw %>%
    select(-any_of("Entrez_Gene_Id")) %>%
    group_by(Hugo_Symbol) %>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    column_to_rownames("Hugo_Symbol")
  
  # 5. Protein Data (RPPA)
  rppa_raw <- read.table(rppa_path, header = TRUE, stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)
  colnames(rppa_raw)[1] <- "Protein_ID"
  protein_data <- rppa_raw %>%
    mutate(Hugo_Symbol = gsub("\\|.*", "", Protein_ID)) %>%
    separate_rows(Hugo_Symbol, sep = " ") %>%
    mutate(Hugo_Symbol = gsub("_p[sStTyY].*|_phospho|_.*", "", Hugo_Symbol)) %>%
    mutate(Hugo_Symbol = case_when(
      Hugo_Symbol == "4E-BP1" ~ "EIF4EBP1",
      Hugo_Symbol == "ACC" ~ "ACACA",
      Hugo_Symbol == "AMPKalpha" ~ "PRKAA1",
      Hugo_Symbol == "c-Kit" ~ "KIT",
      Hugo_Symbol == "c-Met" ~ "MET",
      Hugo_Symbol == "eIF4G" ~ "EIF4G1",
      Hugo_Symbol == "GSK-3alpha-beta" ~ "GSK3B",
      Hugo_Symbol == "IGF-IR" ~ "IGF1R",
      Hugo_Symbol == "Lck" ~ "LCK",
      Hugo_Symbol == "MEK1" ~ "MAP2K1",
      Hugo_Symbol == "mTOR" ~ "MTOR",
      Hugo_Symbol == "p70S6K" ~ "RPS6KB1",
      Hugo_Symbol == "PDK1" ~ "PDPK1",
      Hugo_Symbol == "S6" ~ "RPS6",
      TRUE ~ Hugo_Symbol
    )) %>%
    group_by(Hugo_Symbol) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    column_to_rownames("Hugo_Symbol")
  
  # 6. M-Stage Simplification
  m_stage_simplified <- clinical_meta %>% 
    select(PATIENT_ID, PATH_M_STAGE, PATH_T_STAGE) %>%
    mutate(M_Stage_Simplified = case_when(
      grepl("m1", PATH_M_STAGE, ignore.case = TRUE) ~ "M1",
      grepl("m0", PATH_M_STAGE, ignore.case = TRUE) ~ "M0",
      TRUE ~ "Unknown"
    ))
  
  list(
    cna = cna_data,
    mut = mut_matrix,
    rna = rna_data,
    protein = protein_data,
    tmb = tmb_data,
    os = os_data,
    sample_to_patient = sample_to_patient,
    m_stage = m_stage_simplified
  )
}

# =============================================================================
# Module 2: Helper Functions
# =============================================================================

get_genes_of_interest <- function() {
  unique(c(
    "BRAF", "RET", "NTRK1", "NTRK2", "NTRK3",
    "RET", "KIT", "MTAP", "MAP2K1", "NRG1", "NRAS",
    "ERBB2", "TP53", "MDM2", "MTOR", "CCNE1", "CDKN2A",
    "MDM2", "KRAS", "NF1", "FGFR1", "FGFR2", "MET",
    "PIK3CA", "FBXW7", "CDK12", "ARID1A", "PPP2R1A", "FGFR3", "PTEN"
  ))
}

prepare_matrices_for_modules <- function(data_list, genes) {
  # Find common samples across Mutation, CNA, and RNA
  common_samples <- Reduce(intersect, list(
    colnames(data_list$mut),
    colnames(data_list$cna),
    colnames(data_list$rna)
  ))
  
  # Ensure genes exist in data
  valid_genes <- intersect(genes, rownames(data_list$rna))
  
  # Prepare aligned Clinical Data
  clinical_aligned <- data_list$sample_to_patient %>%
    filter(Sample_ID %in% common_samples) %>%
    left_join(data_list$m_stage, by = c("Patient_ID" = "PATIENT_ID")) %>%
    mutate(
      M_Stage = factor(M_Stage_Simplified, levels = c("M0", "M1", "Unknown")),
      T_Stage = ifelse(is.na(PATH_T_STAGE) | PATH_T_STAGE == "", "Unknown", PATH_T_STAGE)
    ) %>%
    column_to_rownames("Sample_ID")
  
  # Sort by M-Stage
  clinical_aligned <- clinical_aligned[order(clinical_aligned$M_Stage), ]
  sorted_samples <- rownames(clinical_aligned)
  
  list(
    mut = data_list$mut[valid_genes, sorted_samples, drop = FALSE],
    cna = data_list$cna[valid_genes, sorted_samples, drop = FALSE],
    rna = data_list$rna[valid_genes, sorted_samples, drop = FALSE],
    clinical = clinical_aligned,
    samples = sorted_samples,
    genes = valid_genes
  )
}

# =============================================================================
# Part 1: Genomic OncoPrint
# =============================================================================

run_oncoprint_module <- function(matrices, data_list, output_path) {
  message("Generating OncoPrint...")
  
  # 1. Prepare OncoPrint Matrix
  op_matrix <- matrix("", nrow = nrow(matrices$mut), ncol = ncol(matrices$mut),
                      dimnames = dimnames(matrices$mut))
  
  for (g in rownames(op_matrix)) {
    for (s in colnames(op_matrix)) {
      m <- matrices$mut[g, s]
      c <- matrices$cna[g, s]
      # Recode CNA numeric to string if needed (assuming input is already recoded or raw)
      # In load_and_process_data, we loaded raw. We need to recode here or in load.
      # Let's handle recoding quickly here for safety or assume 'process_data' did it.
      # Wait, valid 'cna' in previous script was recoded. Let's add that logic.
      
      c_str <- case_when(
        c == 2 ~ "AMP",
        c == -2 ~ "HOMDEL",
        TRUE ~ ""
      )
      
      # Combined
      vals <- c(m, c_str)
      op_matrix[g, s] <- paste(vals[vals != ""], collapse = ";")
    }
  }
  
  # 2. Annotations
  anno_df <- data.frame(Sample_ID = matrices$samples) %>%
    left_join(data_list$tmb, by = "Sample_ID") %>%
    left_join(data_list$sample_to_patient, by = "Sample_ID") %>%
    left_join(data_list$os, by = "Patient_ID")
  
  ha_bottom <- HeatmapAnnotation(
    TMB = anno_points(anno_df$TMB, pch = 16, size = unit(1, "mm"), 
                      gp = gpar(col = "black"), axis_param = list(side = "right"),
                      ylim = c(0, 160)),
    OS_Months = anno_points(anno_df$OS_Months, pch = 16, size = unit(1, "mm"), 
                            gp = gpar(col = "black"), axis_param = list(side = "right")),
    annotation_name_side = "left"
  )
  
  # 3. Plot
  col <- c("MUT" = "#008000", "AMP" = "red", "HOMDEL" = "blue")
  alter_fun <- list(
    background = function(x, y, w, h) grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA)),
    AMP = function(x, y, w, h) grid.rect(x, y, w - unit(1, "mm"), h - unit(1, "mm"), gp = gpar(fill = col["AMP"], col = NA)),
    HOMDEL = function(x, y, w, h) grid.rect(x, y, w - unit(1, "mm"), h - unit(1, "mm"), gp = gpar(fill = col["HOMDEL"], col = NA)),
    MUT = function(x, y, w, h) grid.rect(x, y, w - unit(1, "mm"), h * 0.66, gp = gpar(fill = col["MUT"], col = NA))
  )
  
  op <- oncoPrint(op_matrix,
                  alter_fun = alter_fun,
                  alter_fun_is_vectorized = FALSE,
                  col = col,
                  bottom_annotation = ha_bottom,
                  row_names_gp = gpar(fontsize = 8),
                  pct_gp = gpar(fontsize = 8),
                  remove_empty_columns = FALSE,
                  get_type = function(x) strsplit(x, ";")[[1]],
                  heatmap_legend_param = list(title = "Alterations", at = c("AMP", "HOMDEL", "MUT"), labels = c("Amplification", "Deep Deletion", "Mutation"))
  )
  
  png(output_path, width = 12, height = 8, units = "in", res = 300)
  draw(op, heatmap_legend_side = "right", annotation_legend_side = "right", column_title = "TCGA Melanoma")
  dev.off()
  message(paste("OncoPrint saved to", output_path))
}

# =============================================================================
# Part 2: Co-occurrence Heatmap
# =============================================================================

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
  
  png(output_plot, width = 10, height = 11, units = "in", res = 300)
  draw(ht, heatmap_legend_side = "bottom", column_title = "TCGA Melanoma Mutation Co-occurrence")
  dev.off()
  message(paste("Co-occurrence heatmap saved to", output_plot))
}

# =============================================================================
# Part 3: RNA Expression Heatmap
# =============================================================================

run_rna_module <- function(matrices, out_img_dir, out_res_dir) {
  message("Generating Staged RNA Expression Heatmap and Statistical Analysis...")
  
  # 1. Prepare Data
  rna_subset <- matrices$rna
  clinical_subset <- matrices$clinical
  ordered_samples <- matrices$samples
  valid_genes <- matrices$genes
  
  # 2. Statistical Analysis: M0 vs M1
  message("Performing statistical comparison: M0 vs M1...")
  m0_samples <- rownames(clinical_subset[clinical_subset$M_Stage == "M0", ])
  m1_samples <- rownames(clinical_subset[clinical_subset$M_Stage == "M1", ])
  
  stats_results <- data.frame()
  for (g in valid_genes) {
    val_m0 <- as.numeric(rna_subset[g, m0_samples])
    val_m1 <- as.numeric(rna_subset[g, m1_samples])
    val_m0 <- val_m0[!is.na(val_m0)]; val_m1 <- val_m1[!is.na(val_m1)]
    
    if (length(val_m0) > 1 & length(val_m1) > 1) {
      wt <- wilcox.test(val_m0, val_m1)
      mean_m0 <- mean(val_m0); mean_m1 <- mean(val_m1)
      log2fc <- log2((mean_m1 + 1) / (mean_m0 + 1))
      stats_results <- rbind(stats_results, data.frame(Gene = g, Mean_M0 = mean_m0, Mean_M1 = mean_m1, Log2FC = log2fc, P_Value = wt$p.value))
    }
  }
  stats_results$FDR <- p.adjust(stats_results$P_Value, method = "fdr")
  stats_results <- stats_results[order(stats_results$P_Value), ]
  write.csv(stats_results, file.path(out_res_dir, "m0_vs_m1_expression_stats.csv"), row.names = FALSE)
  
  # 3. Visualization
  # Z-score and Color
  expr_z <- t(scale(t(rna_subset)))
  expr_z[is.na(expr_z)] <- 0
  col_fun_intense = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  
  # Annotations
  ha_stage <- HeatmapAnnotation(
    M_Stage = clinical_subset$M_Stage,
    T_Stage = clinical_subset$T_Stage,
    col = list(
      M_Stage = structure(circlize::rand_color(length(levels(clinical_subset$M_Stage))), names = levels(clinical_subset$M_Stage)),
      T_Stage = structure(circlize::rand_color(length(unique(clinical_subset$T_Stage))), names = unique(clinical_subset$T_Stage))
    )
  )
  
  ht <- Heatmap(expr_z,
                name = "mRNA Expr (Z-score)",
                col = col_fun_intense,
                top_annotation = ha_stage,
                show_column_names = FALSE,
                cluster_rows = TRUE,
                cluster_columns = FALSE,
                column_split = clinical_subset$M_Stage,
                column_title = "mRNA expression by metastasis status",
                row_names_gp = gpar(fontsize = 8)
  )
  
  output_path <- file.path(out_img_dir, "tcga_melanoma_rna_heatmap_staged.png")
  png(output_path, width = 12, height = 8, units = "in", res = 300)
  draw(ht, merge_legend = TRUE)
  dev.off()
  message(paste("RNA heatmap and stats completed and saved to", out_img_dir))
}

# =============================================================================
# Part 4: Protein Heatmap
# =============================================================================

run_protein_module <- function(protein_data, output_path) {
  message("Generating Protein Expression Heatmap...")
  
  # Z-score
  prot_z <- t(scale(t(protein_data)))
  prot_z[is.na(prot_z)] <- 0
  
  col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  
  ht <- Heatmap(prot_z,
                name = "Protein (Z-score)",
                col = col_fun,
                show_row_names = TRUE,
                show_column_names = FALSE,
                cluster_rows = TRUE,
                cluster_columns = TRUE,
                row_names_gp = gpar(fontsize = 4),
                column_title = "Protein Expression Landscape (RPPA)"
  )
  
  png(output_path, width = 12, height = 10, units = "in", res = 300)
  draw(ht)
  dev.off()
  message(paste("Protein heatmap saved to", output_path))
}

# =============================================================================
# Main Execution Block
# =============================================================================

main <- function(
  run_oncoprint = TRUE,
  run_co_occurrence = TRUE,
  run_rna = TRUE,
  run_protein = TRUE,
  data = NULL,
  matrices = NULL
) {
  # Paths
  base_dir <- "./Data/skcm_tcga_pan_can_atlas_2018"
  out_img_dir <- "./Images/TCGA"
  out_res_dir <- "./Data/Results/TCGA"
  dir.create(out_img_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_res_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 1. Load Data
  if (is.null(data)) {
    message("Loading and processing data (this may take a moment)...")
    data <- load_and_process_data(
      clinical_pat_path = file.path(base_dir, "data_clinical_patient.txt"),
      clinical_samp_path = file.path(base_dir, "data_clinical_sample.txt"),
      cna_path = file.path(base_dir, "data_cna.txt"),
      mut_path = file.path(base_dir, "data_mutations.txt"),
      rna_path = file.path(base_dir, "data_mrna_seq_v2_rsem.txt"),
      rppa_path = file.path(base_dir, "data_rppa.txt")
    )
  } else {
    message("Using provided 'data' object. Skipping loading.")
  }
  
  # 2. Define Genes
  genes <- get_genes_of_interest()
  
  # 3. Prepare Common Matrics
  if (is.null(matrices)) {
    message("Preparing common matrices...")
    matrices <- prepare_matrices_for_modules(data, genes)
  } else {
    message("Using provided 'matrices' object.")
  }
  
  # 4. Execute Modules
  if (run_oncoprint) {
    run_oncoprint_module(
      matrices = matrices, 
      data_list = data, 
      output_path = file.path(out_img_dir, "tcga_melanoma_oncoprint_annotated.png")
    )
  }
  
  if (run_co_occurrence) {
    run_co_occurrence_module(
      matrices = matrices,
      output_csv = file.path(out_res_dir, "mutation_co_occurrence_tcga_melanoma.csv"),
      output_plot = file.path(out_img_dir, "mutation_co_occurrence_heatmap_tcga_melanoma_filtered.png")
    )
  }
  
  if (run_rna) {
    run_rna_module(
      matrices = matrices,
      out_img_dir = out_img_dir,
      out_res_dir = out_res_dir
    )
  }
  
  if (run_protein) {
    run_protein_module(
      protein_data = data$protein,
      output_path = file.path(out_img_dir, "tcga_melanoma_protein_heatmap.png")
    )
  }
  
  message("Selected modules completed successfully.")
  
  # Return the data and matrices invisibly so they can be reused
  invisible(list(data = data, matrices = matrices))
}

# =============================================================================
# Execution Instructions (Run these lines in your R console)
# =============================================================================

if (interactive()) {
  # ---------------------------------------------------------------------------
  # OPTION 1: First Run (Loads all data - takes time)
  # ---------------------------------------------------------------------------
  # Run this line first to load data and save it into 'res'.
  # res <- main()
  
  # ---------------------------------------------------------------------------
  # OPTION 2: Reuse Data (Fast - only runs specific modules)
  # ---------------------------------------------------------------------------
  # Once 'res' is created, you can re-run specific parts without reloading data.
  # Example: Run ONLY Co-occurrence Analysis
  # main(
  #   data = res$data, 
  #   matrices = res$matrices, 
  #   run_oncoprint = FALSE, 
  #   run_co_occurrence = TRUE, 
  #   run_rna = FALSE, 
  #   run_protein = FALSE
  # )
  
  # Example: Run ONLY OncoPrint
  # main(data = res$data, matrices = res$matrices, run_oncoprint = TRUE, run_co_occurrence = FALSE, run_rna = FALSE, run_protein = FALSE)

} else {
  # Command Line Execution (Rscript): Runs everything by default
  main()
}


# Example: Only re-run the Co-occurrence analysis, reusing the loaded 'res' data
main(data = res$data, matrices = res$matrices, 
     run_oncoprint = FALSE, 
     run_co_occurrence = TRUE, 
     run_rna = FALSE, 
     run_protein = FALSE)
