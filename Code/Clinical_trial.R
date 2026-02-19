# --- Step 1: Load Clinical Trial Data ---

# Set path to the clinical trials file
clinical_data_path <- "./Data/clinical_trials/melanoma clinical trials.xlsx"

# Check if file exists
if (!file.exists(clinical_data_path)) {
  stop("File not found: ", clinical_data_path)
}

# Load libraries
library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
library(pheatmap)
library(RColorBrewer)

message("Loading clinical trial data...")
melanoma_clinical_data <- read_excel(clinical_data_path)

# Ensure column names are clean (remove spaces/special char) for easy handling
colnames(melanoma_clinical_data) <- make.names(colnames(melanoma_clinical_data))
# Likely columns based on previous check: Trial.Title, Trial.Drug.s., Drug.Category

# --- Step 2: Define Drug Dictionary ---
# Common drugs used in melanoma treatment
drug_dict <- list(
  "Pembrolizumab" = c("pembrolizumab", "keytruda"),
  "Nivolumab" = c("nivolumab", "opdivo"),
  "Ipilimumab" = c("ipilimumab", "yervoy"),
  "Dabrafenib" = c("dabrafenib", "tafinlar"),
  "Trametinib" = c("trametinib", "mekinist"),
  "Vemurafenib" = c("vemurafenib", "zelboraf"),
  "Cobimetinib" = c("cobimetinib", "cotellic"),
  "Encorafenib" = c("encorafenib", "braftovi"),
  "Binimetinib" = c("binimetinib", "mektovi"),
  "Atezolizumab" = c("atezolizumab", "tecentriq"),
  "Relatlimab" = c("relatlimab"),
  "Temozolomide" = c("temozolomide"),
  "Dacarbazine" = c("dacarbazine"),
  "Interleukin-2" = c("interleukin-2", "il-2", "aldesleukin"),
  "TVEC" = c("talimogene laherparepvec", "t-vec", "imlygic")
)

# Function to extract standardized drug names from text
extract_drugs <- function(text, dict) {
  text <- tolower(text)
  found_drugs <- c()
  
  for (drug_std in names(dict)) {
    aliases <- dict[[drug_std]]
    # Check if any alias is present in the text
    if (any(str_detect(text, fixed(aliases)))) {
      found_drugs <- c(found_drugs, drug_std)
    }
  }
  
  if (length(found_drugs) == 0) return(NA)
  return(paste(sort(unique(found_drugs)), collapse = ";"))
}

# --- Step 3: Text Mining ---
message("Mining drug information...")

# Combine Title and Drug(s) column for text mining clarity if needed, 
# but rely mainly on Trial.Drug.s.
# Using 'Trial.Drug.s.' as primary source based on colnames.txt
drug_col <- grep("Drug", colnames(melanoma_clinical_data), value = TRUE)[1] 

processed_trials <- melanoma_clinical_data %>%
  rowwise() %>%
  mutate(
    Cleaned_Drugs = extract_drugs(!!sym(drug_col), drug_dict),
    Num_Drugs = ifelse(is.na(Cleaned_Drugs), 0, str_count(Cleaned_Drugs, ";") + 1),
    Therapy_Type = case_when(
      Num_Drugs == 0 ~ "Unknown/Other",
      Num_Drugs == 1 ~ "Monotherapy",
      Num_Drugs > 1 ~ "Combination"
    )
  ) %>%
  ungroup()

# Filter out trials with no identified main drugs for plotting focus
analyzed_trials <- processed_trials %>% 
  ungroup() %>%
  as.data.frame() %>% 
  filter(Num_Drugs > 0)

# --- Step 4: Visualizations ---

# A. Top Frequent Drugs (Bar Plot)
drug_counts <- analyzed_trials %>%
  separate_rows(Cleaned_Drugs, sep = ";") %>%
  dplyr::count(Cleaned_Drugs, sort = TRUE) %>%
  top_n(15, n)

p1 <- ggplot(drug_counts, aes(x = reorder(Cleaned_Drugs, n), y = n, fill = Cleaned_Drugs)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 15 Most Frequent Medical Interventions in Melanoma Trials",
       x = "Drug Name", y = "Number of Trials") +
  theme(legend.position = "none")

ggsave("./Images/clinical_trials_top_drugs_bar.png", p1, width = 8, height = 6)
message("Saved top drugs bar plot.")


# B. Monotherapy vs Combination (Pie Chart)
therapy_counts <- analyzed_trials %>%
  dplyr::count(Therapy_Type) %>%
  mutate(prop = n / sum(n) * 100,
         label = paste0(Therapy_Type, "\n", round(prop, 1), "%"))

p2 <- ggplot(therapy_counts, aes(x = "", y = prop, fill = Therapy_Type)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0) +
  theme_void() +
  theme(legend.position = "right") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5)) +
  labs(title = "Trial Distribution: Monotherapy vs Combination", fill = "Therapy Type") +
  scale_fill_brewer(palette = "Set3")

ggsave("./Images/clinical_trials_type_pie.png", p2, width = 7, height = 5)
message("Saved therapy type pie chart.")


# C. Specific Combinations (Heatmap of top drugs)
# Create a co-occurrence matrix for top 10 drugs
top_drugs <- drug_counts$Cleaned_Drugs[1:10]
n_top <- length(top_drugs)
combo_mat <- matrix(0, nrow = n_top, ncol = n_top, dimnames = list(top_drugs, top_drugs))

for (drugs_str in analyzed_trials$Cleaned_Drugs) {
  drugs <- unlist(str_split(drugs_str, ";"))
  # Only consider drugs in our top list
  drugs <- intersect(drugs, top_drugs)
  
  if (length(drugs) >= 2) {
    for (i in 1:(length(drugs)-1)) {
      for (j in (i+1):length(drugs)) {
        d1 <- drugs[i]
        d2 <- drugs[j]
        combo_mat[d1, d2] <- combo_mat[d1, d2] + 1
        combo_mat[d2, d1] <- combo_mat[d2, d1] + 1
      }
    }
  }
}

# Plot Heatmap (Lower Triangle Only)
library(corrplot)

# Re-fill the matrix for corrplot (it handles lower.tri itself)
combo_mat <- matrix(0, nrow = n_top, ncol = n_top, dimnames = list(top_drugs, top_drugs))
for (drugs_str in analyzed_trials$Cleaned_Drugs) {
  drugs <- unlist(str_split(drugs_str, ";"))
  drugs <- intersect(drugs, top_drugs)
  if (length(drugs) >= 2) {
    for (i in 1:(length(drugs)-1)) {
      for (j in (i+1):length(drugs)) {
        d1 <- drugs[i]; d2 <- drugs[j]
        combo_mat[d1, d2] <- combo_mat[d1, d2] + 1
        combo_mat[d2, d1] <- combo_mat[d2, d1] + 1
      }
    }
  }
}

png("./Images/clinical_trials_drug_combinations_corrplot.png", width = 8, height = 7, units = "in", res = 300)
# Use 'is.corr = FALSE' since we are plotting counts, not correlations
corrplot(combo_mat, 
         method = "color", 
         type = "lower", 
         diag = FALSE,
         is.corr = FALSE,
         col = colorRampPalette(c("white", "blue"))(200),
         addCoef.col = "black", # Add counts as text
         tl.col = "black",      # Text label color
         tl.srt = 90,           # Text label rotation
         cl.pos = "r",          # Color legend on right
         title = "Top Drug Co-occurrence in Combination Trials",
         mar = c(0,0,2,0))      # Margin for title
dev.off()
message("Saved combination corrplot.")

message("Clinical trial analysis complete.")

