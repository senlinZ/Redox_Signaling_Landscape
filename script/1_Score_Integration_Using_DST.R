# Set the working directory (recommend using relative paths)
setwd('path/to/your/project')

set.seed(123)

# Read the data
matrix_data <- read.csv('Supplementary_table_8.csv', row.names = 1)

# Extract group information
group <- matrix_data[,1:6]
matrix_data <- matrix_data[,-c(1:6)]

# Load required R packages
library(stringr)
library(tools)
library(reshape2)
library(dplyr)
library(readxl)

# Load RDS file
matrix_data <- readRDS('all_pathway_enrichment_scores.rds')

# Process antioxidant gene names
antioxidant_names <- str_replace_all(str_remove(dir('antioxidant_goterm/'),'.txt'), ' ', '.')
antioxidant_names <- str_replace_all(antioxidant_names, '-', '.')
antioxidant_names <- toTitleCase(antioxidant_names)

# Process response gene names
response_names <- str_replace_all(str_remove(dir('response_goterm/'),'.txt'), ' ', '.')
response_names <- str_replace_all(response_names, '-', '.')
response_names <- toTitleCase(response_names)

# Process oxidase generation gene names
generation_names <- str_replace_all(str_remove(dir('generation_goterm/'),'.txt'), ' ', '.')
generation_names <- str_replace_all(generation_names, '-', '.')
generation_names <- toTitleCase(generation_names)

# Check data matching conditions
table(colnames(matrix_data) %in% antioxidant_names)
table(colnames(matrix_data) %in% response_names)
table(colnames(matrix_data) %in% generation_names)

# Split the dataset into three categories
antioxidant <- matrix_data[,2:20]
response <- matrix_data[,21:33]
generation <- matrix_data[,34:49]

# Normalization function: scale values to [0,1] and normalize by column sum
normalize_combined <- function(mat) {
  min_vals <- apply(mat, 2, min)
  max_vals <- apply(mat, 2, max)
  scaled_mat <- sweep(mat, 2, min_vals, FUN = "-") 
  scaled_mat <- sweep(scaled_mat, 2, max_vals - min_vals, FUN = "/") 
  
  # Add a small bias to avoid division by zero
  bias <- 1e-10
  col_sums <- colSums(scaled_mat) + bias
  normalized_mat <- sweep(scaled_mat, 2, col_sums, FUN = "/")
  
  return(normalized_mat)
}

# Apply normalization
normalized_matrix <- normalize_combined(antioxidant)
normalized_matrix <- normalize_combined(generation)
normalized_matrix <- normalize_combined(response)

# Define a fusion function to combine belief scores
fusion <- function(a, b) {
  m1 <- as.numeric(a)
  m2 <- as.numeric(b)
  
  K <- sum(m1 * m2) 
  res <- sum(m1 * m2) 
  K <- K - res
  combined <- (m1 * m2) / (1 - K)
  normalized_combined <- combined / sum(combined)
  
  return(normalized_combined)
}

# Compute the fused score across columns
combined_belief <- normalized_matrix[, 1]
for (i in 2:ncol(normalized_matrix)) {
  combined_belief <- fusion(combined_belief, normalized_matrix[, i])
}

# Store the fused scores in respective datasets
antioxidant$Antioxidant_fusion <- combined_belief
generation$Generation_fusion <- combined_belief
response$Response_fusion <- combined_belief

# Save the computed scores
all_score <- list(antioxidant = antioxidant, generation = generation, response = response)
saveRDS(all_score, 'all_scores.rds')

# Load system classification data
system <- as.data.frame(read_xlsx('../tissue_clusters.xlsx', sheet = 'All_summary_three_core_pathways'))
a1$tissue_cluster <- a1$tissue
a1 <- a1[, -1]
a1 <- merge(a1, system[, -1], by = 'tissue_cluster')

# Save final processed data
write.csv(a1, '../mean_of_DST_value.csv', row.names = TRUE)
