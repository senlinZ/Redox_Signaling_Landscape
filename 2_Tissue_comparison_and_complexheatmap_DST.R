### Load required libraries
library(reshape2)
library(dplyr)
library(ComplexHeatmap)
library(RColorBrewer)

### Load and process the matrix data
matrix_data = readRDS('matrix_data.rds')
matrix_data_melt = melt(matrix_data)

# Summarize the average expression by tissue and variable
all_oxygen_signal_tissue_cluster <- matrix_data_melt %>% 
  group_by(tissue, variable) %>% 
  summarise(sum_avg_exp_scaled = mean(value)) %>% 
  arrange(desc(sum_avg_exp_scaled))

# Convert from long to wide format
all_oxygen_signal_tissue_cluster = dcast(all_oxygen_signal_tissue_cluster, tissue ~ variable)
rownames(all_oxygen_signal_tissue_cluster) = all_oxygen_signal_tissue_cluster[,1]
all_oxygen_signal_tissue_cluster = all_oxygen_signal_tissue_cluster[,-1]

# Remove unnecessary columns
all_oxygen_signal_tissue_cluster = all_oxygen_signal_tissue_cluster[,-c(49:72)]

# Load tissue annotation data
all_oxygen_signal_tissue_cluster_annotation = readRDS('all_oxygen_signal_tissue_cluster_annotation.rds')

# Select digestive tissues
digestive = all_oxygen_signal_tissue_cluster_annotation$tissue %in% c(
  'Rumen','Reticulum','Omasum','Abomasum','Duodenum','Jejunum','Ileum','Cecum','Colon','Rectum'
)

# Subset data for digestive tissues
use = all_oxygen_signal_tissue_cluster[digestive,]
use_annotation = all_oxygen_signal_tissue_cluster_annotation[digestive,]

# Rename some columns for clarity
colnames(use)[20:32] = c('re7','re8','re9','re10','re11','re12','re13','re1','re3','re4','re5','re6','re2')

# Get unique tissues and generate color mapping
unique_tissues <- unique(use_annotation$tissue)
tissue_colors <- setNames(colorRampPalette(brewer.pal(9, "Set1"))(length(unique_tissues)), unique_tissues)
tissue_color_vector <- tissue_colors[use_annotation$tissue]

tissue_annotation <- rowAnnotation(tissue = anno_simple(use_annotation$tissue, col = tissue_colors))

# Compute column means for barplot annotation
col_sums <- colMeans(use)

top_annotation <- HeatmapAnnotation(
  barplot = anno_barplot(col_sums, gp = gpar(fill = "navy"), border = TRUE),
  annotation_height = unit(2, "cm")
)

tissue_factor <- factor(use_annotation$tissue)

### Separate immune and epithelial cells
table(use_annotation$`Cell type`)

immune_cells = use_annotation$`Cell type` %in% c(
  'B cells','CD4+ T cells','Conventional dendritic cells type 1','Conventional dendritic cells type 2',
  'Cytotoxic CD8+ T cells','Macrophages','Mast cells','Neutrophils','NK T cells','Other immune cells',
  'Plasma cells','Proliferative CD4+ T cells','Proliferative T cells','T cells'
)

epithelial_cells = use_annotation$`Cell type` %in% c(
  'Basal cells','Chief cells','Enterocytes','Epithelial stem cells','Goblet cells','Intestinal stem cells',
  'Isthmus cells','Keratinocytes','Mesenchymal stem cells','Other epithelial cells','Mucous neck cells',
  'Parietal cells','Pit cells','Progenitor cells','Proliferative basal cells','Spinous cells',''
)

# Create heatmaps for epithelial cells
p1 <- Heatmap(use[epithelial_cells,1:19], name = "mat1",
              col = colorRampPalette(c("navy", "white", "firebrick3"))(50),
              left_annotation = tissue_annotation,
              split = tissue_factor[epithelial_cells],
              cluster_rows = TRUE,
              show_row_names = TRUE,
              show_column_names = TRUE)

p2 <- Heatmap(use[epithelial_cells,20:32], name = "mat1",
              col = colorRampPalette(c("navy", "white", "firebrick3"))(50),
              split = tissue_factor[epithelial_cells],
              cluster_rows = TRUE,
              show_row_names = TRUE,
              show_column_names = TRUE)

p3 <- Heatmap(use[epithelial_cells,33:48], name = "mat1",
              col = colorRampPalette(c("navy", "white", "firebrick3"))(50),
              split = tissue_factor[epithelial_cells],
              cluster_rows = TRUE,
              show_row_names = TRUE,
              show_column_names = TRUE)

# Row annotation for antioxidant responses
p4 = rowAnnotation(Anti = anno_barplot(use_annotation$Antioxidant_fusion[epithelial_cells]))
p5 = rowAnnotation(Res = anno_barplot(use_annotation$Response_fusion[epithelial_cells]))
p6 = rowAnnotation(Gen = anno_barplot(use_annotation$Generation_fusion[epithelial_cells]))

dev.off()
# Combine all heatmaps
p4 + p5 + p6 + p1 + p2 + p3

# Repeat the process for immune cells
tissue_color_vector <- tissue_colors[use_annotation$tissue[immune_cells]]
tissue_annotation <- rowAnnotation(tissue = anno_simple(use_annotation$tissue[immune_cells], col = tissue_colors))

p1 <- Heatmap(use[immune_cells,1:19], name = "mat1",
              col = colorRampPalette(c("navy", "white", "firebrick3"))(50),
              left_annotation = tissue_annotation,
              split = tissue_factor[immune_cells],
              cluster_rows = TRUE,
              show_row_names = TRUE,
              show_column_names = TRUE)

p2 <- Heatmap(use[immune_cells,20:32], name = "mat1",
              col = colorRampPalette(c("navy", "white", "firebrick3"))(50),
              split = tissue_factor[immune_cells],
              cluster_rows = TRUE,
              show_row_names = TRUE,
              show_column_names = TRUE)

p3 <- Heatmap(use[immune_cells,33:48], name = "mat1",
              col = colorRampPalette(c("navy", "white", "firebrick3"))(50),
              split = tissue_factor[immune_cells],
              cluster_rows = TRUE,
              show_row_names = TRUE,
              show_column_names = TRUE)

p4 = rowAnnotation(Anti = anno_barplot(use_annotation$Antioxidant_fusion[immune_cells]))
p5 = rowAnnotation(Res = anno_barplot(use_annotation$Response_fusion[immune_cells]))
p6 = rowAnnotation(Gen = anno_barplot(use_annotation$Generation_fusion[immune_cells]))

dev.off()
p4 + p5 + p6 + p1 + p2 + p3














# Comparison
# Load required packages
library(dplyr)
library(limma)
library(tidyr)
library(tibble)
library(stringr)

# Load dataset containing differential expression analysis data for the gastrointestinal tract
matrix_data_melt = readRDS('/Users/senzhu/Documents/ZJU_Project/Single_cell_antioxidant/Manuscript Submission NC 2024.4.11/New FIgure/Fig1/Figure supplementary table 1/胃肠道差异分析数据.rds')

# Define tissue categories
forestomach = c(
  unique(matrix_data_melt[str_detect(matrix_data_melt$tissue, 'Rumen'), "tissue"]),
  unique(matrix_data_melt[str_detect(matrix_data_melt$tissue, 'Reticulum'), "tissue"]),
  unique(matrix_data_melt[str_detect(matrix_data_melt$tissue, 'Omasum'), "tissue"])
)

Abomasum = c(
  unique(matrix_data_melt[str_detect(matrix_data_melt$tissue, 'Abomasum'), "tissue"])
)

Intestine = c(
  unique(matrix_data_melt[str_detect(matrix_data_melt$tissue, 'Duodenum'), "tissue"]),
  unique(matrix_data_melt[str_detect(matrix_data_melt$tissue, 'Jejunum'), "tissue"]),
  unique(matrix_data_melt[str_detect(matrix_data_melt$tissue, 'Ileum'), "tissue"]),
  unique(matrix_data_melt[str_detect(matrix_data_melt$tissue, 'Cecum'), "tissue"]),
  unique(matrix_data_melt[str_detect(matrix_data_melt$tissue, 'Colon'), "tissue"]),
  unique(matrix_data_melt[str_detect(matrix_data_melt$tissue, 'Rectum'), "tissue"])
)

# Define epithelial and immune cell groups
forestomach_epi = c('Rumen_0', 'Rumen_6', 'Rumen_7', 'Rumen_10', 'Rumen_13', 'Rumen_15',
                    'Reticulum_2', 'Reticulum_8', 'Reticulum_16', 'Reticulum_17',
                    'Omasum_3', 'Omasum_5', 'Omasum_6', 'Omasum_14')

Other_epi = c('Abomasum_0', 'Abomasum_1', 'Abomasum_2', 'Abomasum_3', 'Abomasum_4', 'Abomasum_5',
              'Abomasum_6', 'Abomasum_12', 'Duodenum_3', 'Duodenum_11', 'Duodenum_12',
              'Duodenum_13', 'Duodenum_17', 'Duodenum_19', 'Duodenum_22', 'Jejunum_3',
              'Jejunum_11', 'Jejunum_16', 'Jejunum_18', 'Jejunum_27', 'Jejunum_28',
              'Ileum_8', 'Ileum_11', 'Ileum_20', 'Ileum_21', 'Ileum_23', 'Colon_9',
              'Colon_10', 'Rectum_4', 'Rectum_8', 'Rectum_16', 'Rectum_22')

forestomach_immune = c('Rumen_18', 'Rumen_16', 'Rumen_8', 'Rumen_1',
                       'Reticulum_13', 'Reticulum_15', 'Reticulum_19', 'Reticulum_18',
                       'Reticulum_0', 'Reticulum_3', 'Reticulum_12', 'Omasum_2', 'Omasum_9',
                       'Omasum_15', 'Omasum_12', 'Omasum_13', 'Omasum_16')

Other_immune = c('Abomasum_9', 'Abomasum_8', 'Abomasum_10', 'Abomasum_15', 'Abomasum_17',
                 'Duodenum_10', 'Duodenum_5', 'Duodenum_15', 'Duodenum_18', 'Duodenum_0',
                 'Duodenum_6', 'Duodenum_1', 'Duodenum_4', 'Duodenum_2', 'Duodenum_9',
                 'Duodenum_8', 'Jejunum_0', 'Jejunum_1', 'Jejunum_2', 'Jejunum_5', 'Jejunum_6',
                 'Jejunum_7', 'Jejunum_8', 'Jejunum_10', 'Jejunum_14', 'Jejunum_15', 'Jejunum_19',
                 'Jejunum_20', 'Jejunum_21', 'Jejunum_22', 'Jejunum_26', 'Jejunum_29',
                 'Ileum_28', 'Ileum_5', 'Ileum_7', 'Ileum_14', 'Ileum_26', 'Ileum_18',
                 'Ileum_22', 'Ileum_19', 'Ileum_9', 'Ileum_17', 'Ileum_16', 'Ileum_3', 'Ileum_1',
                 'Ileum_0', 'Ileum_2', 'Colon_18', 'Colon_8', 'Colon_7', 'Colon_12', 'Colon_16',
                 'Colon_6', 'Colon_4', 'Colon_5', 'Colon_1', 'Colon_0', 'Colon_2', 'Colon_3',
                 'Colon_14', 'Rectum_12', 'Rectum_15', 'Rectum_7', 'Rectum_10', 'Rectum_13',
                 'Rectum_24', 'Rectum_14', 'Rectum_3', 'Rectum_18', 'Rectum_2', 'Rectum_0', 'Rectum_19')

# Define groups for comparison
group1 <- forestomach_epi  # Group 1
group2 <- Other_epi  # Group 2

# Get all unique variables
all_variables <- unique(matrix_data_melt$variable)

# Initialize results storage
results_list <- list()

# Loop through all variables
for (selected_variable in all_variables) {
  
  # Filter data for the selected variable within group1 and group2
  filtered_data <- matrix_data_melt %>%
    filter(variable == selected_variable & tissue %in% c(group1, group2))
  
  # Extract values for both groups
  group1_values <- filtered_data %>% filter(tissue %in% group1) %>% pull(value)
  group2_values <- filtered_data %>% filter(tissue %in% group2) %>% pull(value)
  
  # Calculate logFC (log fold change), adding a pseudocount to avoid division by zero
  logFC <- log2((mean(group1_values, na.rm = TRUE) + 1) / (mean(group2_values, na.rm = TRUE) + 1))
  
  # Perform t-test to calculate p-value
  if (length(group1_values) > 1 & length(group2_values) > 1) {
    t_test_result <- t.test(group1_values, group2_values, var.equal = TRUE)
    p_value <- t_test_result$p.value
  } else {
    p_value <- NA  # Avoid errors due to insufficient data
  }
  
  # Store results
  results_list[[selected_variable]] <- data.frame(
    Variable = selected_variable,
    logFC = logFC,
    P_value = p_value
  )
}

# Combine all results into a single dataframe
results_df <- do.call(rbind, results_list)

# Print results
print(results_df)
