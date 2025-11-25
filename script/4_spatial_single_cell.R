library(Seurat)  # Load Seurat for single-cell and spatial transcriptomics analysis

# Load and process spatial transcriptomic data for multiple samples

# Define the paths for count matrices and images
counts_path1 <- "LC-3-4/"  # Directory containing the filtered_feature_bc_matrix.h5
image_path1 <- "LC-3-4/spatial"  # Directory containing spatial images

# Load spatial image
image1 <- Read10X_Image(image.dir = image_path1)

# Load count matrix and create Seurat object
spatial_data1 <- Load10X_Spatial(
  data.dir = counts_path1,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  image = image1
)

# Repeat for the second sample
tcounts_path2 <- "LC-1-2/"
image_path2 <- "LC-1-2/spatial/"
image2 <- Read10X_Image(image.dir = image_path2)
spatial_data2 <- Load10X_Spatial(
  data.dir = counts_path2,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  image = image2
)

# Repeat for the third sample
counts_path3 <- "REC/"
image_path3 <- "REC/spatial/"
image3 <- Read10X_Image(image.dir = image_path3)
spatial_data3 <- Load10X_Spatial(
  data.dir = counts_path3,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  image = image3
)

# Normalize the data using SCTransform
spatial_data1 <- SCTransform(spatial_data1, assay = "Spatial", verbose = FALSE)
spatial_data2 <- SCTransform(spatial_data2, assay = "Spatial", verbose = FALSE)
spatial_data3 <- SCTransform(spatial_data3, assay = "Spatial", verbose = FALSE)

# Assign sample identities
spatial_data1@meta.data$orig.ident <- 'spatial_data1'
spatial_data2@meta.data$orig.ident <- 'spatial_data2'
spatial_data3@meta.data$orig.ident <- 'spatial_data3'

# Merge datasets into one object
spatial_data.merge <- merge(spatial_data1, list(spatial_data2, spatial_data3))

# Perform dimensionality reduction and clustering
DefaultAssay(spatial_data.merge) <- "SCT"
VariableFeatures(spatial_data.merge) <- c(
  VariableFeatures(spatial_data1), 
  VariableFeatures(spatial_data2), 
  VariableFeatures(spatial_data3)
)
spatial_data.merge <- RunPCA(spatial_data.merge, verbose = FALSE)
spatial_data.merge <- FindNeighbors(spatial_data.merge, dims = 1:30)
spatial_data.merge <- FindClusters(spatial_data.merge, verbose = FALSE, resolution = 0.1)
spatial_data.merge <- RunUMAP(spatial_data.merge, dims = 1:30)

# Generate cluster and spatial plots
DimPlot(spatial_data.merge, reduction = "umap", label = TRUE, group.by = c('ident', 'orig.ident'))
SpatialDimPlot(spatial_data.merge, group.by = 'seurat_clusters')

# Define cluster-specific colors
my_cols <- c(
  '0' = '#FFF9C4',   # Cluster 0
  '1' = '#A7FFEB',   # Cluster 1
  '2' = '#FF6B6B',   # Cluster 2
  '3' = '#8E24AA'    # Cluster 3
)

spatial_data.merge@meta.data$color_cluster <- factor(spatial_data.merge$seurat_clusters, levels = names(my_cols))

# Visualize spatial clustering
SpatialDimPlot(spatial_data.merge, group.by = "color_cluster", pt.size.factor = 1.6, cols = my_cols)
SpatialDimPlot(spatial_data.merge, label = TRUE, label.size = 3, group.by = 'orig.ident')

# Highlight specific clusters
SpatialDimPlot(spatial_data.merge, cells.highlight = CellsByIdentities(object = spatial_data.merge, idents = c(0, 1, 2, 3)), facet.highlight = TRUE, ncol = 3)

# Interactive visualization
SpatialDimPlot(spatial_data.merge, interactive = TRUE)
SpatialFeaturePlot(spatial_data.merge, features = "GPX1", interactive = TRUE)
LinkedDimPlot(spatial_data.merge, group.by = 'orig.ident')
dev.off()

# Identify marker genes for all clusters
spatial_data.merge <- PrepSCTFindMarkers(spatial_data.merge)
all_markers <- FindAllMarkers(spatial_data.merge, 
                              only.pos = TRUE,  # Only positive markers
                              min.pct = 0.25,   # Expressed in at least 25% of cells
                              logfc.threshold = 0.25)  # Minimum log2 fold-change threshold

# Save results
write.csv(all_markers, '/Users/zhusenlin/Desktop/Article/Spatial/summary/1_Spaceranger_result/custom_analysis/cluster_markers.csv')

# Select top marker genes for visualization
top_markers <- all_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

# Dot plot visualization of key markers
DotPlot(spatial_data.merge, features = c('GJA1', 'TGM3', 'TGM1', 'KRT6B', 'A2ML1', 'LYZ1', 'FTH1', 'MT1E.1', 'BOLA', 'COL1A2')) + RotatedAxis()

# Spatial visualization of key genes
SpatialFeaturePlot(spatial_data.merge, features = c('ND1'), pt.size.factor = 1.5)

# Heatmap of top marker genes
top_markers <- all_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(spatial_data.merge, features = top_markers$gene) + NoLegend()

# Filter out ribosomal genes for better visualization
filtered_genes <- top_markers$gene[!grepl("^RPL", top_markers$gene)]

# Dot plot of selected marker genes across clusters
DotPlot(spatial_data.merge, features = unique(c('DCN', 'LUM', 'COL1A2', 'ACTA2', 'MYH11', 'MYL9', 'KRT5', 'KRT15', 'MT1E.1', 'BOLA', 'ZC3H10', 'KRT6A', 'GJA1', 'TGM3', 'TGM1', 'KRT6B', 'A2ML1', 'LYZ1', 'FTH1', 'COL1A2', 'S100A8', 'S100A9', 'CD3E', 'CD3D', 'CD4', 'PECAM1', 'TOP2A', 'VWF', 'KDR')), group.by = 'seurat_clusters') + RotatedAxis()






# Set cluster identities for spatial_data.merge
Idents(spatial_data.merge) = 'seurat_clusters'

# Generate a DotPlot for selected genes across specific clusters
dotplot_data <- DotPlot(spatial_data.merge, 
                        features = c('COX1','COX2','COX3','TXN','GSR','GLRX','MGST1','PRDX1','ESRRA'), 
                        idents = c(0, 1, 2)) + RotatedAxis()

dotplot_data  # Display the DotPlot

# Convert DotPlot data into a data frame format for further customization
dotplot_df <- dotplot_data$data

# Create a custom ggplot visualization to show temporal changes in gene expression
ggplot(dotplot_df, aes(x = id, y = avg.exp, color = features.plot, group = features.plot)) +
  geom_line(size = 1) +                              # Add line plot to show trends
  geom_point(size = dotplot_df$pct.exp/20,          # Adjust point size based on percentage expression
             alpha = dotplot_df$pct.exp/20) +       # Adjust transparency accordingly
  labs(x = "Cluster Idents", y = "Average Expression", 
       title = "Temporal Change of Gene Expression Across Clusters") +
  theme_minimal() +
  theme(legend.title = element_blank()) +           # Remove legend title for cleaner visualization
  facet_wrap(facets = ~ features.plot, scales = 'free') +  # Create separate plots for each gene
  scale_color_viridis_d(option = "plasma")          # Apply Viridis plasma color scheme
