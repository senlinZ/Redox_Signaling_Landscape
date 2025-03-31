# Load the CSV file and set the first column as row names
data.use <- read.csv('AD_afterrun_palantir_file.csv', row.names = 1)

# Filter out specific cell types (i.e., removing unwanted samples from the dataset)
data.use = data.use[!data.use$orig.ident %in% c('Public_sc1','Public_sc2','CAU_Muscle_sc1','CAU_Muscle_sc2'),]

# Display the column names of the dataset
colnames(data.use)

# Load the dplyr library for data manipulation
library(dplyr)

# Arrange the data based on Pseudotime (sorting cells by their pseudotime values)
data.use <- arrange(data.use, Pseudotime)

# Load ggplot2 for plotting
library(ggplot2)

# Plot entropy changes along the pseudotime for different cell types (AD_NB)
ggplot(data.use, aes(Pseudotime, Entropy, color = AD_NB, group = AD_NB)) +
  geom_point(size = 0.5) +
  ggsci::scale_color_d3() +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme_classic(base_size = 15)

# Comment on the plot: "HSC" and "precursor" cell types show high entropy at the start of pseudotime, indicating rapid differentiation, and the rate of decrease is slow.

# Load the mgcv package for Generalized Additive Models (GAM)
library(mgcv)

# Fit a GAM model to predict entropy as a smooth function of pseudotime
model <- gam(Entropy ~ s(Pseudotime), data = data.use)

# Extract fitted values from the model and add them to the dataset
data.use$fitted <- model$fitted.values

# Compute the first derivative (change rate) of fitted values to find where the change is greatest
data.use$diff <- abs(c(NA, diff(data.use$fitted) / diff(data.use$Pseudotime)))

# Normalize the change rate to a 0-1 scale (Max-Min normalization)
data.use$diff <- (data.use$diff - min(data.use$diff, na.rm = T)) / (max(data.use$diff, na.rm = T) - min(data.use$diff, na.rm = T))

# Compute the second derivative (to find the extrema or peaks in the curve)
data.use$diff.diff <- abs(c(NA, diff(data.use$diff) / diff(data.use$Pseudotime)))

# Plot to visualize the entropy change along pseudotime, highlighting the regions of maximum rate of change
ggplot(data.use, aes(Pseudotime, Entropy, color = AD_NB, group = AD_NB)) +
  geom_point(size = 0.5, color = 'grey') +
  geom_vline(xintercept = c(0.25, 0.68)) +  # Vertical lines at points where the first and second derivatives are max
  geom_point(aes(Pseudotime, fitted), color = "#BE1E6F", size = 0.9) +
  geom_point(aes(Pseudotime, diff), color = "#104676", size = 0.9) +
  ggsci::scale_color_d3() +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme_classic(base_size = 15)

# Pathway scoring for antioxidant activity and related processes
library(Seurat)

# Load the Seurat object containing single-cell data
pbmc = readRDS('Rumen_epithelial.rds')

# Subset the data for specific samples
pbmc = subset(pbmc, subset = orig.ident %in% c('ZJU_AD_sc1', 'ZJU_AD_sc2', 'ZJU_AD_sc3', 'NWAFU_sc',
                                               'CAU_Papillary_sc2', 'CAU_Papillary_sc1'))

# Display the dimensions of the subsetted pbmc object
dim(pbmc)

# Define the GO terms related to antioxidant and metabolic pathways
path = c('GO:0009061', 'GO:0009060', 'GO:0006631', 'GO:0006119', 'GO:0006099', 'GO:0061621', 'GO:0016209',
         'GO:0006749', 'GO:0042182', 'GO:0004601', 'GO:0004459', 'GO:0004784')

# Corresponding pathway names
names = c('Anaerobic respiration', 'Aerobic respiration', 'Fatty acid metabolic process', 'Oxidative phosphorylation',
          'Tricarboxylic acid cycle', 'Canonical glycolysis', 'Antioxidant activity', 'Glutathione metabolic process',
          'Ketone catabolic process', 'Peroxidase activity', 'Lactate dehydrogenase activity', 'Superoxide dismutase activity')

# Connect to the Ensembl database
mart <- useMart("ensembl", dataset = "btaurus_gene_ensembl")

# Iterate over the defined GO terms to retrieve associated genes
for (i in 2:length(names)) {
  genes <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name", "go_id"),
    filters = "go",
    values = path[i],
    mart = mart
  )
  
  # Add module scores to the Seurat object based on the gene sets for each pathway
  pbmc = AddModuleScore(
    pbmc,
    features = list(unique(genes$external_gene_name)),
    name = names[i],
    assay = 'RNA'
  )
}

# Extract metadata from the Seurat object
metadata = pbmc@meta.data

# Combine the metadata with the data.use dataset for further analysis
data.use1 = cbind(data.use, metadata[, c(43:53)])
data.use1 = data.use1[, c(16, 17, 21:33)]

# Add specific gene expressions (e.g., PPARG, PPARD, ESRRA) to the dataset
data.use1$PPARG = pbmc@assays$RNA@data['PPARG',]
data.use1$PPARD = pbmc@assays$RNA@data['PPARD',]
data.use1$ESRRA = pbmc@assays$RNA@data['ESRRA',]

# Transform the dataset to long format for plotting
library(tidyr)
data_long <- data.use1 %>%
  pivot_longer(cols = c(`oxidative phosphorylation1`, `fatty acid metabolic process1`, `tricarboxylic acid cycle1`,
                        `canonical glycolysis1`, `Antioxidant activity1`, `glutathione metabolic process1`,
                        `ketone catabolic process1`, `peroxidase activity1`, `aerobic respiration1`,
                        `lactate dehydrogenase activity1`, `superoxide dismutase activity1`, PPARG, PPARD, ESRRA),
               names_to = "MetabolicProcess", values_to = "Expression")

# Plot a facet grid of metabolic process expressions along pseudotime
ggplot(data_long, aes(Pseudotime, Expression)) +
  geom_point(size = 0.5, color = "grey") +
  geom_point(aes(Pseudotime, diff / 3), size = 0.2, color = "#104676", inherit.aes = FALSE) +
  geom_smooth(color = "red") +
  geom_vline(xintercept = c(0.25, 0.6)) +
  facet_wrap(~ MetabolicProcess, scales = "free_y") +
  theme_classic(base_size = 15) +
  labs(y = "Expression Level", x = "Pseudotime")

# Create a separate plot for a specific metabolic process, such as canonical glycolysis
ggplot(data.use1, aes(Pseudotime, `canonical glycolysis1`)) +
  geom_point(size = 0.5, color = "grey") +
  geom_point(inherit.aes = FALSE, data = data.use1, aes(Pseudotime, diff / 3), size = 0.2, color = "red") +
  geom_smooth() +
  geom_vline(xintercept = c(0.069, 0.156)) +
  theme_classic(base_size = 15)

# Comment: After a sharp change, new metabolic patterns begin to emerge. 
# These plots, when combined with gene expression data, help describe the metabolic reprogramming.
