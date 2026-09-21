############################################################
# 1. Load required package
############################################################

# Load Seurat, the R package used for single-cell RNA-seq analysis.
# Seurat provides functions for working with the Seurat object,
# metadata, gene expression, dimensionality reductions and plots.
library(Seurat)


############################################################
# 2. Set the working directory
############################################################

# Set the working directory to the folder containing the RData file.
# This means R will look inside this folder when we load the dataset.
setwd("/Users/mdreazmorshed/Desktop/Thesis-project/data")

# Print the current working directory to confirm that it is correct.
getwd()

# List the files in the working directory.
# This confirms that the RData and module reference files are present.
list.files()


############################################################
# 3. Load the Seurat object
############################################################

# Load the RData file into the current R session.
# This file should contain a Seurat object called "obj_small".
load("ICI_cRA_joint_small_for_module_pipeline.RData")

# Show all R objects currently present in the environment.
# After loading, we expect to see "obj_small".
ls()


############################################################
# 4. Inspect the basic structure of obj_small
############################################################

# Check the class of the object.
# We expect this to return "Seurat".
class(obj_small)

# Print a summary of the Seurat object.
# This gives information about the number of features (genes),
# cells, assays and dimensionality reductions.
obj_small

# Return the dimensions of the expression dataset.
# For a Seurat object:
# first number  = number of genes/features
# second number = number of cells
dim(obj_small)

# Print the number of genes and cells with clear labels.
cat(
  "Genes:", nrow(obj_small),
  "\nCells:", ncol(obj_small),
  "\n"
)


############################################################
# 5. Inspect available assays
############################################################

# Show which assays are stored in the Seurat object.
# In this object we expect to see the SCT assay.
Assays(obj_small)

# Show which assay is currently being used as the default assay.
DefaultAssay(obj_small)

# Inspect the structure of the default assay.
# This provides information such as the number of features,
# cells, variable features and SCT models.
obj_small[[DefaultAssay(obj_small)]]


############################################################
# 6. Inspect cell metadata
############################################################

# Metadata contains descriptive information about each cell.
# Each row represents one cell.
# Each column represents one piece of information about that cell.

# Show the first six cells and their metadata.
head(obj_small@meta.data)

# Show the names of all metadata variables.
colnames(obj_small@meta.data)

# Check the dimensions of the metadata table.
# We expect 7,788 rows because there are 7,788 cells.
dim(obj_small@meta.data)

# Examine the structure and data type of each metadata variable.
str(obj_small@meta.data)


############################################################
# 7. Examine available metadata variables
############################################################

# Print the metadata column names again so that we can identify
# variables related to sample, cohort, cluster and clinical group.
colnames(obj_small@meta.data)


############################################################
# 8. Inspect dimensionality reductions
############################################################

# Show dimensionality reductions that have already been
# calculated and stored in the Seurat object.
#
# In this dataset we expect:
# integrated.rpca = integrated reduced representation
# umap            = 2D representation used for visualization
Reductions(obj_small)


############################################################
# 9. Visualize cells by cohort using the existing UMAP
############################################################

# Create a UMAP plot where each dot represents one cell.
#
# The position of each cell comes from the UMAP coordinates
# already stored in obj_small.
#
# Cells are colored according to the "cohort" metadata column,
# allowing us to see the distribution of cRA and ICI cells.
#
# IMPORTANT:
# This does NOT calculate a new UMAP.
# It only visualizes the UMAP already stored in the object.

p_cohort <- DimPlot(
  obj_small,
  reduction = "umap",
  group.by = "cohort"
)


############################################################
# 10. Visualize cells by existing Seurat clusters
############################################################

# Create another UMAP using exactly the same cells and
# coordinates, but this time color cells according to their
# existing Seurat cluster assignment.
#
# "integrated.rpca_snn_res.0.8" contains the 29 existing
# clusters (0-28).
#
# label = TRUE adds cluster numbers to the plot.
# repel = TRUE helps prevent labels from overlapping.

p_cluster <- DimPlot(
  obj_small,
  reduction = "umap",
  group.by = "integrated.rpca_snn_res.0.8",
  label = TRUE,
  repel = TRUE
)


############################################################
# 11. Save the UMAP plots
############################################################

# Save the cohort UMAP as a PNG file in the figures folder.
ggplot2::ggsave(
  "../figures/UMAP_by_cohort.png",
  plot = p_cohort,
  width = 8,
  height = 6,
  dpi = 300
)

# Save the cluster UMAP as a PNG file in the figures folder.
ggplot2::ggsave(
  "../figures/UMAP_by_cluster.png",
  plot = p_cluster,
  width = 9,
  height = 7,
  dpi = 300
)


############################################################
# 12. Define a basic immune-cell marker panel
############################################################

# Create a list of well-known immune-cell marker genes.
#
# These markers will help us explore which broad immune
# populations may be represented in the existing clusters.
#
# IMPORTANT:
# We are NOT assigning final cell types yet.
# We are only examining gene-expression patterns.

markers <- c(

  # T-cell-associated markers
  "CD3D",
  "CD3E",

  # Markers useful when examining CD4/CD8 T-cell populations
  "CD4",
  "CD8A",

  # NK/cytotoxic-associated markers
  "NKG7",
  "GNLY",

  # B-cell-associated markers
  "MS4A1",
  "CD79A",

  # Plasma-cell-associated markers
  "MZB1",
  "JCHAIN",

  # Myeloid/monocyte-associated markers
  "LYZ",
  "S100A8",

  # Marker that can be informative for FCGR3A+ monocytes
  # and is also expressed in some NK populations
  "FCGR3A",

  # Broad myeloid/APC-associated marker
  "CST3",

  # Dendritic-cell-associated markers
  "FCER1A",
  "CD1C"
)


############################################################
# 13. Check whether the marker genes exist in the dataset
############################################################

# Check each marker against the genes present in obj_small.
# TRUE means the gene is present.
# FALSE means it is not available in this dataset.
markers %in% rownames(obj_small)

# Print any marker genes that are missing.
# If character(0) is returned, all markers are present.
setdiff(markers, rownames(obj_small))


############################################################
# 14. Examine marker-gene expression across clusters
############################################################

# DotPlot summarizes expression of the selected genes
# across the 29 existing clusters.
#
# Rows = clusters
# Columns = marker genes
#
# Dot size generally represents the percentage of cells in
# the cluster expressing the gene.
#
# Dot color represents average expression of the gene
# within that cluster.
#
# This helps us investigate the biological identity of
# clusters using multiple markers rather than relying on
# only one gene.

p_markers <- DotPlot(
  obj_small,
  features = markers,
  group.by = "integrated.rpca_snn_res.0.8"
) + RotatedAxis()


############################################################
# 15. Save the marker DotPlot
############################################################

# Save the marker-expression DotPlot in the figures folder.
ggplot2::ggsave(
  "../figures/DotPlot_basic_immune_markers.png",
  plot = p_markers,
  width = 12,
  height = 9,
  dpi = 300
)

############################################################
# 16. Define the Cytotoxic gene module
############################################################

# We will begin module scoring using ONE of the supplied
# 62 immune modules: the Cytotoxic module.
#
# The Cytotoxic module contains 8 genes in the supplied
# module library.

cytotoxic_genes <- c(
  "NKG7",
  "GNLY",
  "PRF1",
  "GZMB",
  "GZMH",
  "CTSW",
  "KLRD1",
  "FGFBP2"
)

# Print the module genes.
cytotoxic_genes


############################################################
# 17. Check whether all Cytotoxic module genes are present
############################################################

# Before calculating a module score, check whether each
# module gene exists in our Seurat object.
#
# TRUE  = gene is present
# FALSE = gene is missing

cytotoxic_genes %in% rownames(obj_small)


# Print only genes that are missing from the dataset.
#
# If this returns:
#
# character(0)
#
# then none of the genes are missing.

setdiff(
  cytotoxic_genes,
  rownames(obj_small)
)


############################################################
# 18. Calculate the Cytotoxic module score
############################################################

# AddModuleScore expects the gene modules as a LIST.
#
# Here we have only one module, so the list contains
# only our Cytotoxic gene vector.

cytotoxic_module <- list(
  cytotoxic_genes
)


# Calculate the Cytotoxic module activity for EVERY CELL.
#
# AddModuleScore compares the average expression of our
# Cytotoxic genes with matched control genes.
#
# IMPORTANT:
# This does NOT calculate one score for the whole dataset.
#
# Each of the 7,788 cells receives its OWN Cytotoxic score.

obj_small <- AddModuleScore(
  object = obj_small,
  features = cytotoxic_module,
  name = "Cytotoxic_Score"
)

############################################################
# 19. Inspect the new Cytotoxic module-score column
############################################################

# AddModuleScore stores the calculated score in the
# metadata of each cell.
#
# First, print all metadata column names again.
# We expect to see a new column called:
# Cytotoxic_Score1

colnames(obj_small@meta.data)


# Print only metadata columns containing the word "Cytotoxic".
# This makes it easier to confirm the exact name created
# by AddModuleScore.

grep(
  "Cytotoxic",
  colnames(obj_small@meta.data),
  value = TRUE
)

############################################################
# 20. Examine Cytotoxic scores for individual cells
############################################################

# Show the first 10 Cytotoxic module scores.
#
# IMPORTANT:
# Each number belongs to ONE CELL.
#
# Therefore, these are NOT cluster scores and NOT one score
# for the entire Cytotoxic module.

head(
  obj_small@meta.data$Cytotoxic_Score1,
  10
)

############################################################
# 21. Summarize Cytotoxic scores across all cells
############################################################

# We have 7,788 cells, so we also have 7,788 Cytotoxic scores.
#
# summary() gives us an overview of these scores:
#
# Min.    = lowest Cytotoxic score
# 1st Qu. = 25% of cells are below this value
# Median  = middle score
# Mean    = average score
# 3rd Qu. = 75% of cells are below this value
# Max.    = highest Cytotoxic score
#
# This helps us understand the score distribution before
# deciding which cells have relatively high or low activity.

summary(
  obj_small@meta.data$Cytotoxic_Score1
)
############################################################
# 22. Visualize Cytotoxic module activity on the UMAP
############################################################

# FeaturePlot places the Cytotoxic module score onto the
# existing UMAP.
#
# Each dot is still ONE CELL.
#
# The UMAP coordinates have NOT changed.
# Only the color of each cell now represents its
# Cytotoxic AddModuleScore.

p_cytotoxic_umap <- FeaturePlot(
  obj_small,
  features = "Cytotoxic_Score1",
  reduction = "umap"
)

# Display the plot.
print(p_cytotoxic_umap)


############################################################
# 23. Save the Cytotoxic-score UMAP
############################################################

ggplot2::ggsave(
  "../figures/UMAP_Cytotoxic_AddModuleScore.png",
  plot = p_cytotoxic_umap,
  width = 8,
  height = 6,
  dpi = 300
)