# Task 1 — exploring  WuImmune data cells, for meeting 3
# Load Seurat
library(Seurat)

# Load the WuImmune dataset
load("data/WuImmune.RData")

# Inspect the Seurat object
WuImmune

# Check number of genes and cells
dim(WuImmune)

# Inspect cell metadata
head(WuImmune@meta.data)

# Check available assays
Assays(WuImmune)
DefaultAssay(WuImmune)

# Check the major immune-cell annotations
table(WuImmune$celltype_major)

# Check the detailed immune-cell annotations
table(WuImmune$celltype_minor)