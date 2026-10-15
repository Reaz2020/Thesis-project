# ================================================================
# Task 2 — Recluster WuImmune immune cells for meeting 3
# Run from the project root:
#   Rscript scripts/07_recluster_WuImmune.R
#
# Input:  data/WuImmune.RData
# Output: data/WuImmune_reclustered.rds
# Plots:  figures/WuImmune_*.png
#
# IMPORTANT: Original cell-type annotations and original RData
# are preserved. This script creates a NEW processed object.
# ================================================================

# 1. Load Seurat: the R package used to analyze single-cell data.
suppressPackageStartupMessages(library(Seurat))

# 2. Create the figures folder if it does not exist.
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

# 3. Check that the input file exists before trying to load it.
input_file <- "data/WuImmune.RData"
if (!file.exists(input_file)) {
  stop("Cannot find ", input_file, ". Run this script from the thesis project root.")
}

# 4. Load the RData file; it should create an object named WuImmune.
load(input_file)

# 5. Confirm that the expected Seurat object was loaded.
if (!exists("WuImmune") || !inherits(WuImmune, "Seurat")) {
  stop("The RData file did not contain a Seurat object named 'WuImmune'.")
}

# 6. Confirm that the RNA assay exists, because Salim specified it.
if (!"RNA" %in% Assays(WuImmune)) {
  stop("RNA assay not found in WuImmune.")
}

# 7. Confirm that both original reference annotation columns exist.
reference_columns <- c("celltype_major", "celltype_minor")
if (!all(reference_columns %in% colnames(WuImmune@meta.data))) {
  stop("Missing original cell-type annotation column(s).")
}

# 8. Store original annotations so we can verify they stay unchanged.
original_major <- WuImmune$celltype_major
original_minor <- WuImmune$celltype_minor
original_cell_ids <- colnames(WuImmune)

# 9. Use the RNA assay for this preprocessing workflow.
DefaultAssay(WuImmune) <- "RNA"

# 10. Normalize gene-expression measurements across cells.
#     This creates/updates normalized expression in the RNA data layer.
message("Step 1/8: Normalizing RNA data...")
WuImmune <- NormalizeData(WuImmune, verbose = FALSE)

# 11. Select the 3,000 genes showing high variability across cells.
#     These are used to focus the following PCA on informative variation.
message("Step 2/8: Finding 3,000 variable genes...")
WuImmune <- FindVariableFeatures(
  WuImmune,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = FALSE
)

# 12. Scale the selected variable genes (center and standardize).
#     This is a common preparation step for PCA.
message("Step 3/8: Scaling variable genes...")
WuImmune <- ScaleData(
  WuImmune,
  features = VariableFeatures(WuImmune),
  verbose = FALSE
)

# 13. Recalculate PCA using the selected variable genes.
message("Step 4/8: Running PCA...")
WuImmune <- RunPCA(
  WuImmune,
  features = VariableFeatures(WuImmune),
  npcs = 50,
  verbose = FALSE
)

# 14. Save the elbow plot to inspect how much variation PCs capture.
#     Salim suggested starting with PCs 1–30 for later steps.
png("figures/WuImmune_elbow_plot.png", width = 1200, height = 800, res = 140)
print(ElbowPlot(WuImmune, ndims = 50))
dev.off()

# 15. Find cells with similar expression profiles in PCA space.
message("Step 5/8: Finding neighbors...")
WuImmune <- FindNeighbors(WuImmune, dims = 1:30, verbose = FALSE)

# 16. Group similar cells into new computational clusters.
#     Cluster numbers are NOT automatically biological cell-type names.
message("Step 6/8: Finding clusters...")
WuImmune <- FindClusters(WuImmune, resolution = 0.5, verbose = FALSE)

# 17. Calculate a NEW UMAP using the immune-cell-only PCA.
message("Step 7/8: Running UMAP...")
set.seed(42)  # Helps make the result reproducible.
WuImmune <- RunUMAP(WuImmune, dims = 1:30, seed.use = 42, verbose = FALSE)

# 18. Check that the original cell IDs and reference labels are unchanged.
stopifnot(
  identical(colnames(WuImmune), original_cell_ids),
  identical(WuImmune$celltype_major, original_major),
  identical(WuImmune$celltype_minor, original_minor)
)

# 19. Plot UMAP colored by the ORIGINAL major cell-type labels.
message("Step 8/8: Saving UMAP plots...")
p_major <- DimPlot(
  WuImmune, reduction = "umap", group.by = "celltype_major",
  label = TRUE, repel = TRUE
) + ggplot2::ggtitle("WuImmune — major reference cell types")
ggplot2::ggsave(
  "figures/WuImmune_umap_major.png",
  plot = p_major, width = 11, height = 8, dpi = 150
)

# 20. Plot UMAP colored by the ORIGINAL detailed cell-type labels.
p_minor <- DimPlot(
  WuImmune, reduction = "umap", group.by = "celltype_minor",
  label = TRUE, repel = TRUE
) + ggplot2::ggtitle("WuImmune — minor reference cell types")
ggplot2::ggsave(
  "figures/WuImmune_umap_minor.png",
  plot = p_minor, width = 13, height = 9, dpi = 150
)

# 21. Plot UMAP colored by the NEW computational clusters.
p_clusters <- DimPlot(
  WuImmune, reduction = "umap", group.by = "seurat_clusters",
  label = TRUE, repel = TRUE
) + ggplot2::ggtitle("WuImmune — new immune-cell clusters")
ggplot2::ggsave(
  "figures/WuImmune_umap_clusters.png",
  plot = p_clusters, width = 11, height = 8, dpi = 150
)

# 22. Save the processed object separately; do NOT overwrite the RData.
message("Saving reclustered object (this may take a while)...")
saveRDS(WuImmune, file = "data/WuImmune_reclustered.rds")

# 23. Print a short report so you can check the result in the terminal.
message("=== Task 2 completed ===")
print(dim(WuImmune))
cat("Number of new clusters:", length(unique(WuImmune$seurat_clusters)), "\n")
cat("Major reference annotations:\n")
print(table(WuImmune$celltype_major))
cat("Minor reference annotations:\n")
print(table(WuImmune$celltype_minor))
cat("Saved: data/WuImmune_reclustered.rds\n")
cat("Plots saved in: figures/\n")
