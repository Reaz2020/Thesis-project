
# ============================================================
# Script 11: Gene module scoring using Seurat AddModuleScore
# Dataset: WuImmune
# Modules: 62 WuImmune-compatible immune gene modules
# ============================================================

library(Seurat)

# ------------------------------------------------------------
# 1. Set up output folders
# ------------------------------------------------------------

dir.create("results", showWarnings = FALSE)
dir.create("data", showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Load WuImmune dataset and refined gene modules
# ------------------------------------------------------------

wu <- readRDS("data/WuImmune_reclustered.rds")

modules <- readRDS(
  "data/WuImmune_62_modules_available.rds"
)

cat("Number of cells:", ncol(wu), "\n")
cat("Number of modules:", length(modules), "\n")

# ------------------------------------------------------------
# 3. Check the inputs before scoring
# ------------------------------------------------------------

stopifnot(
  inherits(wu, "Seurat"),
  length(modules) == 62,
  !is.null(names(modules)),
  !anyDuplicated(names(modules)),
  all(lengths(modules) > 0)
)

DefaultAssay(wu) <- "RNA"

# Check that the normalized RNA data layer exists
rna_layers <- Layers(wu[["RNA"]])

if (!"data" %in% rna_layers) {
  stop(
    "Normalized RNA data layer is missing. ",
    "Check the preprocessing from Script 07."
  )
}

# Confirm that every module gene exists in the RNA assay
available_genes <- rownames(wu[["RNA"]])

missing_genes <- setdiff(
  unique(unlist(modules, use.names = FALSE)),
  available_genes
)

if (length(missing_genes) > 0) {
  stop(
    "Some module genes are missing from WuImmune: ",
    paste(missing_genes, collapse = ", ")
  )
}

cat("Input checks passed.\n")

# ------------------------------------------------------------
# 4. Calculate scores for all 62 modules
# ------------------------------------------------------------

set.seed(123)

wu <- AddModuleScore(
  object = wu,
  features = unname(modules),
  assay = "RNA",
  name = "AMS_",
  seed = 123
)

cat("AddModuleScore calculation completed.\n")

# ------------------------------------------------------------
# 5. Identify the newly generated score columns
# ------------------------------------------------------------

score_columns <- paste0("AMS_", seq_along(modules))

if (!all(score_columns %in% colnames(wu@meta.data))) {
  stop("Expected AddModuleScore columns were not found.")
}

# ------------------------------------------------------------
# 6. Create a separate score table with module names
# ------------------------------------------------------------

scores <- wu@meta.data[, score_columns, drop = FALSE]

colnames(scores) <- paste0(
  "AddModuleScore__",
  make.unique(names(modules))
)

scores <- data.frame(
  cell_id = rownames(scores),
  scores,
  check.names = FALSE,
  row.names = NULL
)

stopifnot(
  nrow(scores) == ncol(wu),
  ncol(scores) == 63,
  !anyNA(scores[, -1, drop = FALSE]),
  identical(scores$cell_id, colnames(wu))
)

# ------------------------------------------------------------
# 7. Save results
# ------------------------------------------------------------

write.csv(
  scores,
  "results/WuImmune_AddModuleScore_scores.csv",
  row.names = FALSE
)

saveRDS(
  scores,
  "results/WuImmune_AddModuleScore_scores.rds"
)

saveRDS(
  wu,
  "data/WuImmune_AddModuleScore_scored.rds"
)

# Save a mapping between generated score columns and modules
module_mapping <- data.frame(
  score_column = score_columns,
  module_name = names(modules)
)

write.csv(
  module_mapping,
  "results/WuImmune_AddModuleScore_module_mapping.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 8. Print final checks
# ------------------------------------------------------------

cat("\n========== SCORING SUMMARY ==========\n")
cat("Cells scored:", nrow(scores), "\n")
cat("Modules scored:", ncol(scores) - 1, "\n")
cat("Missing score values:", sum(is.na(scores[, -1])), "\n")

cat("\nSaved files:\n")
cat("results/WuImmune_AddModuleScore_scores.csv\n")
cat("results/WuImmune_AddModuleScore_scores.rds\n")
cat("data/WuImmune_AddModuleScore_scored.rds\n")
cat("results/WuImmune_AddModuleScore_module_mapping.csv\n")

cat("\nScript 11 completed successfully!\n")
