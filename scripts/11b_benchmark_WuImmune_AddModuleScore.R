# ============================================================
# Script 11b: Benchmark Seurat AddModuleScore
# Dataset: WuImmune
# Modules: 62 original WuImmune-compatible immune modules
#
# Uses exactly the same scoring settings as Script 11.
# Does NOT overwrite original scoring results.
# ============================================================

library(Seurat)

script_start <- proc.time()[["elapsed"]]

dir.create("results", showWarnings = FALSE)

cat("========================================\n")
cat("ADDMODULESCORE BENCHMARK\n")
cat("========================================\n")

# ------------------------------------------------------------
# 1. Load dataset and modules
# ------------------------------------------------------------

wu <- readRDS("data/WuImmune_reclustered.rds")

modules <- readRDS(
  "data/WuImmune_62_modules_available.rds"
)

cat("Seurat version:",
    as.character(packageVersion("Seurat")), "\n")

cat("Genes:", nrow(wu), "\n")
cat("Cells:", ncol(wu), "\n")
cat("Modules:", length(modules), "\n")

# ------------------------------------------------------------
# 2. Validate inputs
# ------------------------------------------------------------

stopifnot(
  inherits(wu, "Seurat"),
  ncol(wu) == 50063,
  length(modules) == 62,
  !is.null(names(modules)),
  !anyDuplicated(names(modules)),
  all(lengths(modules) > 0)
)

DefaultAssay(wu) <- "RNA"

if (!"data" %in% Layers(wu[["RNA"]])) {
  stop("Normalized RNA data layer is missing.")
}

available_genes <- rownames(wu[["RNA"]])

missing_genes <- setdiff(
  unique(unlist(modules, use.names = FALSE)),
  available_genes
)

if (length(missing_genes) > 0) {
  stop(
    "Missing module genes: ",
    paste(missing_genes, collapse = ", ")
  )
}

cat("Input validation passed.\n")

# ------------------------------------------------------------
# 3. Run AddModuleScore
# Exactly the same call as Script 11
# ------------------------------------------------------------

set.seed(123)

cat("\nStarting AddModuleScore...\n")

gc()

scoring_start <- proc.time()[["elapsed"]]

wu <- AddModuleScore(
  object = wu,
  features = unname(modules),
  assay = "RNA",
  name = "AMS_",
  seed = 123
)

scoring_seconds <- proc.time()[["elapsed"]] - scoring_start

cat("Scoring completed.\n")
cat("Scoring seconds:", round(scoring_seconds, 2), "\n")

# ------------------------------------------------------------
# 4. Extract module scores
# ------------------------------------------------------------

score_columns <- paste0("AMS_", seq_along(modules))

if (!all(score_columns %in% colnames(wu@meta.data))) {
  stop("Expected AddModuleScore columns are missing.")
}

scores <- wu@meta.data[
  ,
  score_columns,
  drop = FALSE
]

colnames(scores) <- paste0(
  "AddModuleScore__",
  make.unique(names(modules))
)

# ------------------------------------------------------------
# 5. Validate score quality
# ------------------------------------------------------------

stopifnot(
  nrow(scores) == 50063,
  ncol(scores) == 62,
  identical(rownames(scores), colnames(wu))
)

score_matrix <- as.matrix(scores)

missing_scores <- sum(is.na(score_matrix))

nonfinite_scores <- sum(!is.finite(score_matrix))

constant_score_modules <- sum(
  apply(
    score_matrix,
    2,
    function(x) length(unique(x)) <= 1
  )
)

cat("\nScore validation:\n")
cat("Cells scored:", nrow(score_matrix), "\n")
cat("Modules scored:", ncol(score_matrix), "\n")
cat("Missing scores:", missing_scores, "\n")
cat("Non-finite scores:", nonfinite_scores, "\n")
cat("Constant-score modules:", constant_score_modules, "\n")

if (
  missing_scores > 0 ||
  nonfinite_scores > 0 ||
  constant_score_modules > 0
) {
  stop("Score quality validation failed.")
}

# ------------------------------------------------------------
# 6. Compare against original Script 11 results
# ------------------------------------------------------------

original_path <- "results/WuImmune_AddModuleScore_scores.rds"

if (!file.exists(original_path)) {
  stop("Original AddModuleScore score file was not found.")
}

original <- readRDS(original_path)

stopifnot(
  identical(original$cell_id, colnames(wu)),
  identical(colnames(original)[-1], colnames(scores))
)

original_matrix <- as.matrix(
  original[, -1, drop = FALSE]
)

max_absolute_difference <- max(
  abs(score_matrix - original_matrix)
)

cat("\nComparison with original scores:\n")
cat(
  "Maximum absolute difference:",
  max_absolute_difference,
  "\n"
)

if (
  !is.finite(max_absolute_difference) ||
  max_absolute_difference > 1e-8
) {
  stop(
    "Benchmark scores do not match original Script 11 scores."
  )
}

# ------------------------------------------------------------
# 7. Calculate total execution time
# ------------------------------------------------------------

total_execution_seconds <-
  proc.time()[["elapsed"]] - script_start

# ------------------------------------------------------------
# 8. Save benchmark CSV
# ------------------------------------------------------------

benchmark <- data.frame(
  method = "AddModuleScore",
  package_version = as.character(packageVersion("Seurat")),
  cells = nrow(score_matrix),
  modules = ncol(score_matrix),
  genes = nrow(wu),

  scoring_seconds = round(scoring_seconds, 2),
  scoring_minutes = round(scoring_seconds / 60, 2),
  total_execution_seconds = round(
    total_execution_seconds, 2
  ),

  installation_difficulty = "Easy",
  installation_seconds = NA_real_,

  peak_ram_gib = NA_real_,
  peak_memory_footprint_gib = NA_real_,

  input_expression = "Seurat RNA data layer",
  seed = 123,

  missing_scores = missing_scores,
  nonfinite_scores = nonfinite_scores,
  constant_score_modules = constant_score_modules,

  max_absolute_difference = max_absolute_difference,

  status = "Completed",
  timestamp = as.character(Sys.time())
)

write.csv(
  benchmark,
  "results/benchmark_AddModuleScore.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 9. Print final summary
# ------------------------------------------------------------

cat("\n========================================\n")
cat("ADDMODULESCORE FINAL SUMMARY\n")
cat("========================================\n")

cat("Cells scored:", nrow(score_matrix), "\n")
cat("Modules scored:", ncol(score_matrix), "\n")
cat("Missing scores:", missing_scores, "\n")
cat("Non-finite scores:", nonfinite_scores, "\n")
cat("Constant-score modules:", constant_score_modules, "\n")

cat("Scoring seconds:", round(scoring_seconds, 2), "\n")
cat(
  "Total script seconds:",
  round(total_execution_seconds, 2),
  "\n"
)

cat(
  "Maximum difference from original:",
  max_absolute_difference,
  "\n"
)

cat("\nBenchmark saved:\n")
cat("results/benchmark_AddModuleScore.csv\n")

cat("\nScript 11b completed successfully!\n")