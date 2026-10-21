
# ============================================================
# Script 17: Full GSVA scoring - WuImmune
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(GSVA)
  library(Matrix)
})

dir.create("results", showWarnings = FALSE, recursive = TRUE)

set.seed(123)
total_start <- proc.time()[["elapsed"]]

cat("========================================\n")
cat("GSVA FULL SCORING\n")
cat("========================================\n")
cat("GSVA version:", as.character(packageVersion("GSVA")), "\n")

# ------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------

wu <- readRDS("data/WuImmune_reclustered.rds")
modules <- readRDS("data/WuImmune_62_modules_available.rds")

expr <- SeuratObject::GetAssayData(
  wu,
  assay = "RNA",
  layer = "data"
)

stopifnot(
  length(modules) == 62,
  ncol(expr) == 50063,
  !anyDuplicated(names(modules)),
  !anyDuplicated(colnames(expr))
)

module_names <- names(modules)
cell_ids <- colnames(expr)

cat("Genes:", nrow(expr), "\n")
cat("Cells:", ncol(expr), "\n")
cat("Modules:", length(modules), "\n")

# Confirm that every original module has at least one gene
# represented in the expression matrix.
available_counts <- vapply(
  modules,
  function(g) sum(g %in% rownames(expr)),
  integer(1)
)

if (any(available_counts == 0)) {
  stop("At least one module has no genes in the expression matrix.")
}

# ------------------------------------------------------------
# 2. Create GSVA parameter object
# ------------------------------------------------------------

cat("\nCreating GSVA parameters...\n")

param <- GSVA::gsvaParam(
  exprData = expr,
  geneSets = modules,
  kcdf = "Gaussian",
  minSize = 1,
  maxSize = Inf,
  sparse = TRUE
)

# ------------------------------------------------------------
# 3. Run GSVA
# ------------------------------------------------------------

cat("\nStarting GSVA scoring...\n")

gc()
scoring_start <- proc.time()[["elapsed"]]

scores <- GSVA::gsva(
  param,
  verbose = TRUE
)

scoring_seconds <- proc.time()[["elapsed"]] - scoring_start

cat("\nGSVA scoring completed.\n")
cat("Scoring seconds:", round(scoring_seconds, 2), "\n")

# ------------------------------------------------------------
# 4. Validate dimensions and identities
# ------------------------------------------------------------

scores <- as.matrix(scores)

cat("\nScore dimensions:\n")
print(dim(scores))

# GSVA returns modules x cells.
if (
  nrow(scores) != 62 ||
  ncol(scores) != length(cell_ids)
) {
  cat("\nMissing modules:\n")
  print(setdiff(module_names, rownames(scores)))
  stop("GSVA did not score all 62 modules.")
}

if (
  !setequal(rownames(scores), module_names) ||
  !setequal(colnames(scores), cell_ids)
) {
  stop("GSVA module names or cell IDs do not match.")
}

scores <- scores[
  module_names,
  cell_ids,
  drop = FALSE
]

missing_scores <- sum(is.na(scores))
nonfinite_scores <- sum(!is.finite(scores))

cat("Missing scores:", missing_scores, "\n")
cat("Non-finite scores:", nonfinite_scores, "\n")

if (missing_scores > 0 || nonfinite_scores > 0) {
  stop("Invalid GSVA scores detected. Results will not be saved.")
}

constant_modules <- apply(
  scores,
  1,
  function(x) all(x == x[1])
)

cat("Constant-score modules:", sum(constant_modules), "\n")

# ------------------------------------------------------------
# 5. Create quality report
# ------------------------------------------------------------

quality_report <- data.frame(
  module = module_names,
  original_genes = lengths(modules),
  available_genes = as.integer(available_counts),
  min_score = apply(scores, 1, min),
  max_score = apply(scores, 1, max),
  mean_score = rowMeans(scores),
  sd_score = apply(scores, 1, sd),
  missing_scores = rowSums(is.na(scores)),
  nonfinite_scores = rowSums(!is.finite(scores)),
  constant_score = as.logical(constant_modules),
  row.names = NULL
)

# ------------------------------------------------------------
# 6. Save validated scores
# ------------------------------------------------------------

cat("\nSaving GSVA scores...\n")

# Convert from modules x cells to cells x modules.
score_df <- data.frame(
  cell_id = cell_ids,
  t(scores),
  check.names = FALSE,
  row.names = NULL
)

colnames(score_df)[-1] <- paste0(
  "GSVA__",
  module_names
)

write.csv(
  score_df,
  "results/WuImmune_GSVA_scores.csv",
  row.names = FALSE
)

saveRDS(
  score_df,
  "results/WuImmune_GSVA_scores.rds"
)

write.csv(
  quality_report,
  "results/WuImmune_GSVA_quality_report.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 7. Benchmark
# ------------------------------------------------------------

total_seconds <- proc.time()[["elapsed"]] - total_start

benchmark <- data.frame(
  method = "GSVA",
  package_version = as.character(packageVersion("GSVA")),
  cells = ncol(expr),
  modules = length(modules),
  genes = nrow(expr),
  scoring_seconds = round(scoring_seconds, 2),
  scoring_minutes = round(scoring_seconds / 60, 2),
  total_execution_seconds = round(total_seconds, 2),
  installation_difficulty = "Easy",
  installation_seconds = NA_real_,
  peak_ram_gib = NA_real_,
  peak_memory_footprint_gib = NA_real_,
  input_expression = "Seurat RNA data layer",
  kcdf = "Gaussian",
  sparse = TRUE,
  min_size = 1,
  max_size = "Inf",
  missing_scores = missing_scores,
  nonfinite_scores = nonfinite_scores,
  constant_score_modules = sum(constant_modules),
  status = "Completed",
  timestamp = as.character(Sys.time()),
  row.names = NULL
)

write.csv(
  benchmark,
  "results/benchmark_GSVA.csv",
  row.names = FALSE
)

cat("\n========================================\n")
cat("GSVA FINAL SUMMARY\n")
cat("========================================\n")
cat("Cells scored:", ncol(scores), "\n")
cat("Modules scored:", nrow(scores), "\n")
cat("Missing scores:", missing_scores, "\n")
cat("Non-finite scores:", nonfinite_scores, "\n")
cat("Constant-score modules:", sum(constant_modules), "\n")
cat("Scoring seconds:", round(scoring_seconds, 2), "\n")
cat("Total script seconds:", round(total_seconds, 2), "\n")
cat("\nGSVA Script 17 completed successfully!\n")
