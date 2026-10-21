
# ============================================================
# Script 15: VAM scoring, diagnostics and benchmarking
# Dataset: WuImmune
# Modules: 62 original immune gene modules
# Method: VAM (Variance-adjusted Mahalanobis)
# Package: VAM 1.1.0
#
# Configuration:
#   center = FALSE
#   gamma = FALSE
#   sample.cov = FALSE
#
# This version investigates missing VAM scores and
# does not save invalid results as completed scores.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(VAM)
  library(Matrix)
})

dir.create("results", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 1. Load dataset and gene modules
# ------------------------------------------------------------

cat("Loading WuImmune dataset...\n")

wu <- readRDS("data/WuImmune_reclustered.rds")
modules <- readRDS("data/WuImmune_62_modules_available.rds")

DefaultAssay(wu) <- "RNA"

stopifnot(
  inherits(wu, "Seurat"),
  is.list(modules),
  length(modules) == 62,
  !is.null(names(modules)),
  !anyDuplicated(names(modules)),
  ncol(wu) == 50063
)

expr <- SeuratObject::GetAssayData(
  wu,
  assay = "RNA",
  layer = "data"
)

cell_ids <- colnames(wu)
module_names <- names(modules)

stopifnot(
  identical(colnames(expr), cell_ids),
  all(lengths(modules) > 0)
)

modules <- lapply(
  modules,
  function(genes) {
    unique(intersect(genes, rownames(expr)))
  }
)

if (any(lengths(modules) == 0)) {
  stop("One or more modules contain no available genes.")
}

cat("Cells:", length(cell_ids), "\n")
cat("Genes:", nrow(expr), "\n")
cat("Modules:", length(modules), "\n")
cat(
  "Module size range:",
  min(lengths(modules)), "to",
  max(lengths(modules)), "\n"
)

# ------------------------------------------------------------
# 2. Check technical variance
# ------------------------------------------------------------

cat("\nChecking technical variance...\n")

tech_var <- VAM:::getTechVarProp(wu)

stopifnot(
  length(tech_var) == nrow(expr),
  all(is.finite(tech_var)),
  all(tech_var > 0)
)

cat("Technical variance: valid\n")
cat("Minimum:", min(tech_var), "\n")
cat("Maximum:", max(tech_var), "\n")

# ------------------------------------------------------------
# 3. Run official VAM Seurat function
# ------------------------------------------------------------

cat("\nStarting VAM scoring...\n")

gc()

scoring_start <- proc.time()[["elapsed"]]

wu_scored <- VAM::vamForSeurat(
  seurat.data = wu,
  gene.set.collection = modules,
  center = FALSE,
  gamma = FALSE,
  sample.cov = FALSE,
  return.dist = FALSE
)

scoring_seconds <- (
  proc.time()[["elapsed"]] - scoring_start
)

cat("\nVAM scoring completed.\n")
cat(
  "Scoring time:",
  round(scoring_seconds, 2),
  "seconds\n"
)

# ------------------------------------------------------------
# 4. Inspect VAM assay
# ------------------------------------------------------------

cat("\n========== VAM ASSAY INSPECTION ==========\n")

if (!"VAMcdf" %in% Assays(wu_scored)) {
  stop("VAMcdf assay was not created.")
}

vam_assay <- wu_scored[["VAMcdf"]]

cat("Assay class:\n")
print(class(vam_assay))

cat("\nAvailable layers:\n")
print(SeuratObject::Layers(vam_assay))

cat("\nAssay dimensions:\n")
print(dim(vam_assay))

cat("\nFirst five feature names:\n")
print(head(rownames(vam_assay), 5))

# ------------------------------------------------------------
# 5. Inspect counts and data layers
# ------------------------------------------------------------

cat("\n========== VAM LAYER INSPECTION ==========\n")

vam_counts <- SeuratObject::GetAssayData(
  wu_scored,
  assay = "VAMcdf",
  layer = "counts"
)

cat("\nCounts layer dimensions:\n")
print(dim(vam_counts))

cat("\nCounts: first module, first 20 cells:\n")
print(as.numeric(vam_counts[1, 1:20]))

cat("\nCounts: missing values:\n")
print(sum(is.na(vam_counts)))

cat("\nCounts: non-finite values:\n")
print(sum(!is.finite(vam_counts)))

if ("data" %in% SeuratObject::Layers(vam_assay)) {

  vam_data <- SeuratObject::GetAssayData(
    wu_scored,
    assay = "VAMcdf",
    layer = "data"
  )

  cat("\nData layer dimensions:\n")
  print(dim(vam_data))

  cat("\nData: first module, first 20 cells:\n")
  print(as.numeric(vam_data[1, 1:20]))

  cat("\nData: missing values:\n")
  print(sum(is.na(vam_data)))

  cat("\nData: non-finite values:\n")
  print(sum(!is.finite(vam_data)))

} else {
  cat("\nNo data layer in VAMcdf assay.\n")
}

# ------------------------------------------------------------
# 6. Directly test first module with VAM::vam()
# ------------------------------------------------------------

cat("\n========== DIRECT VAM TEST ==========\n")

test_genes <- modules[[1]]

test_expr <- Matrix::t(
  expr[test_genes, , drop = FALSE]
)

test_tech_var <- tech_var[
  match(test_genes, rownames(expr))
]

direct_result <- VAM::vam(
  gene.expr = test_expr,
  tech.var.prop = test_tech_var,
  center = FALSE,
  gamma = FALSE
)

cat("Module:", module_names[1], "\n")
cat("Genes:", length(test_genes), "\n")

cat("\nDirect VAM distance summary:\n")
print(summary(direct_result$distance.sq))

cat("\nDirect VAM CDF summary:\n")
print(summary(direct_result$cdf.value))

cat(
  "Direct VAM missing scores:",
  sum(is.na(direct_result$cdf.value)),
  "\n"
)

cat(
  "Direct VAM non-finite scores:",
  sum(!is.finite(direct_result$cdf.value)),
  "\n"
)

# ------------------------------------------------------------
# 7. Align module and cell names
# ------------------------------------------------------------

cat("\n========== SCORE EXTRACTION ==========\n")

stopifnot(
  nrow(vam_counts) == 62,
  ncol(vam_counts) == 50063
)

# Seurat replaces underscores with dashes.

seurat_module_names <- gsub(
  "_", "-", module_names, fixed = TRUE
)

if (anyDuplicated(seurat_module_names)) {
  stop("Module names are not unique after conversion.")
}

if (!setequal(
  rownames(vam_counts),
  seurat_module_names
)) {
  stop("VAM output module names do not match.")
}

if (!setequal(
  colnames(vam_counts),
  cell_ids
)) {
  stop("VAM output cell IDs do not match.")
}

vam_counts <- vam_counts[
  seurat_module_names,
  cell_ids,
  drop = FALSE
]

rownames(vam_counts) <- module_names

score_matrix <- t(as.matrix(vam_counts))

stopifnot(
  identical(rownames(score_matrix), cell_ids),
  identical(colnames(score_matrix), module_names)
)

# ------------------------------------------------------------
# 8. Diagnose missing scores
# ------------------------------------------------------------

cat("\n========== VAM DIAGNOSTICS ==========\n")

na_counts <- colSums(is.na(score_matrix))
nonfinite_counts <- colSums(!is.finite(score_matrix))

diagnostics <- data.frame(
  module = module_names,
  n_genes = as.integer(lengths(modules)),
  missing_scores = as.integer(na_counts),
  nonfinite_scores = as.integer(nonfinite_counts),
  row.names = NULL
)

write.csv(
  diagnostics,
  "results/WuImmune_VAM_diagnostics.csv",
  row.names = FALSE
)

cat(
  "Total missing scores:",
  sum(na_counts), "\n"
)

cat(
  "Modules with missing scores:",
  sum(na_counts > 0), "\n"
)

cat(
  "Total non-finite scores:",
  sum(nonfinite_counts), "\n"
)

cat(
  "Modules with non-finite scores:",
  sum(nonfinite_counts > 0), "\n"
)

if (any(nonfinite_counts > 0)) {

  cat("\nAffected modules:\n")

  print(
    diagnostics[
      diagnostics$nonfinite_scores > 0,
      ,
      drop = FALSE
    ],
    row.names = FALSE
  )

  cat("\nDiagnostic report saved.\n")

  stop(
    "VAM produced non-finite scores. ",
    "Results were not saved as completed outputs."
  )
}

# ------------------------------------------------------------
# 9. Check score quality
# ------------------------------------------------------------

cat("\n========== VAM SCORE QUALITY ==========\n")

score_min <- min(score_matrix)
score_max <- max(score_matrix)

cat("Minimum score:", score_min, "\n")
cat("Maximum score:", score_max, "\n")

if (score_min < 0 || score_max > 1) {
  stop("VAM scores fall outside [0, 1].")
}

constant_modules <- colnames(score_matrix)[
  apply(
    score_matrix,
    2,
    function(x) all(x == x[1])
  )
]

all_zero_modules <- colnames(score_matrix)[
  colSums(score_matrix != 0) == 0
]

cat(
  "Constant-score modules:",
  length(constant_modules), "\n"
)

cat(
  "All-zero modules:",
  length(all_zero_modules), "\n"
)

if (length(constant_modules) > 0) {
  cat("\nConstant modules:\n")
  print(constant_modules)
}

if (length(all_zero_modules) > 0) {
  cat("\nAll-zero modules:\n")
  print(all_zero_modules)
}

# ------------------------------------------------------------
# 10. Save validated scores
# ------------------------------------------------------------

score_df <- data.frame(
  cell_id = cell_ids,
  score_matrix,
  check.names = FALSE,
  row.names = NULL
)

colnames(score_df)[-1] <- paste0(
  "VAM__",
  module_names
)

write.csv(
  score_df,
  "results/WuImmune_VAM_scores.csv",
  row.names = FALSE
)

saveRDS(
  score_df,
  "results/WuImmune_VAM_scores.rds"
)

saveRDS(
  wu_scored,
  "data/WuImmune_VAM_scored.rds"
)

# ------------------------------------------------------------
# 11. Save quality report
# ------------------------------------------------------------

quality_report <- data.frame(
  module = module_names,
  n_genes = as.integer(lengths(modules)),
  min_score = apply(score_matrix, 2, min),
  max_score = apply(score_matrix, 2, max),
  mean_score = colMeans(score_matrix),
  constant_score = module_names %in% constant_modules,
  all_zero = module_names %in% all_zero_modules,
  row.names = NULL
)

write.csv(
  quality_report,
  "results/WuImmune_VAM_quality_report.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 12. Save benchmark
# ------------------------------------------------------------

benchmark <- data.frame(
  method = "VAM",
  package_version = as.character(packageVersion("VAM")),
  cells = length(cell_ids),
  modules = length(modules),
  scoring_seconds = round(scoring_seconds, 2),
  scoring_minutes = round(scoring_seconds / 60, 2),
  installation_difficulty = "Easy",
  installation_seconds = NA_real_,
  peak_ram_gib = NA_real_,
  peak_memory_footprint_gib = NA_real_,
  total_execution_seconds = NA_real_,
  technical_variance = "Gene-specific VST estimates",
  center = FALSE,
  gamma = FALSE,
  sample_cov = FALSE,
  constant_score_modules = length(constant_modules),
  all_zero_modules = length(all_zero_modules),
  timestamp = as.character(Sys.time())
)

write.csv(
  benchmark,
  "results/benchmark_VAM.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 13. Final summary
# ------------------------------------------------------------

cat("\n========== VAM SUMMARY ==========\n")

cat("Cells scored:", nrow(score_matrix), "\n")
cat("Modules scored:", ncol(score_matrix), "\n")
cat("Missing scores:", sum(is.na(score_matrix)), "\n")

cat(
  "Scoring time (seconds):",
  round(scoring_seconds, 2), "\n"
)

cat(
  "Scoring time (minutes):",
  round(scoring_seconds / 60, 2), "\n"
)

cat(
  "Constant-score modules:",
  length(constant_modules), "\n"
)

cat(
  "All-zero modules:",
  length(all_zero_modules), "\n"
)

cat("\nSaved files:\n")
cat("results/WuImmune_VAM_scores.csv\n")
cat("results/WuImmune_VAM_scores.rds\n")
cat("data/WuImmune_VAM_scored.rds\n")
cat("results/benchmark_VAM.csv\n")
cat("results/WuImmune_VAM_diagnostics.csv\n")
cat("results/WuImmune_VAM_quality_report.csv\n")

cat("\nScript 15 completed successfully!\n")
