
# ============================================================
# Script 16: VISION gene-module scoring
#
# Thesis: Explainable Gene Module-Based Framework
# Dataset: WuImmune
# Method: VISION
#
# Inputs:
#   data/WuImmune_reclustered.rds
#   data/WuImmune_62_modules_available.rds
#
# Outputs:
#   results/WuImmune_VISION_scores.csv
#   results/WuImmune_VISION_scores.rds
#   results/WuImmune_VISION_quality_report.csv
#   results/benchmark_VISION.csv
#
# Memory usage is measured externally using /usr/bin/time -l.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(VISION)
  library(Matrix)
})

# ------------------------------------------------------------
# 1. Setup
# ------------------------------------------------------------

dir.create(
  "results",
  recursive = TRUE,
  showWarnings = FALSE
)

set.seed(123)

total_start <- proc.time()[["elapsed"]]

cat("========================================\n")
cat("VISION gene-module scoring\n")
cat("========================================\n")

cat(
  "VISION version:",
  as.character(packageVersion("VISION")),
  "\n"
)

# ------------------------------------------------------------
# 2. Load dataset and modules
# ------------------------------------------------------------

cat("\nLoading WuImmune dataset...\n")

wu <- readRDS(
  "data/WuImmune_reclustered.rds"
)

modules <- readRDS(
  "data/WuImmune_62_modules_available.rds"
)

stopifnot(
  inherits(wu, "Seurat"),
  is.list(modules),
  length(modules) == 62,
  !is.null(names(modules)),
  !anyDuplicated(names(modules)),
  ncol(wu) == 50063
)

module_names <- names(modules)
cell_ids <- colnames(wu)

cat("Cells:", length(cell_ids), "\n")
cat("Modules:", length(modules), "\n")

# ------------------------------------------------------------
# 3. Extract normalized RNA expression
# ------------------------------------------------------------

cat("\nExtracting normalized RNA expression...\n")

expr <- SeuratObject::GetAssayData(
  wu,
  assay = "RNA",
  layer = "data"
)

stopifnot(
  ncol(expr) == length(cell_ids),
  identical(colnames(expr), cell_ids)
)

cat("Expression genes:", nrow(expr), "\n")
cat("Expression cells:", ncol(expr), "\n")

# Keep the same available genes used by other methods.

modules <- lapply(
  modules,
  function(genes) {
    unique(intersect(genes, rownames(expr)))
  }
)

names(modules) <- module_names

module_sizes <- lengths(modules)

if (any(module_sizes < 3)) {
  stop(
    "One or more modules have fewer than 3 genes."
  )
}

cat(
  "Module size range:",
  min(module_sizes),
  "to",
  max(module_sizes),
  "\n"
)

# ------------------------------------------------------------
# 4. Create VISION signatures
# ------------------------------------------------------------

cat("\nCreating 62 VISION signatures...\n")

signatures <- lapply(
  seq_along(modules),
  function(i) {

    genes <- modules[[i]]

    weights <- setNames(
      rep(1, length(genes)),
      genes
    )

    VISION::createGeneSignature(
      name = module_names[i],
      sigData = weights
    )
  }
)

names(signatures) <- module_names

cat(
  "Signatures created:",
  length(signatures),
  "\n"
)

# ------------------------------------------------------------
# 5. Create VISION object
# ------------------------------------------------------------

cat("\nCreating VISION object...\n")

vision <- VISION::Vision(expr)

cat("VISION object created.\n")

# ------------------------------------------------------------
# 6. Add signatures
# ------------------------------------------------------------

cat("\nAdding gene signatures...\n")

vision <- VISION::addSignatures(
  object = vision,
  signatures = signatures,
  min_signature_genes = 3,
  sig_gene_threshold = 0
)

# Check retained signatures before scoring.

retained_names <- names(vision@sigData)

cat(
  "Signatures retained:",
  length(retained_names),
  "\n"
)

if (
  length(retained_names) != 62 ||
  !setequal(retained_names, module_names)
) {
  cat("\nRetained signatures:\n")
  print(retained_names)

  stop(
    "VISION did not retain all 62 original modules. ",
    "Scoring stopped to prevent an incomplete comparison."
  )
}

# ------------------------------------------------------------
# 7. Calculate VISION signature scores
# ------------------------------------------------------------

cat("\nStarting VISION signature scoring...\n")

gc()

scoring_start <- proc.time()[["elapsed"]]

vision <- VISION::calcSignatureScores(
  object = vision,
  sig_norm_method = NULL,
  sig_gene_importance = TRUE
)

scoring_seconds <- (
  proc.time()[["elapsed"]] - scoring_start
)

cat("\nVISION scoring completed.\n")

cat(
  "Scoring time:",
  round(scoring_seconds, 2),
  "seconds\n"
)

# ------------------------------------------------------------
# 8. Extract scores
# ------------------------------------------------------------

cat("\nExtracting VISION scores...\n")

scores <- VISION::getSignatureScores(vision)

scores <- as.matrix(scores)

cat("\nRaw score dimensions:\n")
print(dim(scores))

cat("\nFirst score column names:\n")
print(head(colnames(scores)))

# The successful small test returned:
# rows = cells, columns = signatures.

if (
  nrow(scores) != length(cell_ids) ||
  ncol(scores) != length(module_names)
) {
  stop(
    "Unexpected VISION score dimensions."
  )
}

# ------------------------------------------------------------
# 9. Validate cell and module identities
# ------------------------------------------------------------

cat("\nValidating cell and module identities...\n")

if (
  is.null(rownames(scores)) ||
  !setequal(rownames(scores), cell_ids)
) {
  stop(
    "VISION score cell IDs do not match WuImmune."
  )
}

if (
  is.null(colnames(scores)) ||
  !setequal(colnames(scores), module_names)
) {
  stop(
    "VISION score module names do not match."
  )
}

# Restore exact original ordering.

scores <- scores[
  cell_ids,
  module_names,
  drop = FALSE
]

stopifnot(
  identical(rownames(scores), cell_ids),
  identical(colnames(scores), module_names)
)

cat("Cell IDs: matched\n")
cat("Module names: matched\n")
cat("Original ordering: restored\n")

# ------------------------------------------------------------
# 10. Validate score quality
# ------------------------------------------------------------

cat("\nChecking score quality...\n")

missing_scores <- sum(is.na(scores))

nonfinite_scores <- sum(
  !is.finite(scores)
)

cat(
  "Missing scores:",
  missing_scores,
  "\n"
)

cat(
  "Non-finite scores:",
  nonfinite_scores,
  "\n"
)

if (
  missing_scores > 0 ||
  nonfinite_scores > 0
) {
  stop(
    "VISION produced missing or non-finite scores. ",
    "Invalid outputs will not be saved."
  )
}

cat(
  "Minimum score:",
  min(scores),
  "\n"
)

cat(
  "Maximum score:",
  max(scores),
  "\n"
)

# Negative values are allowed for VISION.

constant_modules <- vapply(
  seq_len(ncol(scores)),
  function(i) {
    x <- scores[, i]
    all(x == x[1])
  },
  logical(1)
)

names(constant_modules) <- module_names

cat(
  "Constant-score modules:",
  sum(constant_modules),
  "\n"
)

if (any(constant_modules)) {
  cat("\nConstant-score module names:\n")
  print(names(constant_modules)[constant_modules])
}

# ------------------------------------------------------------
# 11. Create quality report
# ------------------------------------------------------------

quality_report <- data.frame(
  module = module_names,
  n_genes = as.integer(module_sizes),
  min_score = apply(scores, 2, min),
  max_score = apply(scores, 2, max),
  mean_score = colMeans(scores),
  sd_score = apply(scores, 2, sd),
  missing_scores = as.integer(
    colSums(is.na(scores))
  ),
  nonfinite_scores = as.integer(
    colSums(!is.finite(scores))
  ),
  constant_score = unname(constant_modules),
  row.names = NULL
)

# ------------------------------------------------------------
# 12. Save scores
# ------------------------------------------------------------

cat("\nSaving validated VISION scores...\n")

score_df <- data.frame(
  cell_id = cell_ids,
  scores,
  check.names = FALSE,
  row.names = NULL
)

colnames(score_df)[-1] <- paste0(
  "VISION__",
  module_names
)

write.csv(
  score_df,
  "results/WuImmune_VISION_scores.csv",
  row.names = FALSE
)

saveRDS(
  score_df,
  "results/WuImmune_VISION_scores.rds"
)

write.csv(
  quality_report,
  "results/WuImmune_VISION_quality_report.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 13. Benchmark
# ------------------------------------------------------------

total_seconds <- (
  proc.time()[["elapsed"]] - total_start
)

benchmark <- data.frame(
  method = "VISION",
  package_version = as.character(
    packageVersion("VISION")
  ),
  cells = length(cell_ids),
  modules = length(module_names),
  genes = nrow(expr),
  scoring_seconds = round(
    scoring_seconds, 2
  ),
  scoring_minutes = round(
    scoring_seconds / 60, 2
  ),
  total_execution_seconds = round(
    total_seconds, 2
  ),
  installation_difficulty = "Moderate",
  installation_seconds = NA_real_,
  peak_ram_gib = NA_real_,
  peak_memory_footprint_gib = NA_real_,
  input_expression = "Seurat RNA data layer",
  signature_min_genes = 3,
  sig_gene_threshold = 0.01,
  sig_norm_method = "NULL",
  sig_gene_importance = TRUE,
  constant_score_modules = sum(
    constant_modules
  ),
  timestamp = as.character(Sys.time()),
  row.names = NULL
)

write.csv(
  benchmark,
  "results/benchmark_VISION.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 14. Final summary
# ------------------------------------------------------------

cat("\n========================================\n")
cat("VISION FINAL SUMMARY\n")
cat("========================================\n")

cat(
  "Cells scored:",
  nrow(scores),
  "\n"
)

cat(
  "Modules scored:",
  ncol(scores),
  "\n"
)

cat(
  "Missing scores:",
  missing_scores,
  "\n"
)

cat(
  "Non-finite scores:",
  nonfinite_scores,
  "\n"
)

cat(
  "Constant-score modules:",
  sum(constant_modules),
  "\n"
)

cat(
  "Scoring seconds:",
  round(scoring_seconds, 2),
  "\n"
)

cat(
  "Total script seconds:",
  round(total_seconds, 2),
  "\n"
)

cat("\nSaved files:\n")

cat(
  "results/WuImmune_VISION_scores.csv\n"
)

cat(
  "results/WuImmune_VISION_scores.rds\n"
)

cat(
  "results/WuImmune_VISION_quality_report.csv\n"
)

cat(
  "results/benchmark_VISION.csv\n"
)

cat("\nVISION Script 16 completed successfully!\n")
