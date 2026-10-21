
# ============================================================
# Script 14: JASMINE scoring and benchmarking
# Dataset: WuImmune
# Modules: 62 immune gene modules
# Method: JASMINE oddsratio
# ============================================================

library(Seurat)

source(
  "external/JASMINE/JASMINE_V1_11October2021.r"
)

dir.create("results", showWarnings = FALSE)

# 1. Load data and modules
wu <- readRDS("data/WuImmune_reclustered.rds")

modules <- readRDS(
  "data/WuImmune_62_modules_available.rds"
)

DefaultAssay(wu) <- "RNA"

expr <- SeuratObject::GetAssayData(
  wu,
  assay = "RNA",
  layer = "data"
)

stopifnot(
  ncol(expr) == 50063,
  length(modules) == 62,
  all(unlist(modules, use.names = FALSE) %in% rownames(expr))
)

cell_ids <- colnames(expr)
module_names <- names(modules)

# 2. Prepare raw component matrices
raw_rank <- matrix(
  NA_real_,
  nrow = length(cell_ids),
  ncol = length(modules),
  dimnames = list(cell_ids, module_names)
)

raw_or <- raw_rank

batch_size <- 500
n_batches <- ceiling(length(cell_ids) / batch_size)

cat("Cells:", length(cell_ids), "\n")
cat("Modules:", length(modules), "\n")
cat("Batches:", n_batches, "\n")

gc()
start_time <- Sys.time()

# 3. Calculate raw components in batches
for (b in seq_len(n_batches)) {

  first <- (b - 1) * batch_size + 1
  last <- min(b * batch_size, length(cell_ids))
  idx <- first:last

  cat("Batch", b, "of", n_batches, "\n")

  batch_expr <- as.matrix(
    expr[, idx, drop = FALSE]
  )

  for (m in seq_along(modules)) {

    genes <- modules[[m]]

    if (length(genes) < 2) {
      stop(
        "JASMINE requires at least two genes: ",
        module_names[m]
      )
    }

    rank_values <- apply(
      batch_expr,
      2,
      function(x) {
        names(x) <- rownames(batch_expr)
        RankCalculation(x, genes)
      }
    )

    or_values <- ORCalculation(
      batch_expr,
      genes
    )

    raw_rank[idx, m] <- rank_values
    raw_or[idx, m] <- or_values
  }

  rm(batch_expr)
}

# 4. Normalize once across ALL cells
safe_minmax <- function(x) {

  if (any(!is.finite(x))) {
    stop("Non-finite raw JASMINE component detected.")
  }

  range_x <- range(x)

  if (range_x[1] == range_x[2]) {
    # Constant components have no between-cell variation.
    return(rep(0, length(x)))
  }

  (x - range_x[1]) /
    (range_x[2] - range_x[1])
}

score_matrix <- matrix(
  NA_real_,
  nrow = length(cell_ids),
  ncol = length(modules),
  dimnames = list(cell_ids, module_names)
)

for (m in seq_along(modules)) {

  normalized_rank <- safe_minmax(raw_rank[, m])
  normalized_or <- safe_minmax(raw_or[, m])

  score_matrix[, m] <- (
    normalized_rank + normalized_or
  ) / 2
}

elapsed_seconds <- as.numeric(
  difftime(Sys.time(), start_time, units = "secs")
)

# 5. Validate scores
stopifnot(
  nrow(score_matrix) == 50063,
  ncol(score_matrix) == 62,
  !anyNA(score_matrix),
  all(is.finite(score_matrix))
)

# 6. Save score table
scores <- data.frame(
  cell_id = cell_ids,
  score_matrix,
  check.names = FALSE,
  row.names = NULL
)

colnames(scores)[-1] <- paste0(
  "JASMINE__",
  module_names
)

write.csv(
  scores,
  "results/WuImmune_JASMINE_scores.csv",
  row.names = FALSE
)

saveRDS(
  scores,
  "results/WuImmune_JASMINE_scores.rds"
)

# 7. Save scored Seurat object
metadata_scores <- as.data.frame(
  score_matrix,
  check.names = FALSE
)

colnames(metadata_scores) <- paste0(
  module_names,
  "_JASMINE"
)

rownames(metadata_scores) <- cell_ids

wu <- SeuratObject::AddMetaData(
  wu,
  metadata = metadata_scores
)

saveRDS(
  wu,
  "data/WuImmune_JASMINE_scored.rds"
)

# 8. Save separate benchmark
benchmark <- data.frame(
  method = "JASMINE",
  variant = "oddsratio",
  cells = length(cell_ids),
  modules = length(modules),
  batch_size = batch_size,
  scoring_seconds = round(elapsed_seconds, 2),
  scoring_minutes = round(elapsed_seconds / 60, 2),
  installation_difficulty = "Moderate",
  installation_seconds = NA_real_,
  peak_ram_gib = NA_real_,
  peak_memory_footprint_gib = NA_real_,
  total_execution_seconds = NA_real_,
  timestamp = as.character(Sys.time())
)

write.csv(
  benchmark,
  "results/benchmark_JASMINE.csv",
  row.names = FALSE
)

cat("\n========== JASMINE SUMMARY ==========\n")
cat("Cells scored:", nrow(scores), "\n")
cat("Modules scored:", ncol(scores) - 1, "\n")
cat("Missing scores:", sum(is.na(score_matrix)), "\n")
cat("Scoring time (seconds):", round(elapsed_seconds, 2), "\n")
cat("Scoring time (minutes):", round(elapsed_seconds / 60, 2), "\n")
cat("Method: oddsratio\n")

cat("\nSaved files:\n")
cat("results/WuImmune_JASMINE_scores.csv\n")
cat("results/WuImmune_JASMINE_scores.rds\n")
cat("data/WuImmune_JASMINE_scored.rds\n")
cat("results/benchmark_JASMINE.csv\n")

cat("\nScript 14 completed successfully!\n")
