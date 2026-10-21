
# ============================================================
# Script 13: AUCell scoring and computational benchmarking
# Dataset: WuImmune
# Modules: 62 WuImmune-compatible modules
# ============================================================

library(Seurat)
library(AUCell)
library(BiocParallel)

dir.create("results", showWarnings = FALSE)

# 1. Load inputs
wu <- readRDS("data/WuImmune_reclustered.rds")
modules <- readRDS("data/WuImmune_62_modules_available.rds")

DefaultAssay(wu) <- "RNA"

stopifnot(
  ncol(wu) == 50063,
  length(modules) == 62,
  !anyDuplicated(names(modules)),
  all(lengths(modules) > 0)
)

# 2. Obtain sparse expression matrix
expr <- SeuratObject::GetAssayData(
  wu,
  assay = "RNA",
  layer = "data"
)

stopifnot(
  all(unlist(modules, use.names = FALSE) %in% rownames(expr))
)

# Keep the ranking threshold constant across all batches
auc_max_rank <- ceiling(0.05 * nrow(expr))
batch_size <- 500

cell_ids <- colnames(expr)
n_cells <- length(cell_ids)
n_batches <- ceiling(n_cells / batch_size)

cat("Cells:", n_cells, "\n")
cat("Genes:", nrow(expr), "\n")
cat("Modules:", length(modules), "\n")
cat("AUC maximum rank:", auc_max_rank, "\n")
cat("Batches:", n_batches, "\n")

# 3. Prepare output matrix
score_matrix <- matrix(
  NA_real_,
  nrow = n_cells,
  ncol = length(modules),
  dimnames = list(cell_ids, names(modules))
)

ranking_seconds <- 0
auc_seconds <- 0

gc()
total_start <- Sys.time()

# 4. Process cells in memory-conscious batches
set.seed(123)

for (batch in seq_len(n_batches)) {

  first <- (batch - 1) * batch_size + 1
  last <- min(batch * batch_size, n_cells)
  idx <- first:last

  cat(
    "\nBatch", batch, "of", n_batches,
    "- cells", first, "to", last, "\n"
  )

  batch_expr <- expr[, idx, drop = FALSE]

  # Step A: Build gene rankings
  rank_start <- Sys.time()

  rankings <- AUCell::AUCell_buildRankings(
    batch_expr,
    plotStats = FALSE,
    BPPARAM = BiocParallel::SerialParam(),
    verbose = FALSE
  )

  ranking_seconds <- ranking_seconds +
    as.numeric(difftime(
      Sys.time(), rank_start, units = "secs"
    ))

  # Step B: Calculate module AUC scores
  auc_start <- Sys.time()

  auc_result <- AUCell::AUCell_calcAUC(
    geneSets = modules,
    rankings = rankings,
    aucMaxRank = auc_max_rank,
    normAUC = TRUE,
    nCores = 1,
    verbose = FALSE
  )

  auc_seconds <- auc_seconds +
    as.numeric(difftime(
      Sys.time(), auc_start, units = "secs"
    ))

  auc_values <- AUCell::getAUC(auc_result)

  stopifnot(
    all(names(modules) %in% rownames(auc_values)),
    all(cell_ids[idx] %in% colnames(auc_values))
  )

  score_matrix[idx, ] <- t(
    auc_values[
      names(modules),
      cell_ids[idx],
      drop = FALSE
    ]
  )

  rm(batch_expr, rankings, auc_result, auc_values)
}

total_seconds <- as.numeric(
  difftime(Sys.time(), total_start, units = "secs")
)

# 5. Validate scores
stopifnot(
  nrow(score_matrix) == 50063,
  ncol(score_matrix) == 62,
  !anyNA(score_matrix),
  all(is.finite(score_matrix))
)

cat("\nAUCell scoring validated.\n")

# 6. Save score tables
scores <- data.frame(
  cell_id = cell_ids,
  score_matrix,
  check.names = FALSE,
  row.names = NULL
)

colnames(scores)[-1] <- paste0(
  "AUCell__",
  names(modules)
)

write.csv(
  scores,
  "results/WuImmune_AUCell_scores.csv",
  row.names = FALSE
)

saveRDS(
  scores,
  "results/WuImmune_AUCell_scores.rds"
)

# 7. Add scores to Seurat metadata
metadata_scores <- as.data.frame(
  score_matrix,
  check.names = FALSE
)

colnames(metadata_scores) <- paste0(
  names(modules),
  "_AUCell"
)

rownames(metadata_scores) <- cell_ids

wu <- SeuratObject::AddMetaData(
  object = wu,
  metadata = metadata_scores
)

saveRDS(
  wu,
  "data/WuImmune_AUCell_scored.rds"
)

# 8. Save separate computational benchmark
benchmark <- data.frame(
  method = "AUCell",
  cells = n_cells,
  modules = length(modules),
  genes = nrow(expr),
  batch_size = batch_size,
  batches = n_batches,
  auc_max_rank = auc_max_rank,
  ranking_seconds = round(ranking_seconds, 2),
  auc_seconds = round(auc_seconds, 2),
  scoring_seconds = round(total_seconds, 2),
  scoring_minutes = round(total_seconds / 60, 2),
  package_version = as.character(packageVersion("AUCell")),
  R_version = R.version.string,
  installation_difficulty = NA_character_,
  installation_seconds = NA_real_,
  peak_ram_gib = NA_real_,
  peak_memory_footprint_gib = NA_real_,
  total_execution_seconds = NA_real_,
  timestamp = as.character(Sys.time())
)

write.csv(
  benchmark,
  "results/benchmark_AUCell.csv",
  row.names = FALSE
)

# 9. Print summary
cat("\n========== AUCELL SCORING SUMMARY ==========\n")
cat("Cells scored:", nrow(scores), "\n")
cat("Modules scored:", ncol(scores) - 1, "\n")
cat("Missing scores:", sum(is.na(score_matrix)), "\n")
cat("Ranking time (seconds):", round(ranking_seconds, 2), "\n")
cat("AUC time (seconds):", round(auc_seconds, 2), "\n")
cat("Total scoring time (seconds):", round(total_seconds, 2), "\n")
cat("Total scoring time (minutes):", round(total_seconds / 60, 2), "\n")
cat("AUCell version:", as.character(packageVersion("AUCell")), "\n")

cat("\nSaved files:\n")
cat("results/WuImmune_AUCell_scores.csv\n")
cat("results/WuImmune_AUCell_scores.rds\n")
cat("data/WuImmune_AUCell_scored.rds\n")
cat("results/benchmark_AUCell.csv\n")

cat("\nScript 13 completed successfully!\n")
