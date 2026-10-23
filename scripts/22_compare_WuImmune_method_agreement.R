
# ============================================================
# Script 22: Agreement between six scoring methods
#
# Dataset: WuImmune
# Cells: 50,063
# Gene modules: 62
# Methods: AddModuleScore, UCell, AUCell,
#          JASMINE, VISION, GSVA
#
# Analysis:
# 1. Load previously saved score files
# 2. Match cells using cell_id
# 3. Calculate Spearman correlations for each module
# 4. Summarize method-pair agreement
# 5. Generate comparison figures
#
# No gene-module scoring is repeated.
# ============================================================

cat("\n========================================\n")
cat("WUIMMUNE METHOD AGREEMENT ANALYSIS\n")
cat("========================================\n")

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

methods <- c(
  "AddModuleScore",
  "UCell",
  "AUCell",
  "JASMINE",
  "VISION",
  "GSVA"
)

# ------------------------------------------------------------
# 1. Load reference dataset and module definitions
# ------------------------------------------------------------

wu <- readRDS("data/WuImmune_reclustered.rds")

reference_cells <- colnames(wu)

modules <- readRDS(
  "data/WuImmune_62_modules_available.rds"
)

module_names <- names(modules)

stopifnot(
  length(reference_cells) == 50063,
  length(module_names) == 62,
  !anyDuplicated(reference_cells),
  !anyDuplicated(module_names)
)

cat("\nReference cells:", length(reference_cells), "\n")
cat("Gene modules:", length(module_names), "\n")
cat("Scoring methods:", length(methods), "\n")

# Free the large Seurat object.
rm(wu, modules)
invisible(gc())

# ------------------------------------------------------------
# 2. Load and align score matrices
# ------------------------------------------------------------

score_matrices <- list()

for (method in methods) {

  cat("\nLoading:", method, "\n")

  file_path <- paste0(
    "results/WuImmune_",
    method,
    "_scores.rds"
  )

  if (!file.exists(file_path)) {
    stop("Missing score file: ", file_path)
  }

  score_table <- readRDS(file_path)

  if (!is.data.frame(score_table)) {
    stop(method, ": expected a data.frame")
  }

  if (!"cell_id" %in% names(score_table)) {
    stop(method, ": missing cell_id column")
  }

  ids <- as.character(score_table$cell_id)

  if (anyNA(ids) || anyDuplicated(ids)) {
    stop(method, ": missing or duplicated cell IDs")
  }

  if (!setequal(ids, reference_cells)) {
    stop(method, ": cell IDs differ from WuImmune")
  }

  # Match by cell ID, not by original row position.
  cell_index <- match(reference_cells, ids)

  expected_columns <- paste0(
    method,
    "__",
    module_names
  )

  missing_columns <- setdiff(
    expected_columns,
    names(score_table)
  )

  if (length(missing_columns) > 0) {
    stop(
      method,
      ": missing module score columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  # Extract only the 62 module score columns.
  score_matrix <- as.matrix(
    score_table[
      cell_index,
      expected_columns,
      drop = FALSE
    ]
  )

  storage.mode(score_matrix) <- "double"

  colnames(score_matrix) <- module_names
  rownames(score_matrix) <- reference_cells

  if (any(!is.finite(score_matrix))) {
    stop(method, ": non-finite scores detected")
  }

  # A constant module cannot have Spearman correlation.
  constant_modules <- vapply(
    seq_len(ncol(score_matrix)),
    function(j) {
      length(unique(score_matrix[, j])) < 2
    },
    logical(1)
  )

  if (any(constant_modules)) {
    stop(
      method,
      ": constant-score modules: ",
      paste(
        colnames(score_matrix)[constant_modules],
        collapse = ", "
      )
    )
  }

  score_matrices[[method]] <- score_matrix

  cat(
    "Aligned dimensions:",
    nrow(score_matrix),
    "cells x",
    ncol(score_matrix),
    "modules\n"
  )

  rm(score_table, score_matrix)
  invisible(gc())
}

cat("\nAll six score files loaded and aligned.\n")

# ------------------------------------------------------------
# 3. Define the 15 unique method pairs
# ------------------------------------------------------------

method_pairs <- combn(
  methods,
  2,
  simplify = FALSE
)

cat("Method pairs:", length(method_pairs), "\n")

# ------------------------------------------------------------
# 4. Calculate Spearman rank correlations
#
# Spearman measures similarity of cell rankings,
# even when methods have different score scales.
# ------------------------------------------------------------

agreement_results <- vector(
  "list",
  length(module_names) * length(method_pairs)
)

result_index <- 1L

for (module in module_names) {

  cat("\nProcessing module:", module, "\n")

  # Rank scores once per method for this module.
  # This is more efficient than repeatedly ranking
  # the same scores for every method pair.
  ranked_scores <- lapply(
    methods,
    function(method) {
      rank(
        score_matrices[[method]][, module],
        ties.method = "average"
      )
    }
  )

  names(ranked_scores) <- methods

  for (pair in method_pairs) {

    method_1 <- pair[1]
    method_2 <- pair[2]

    # Pearson correlation of ranks equals
    # Spearman rank correlation.
    rho <- cor(
      ranked_scores[[method_1]],
      ranked_scores[[method_2]],
      method = "pearson"
    )

    agreement_results[[result_index]] <- data.frame(
      module = module,
      method_1 = method_1,
      method_2 = method_2,
      spearman_rho = as.numeric(rho),
      n_cells = length(reference_cells),
      stringsAsFactors = FALSE
    )

    result_index <- result_index + 1L
  }
}

agreement <- do.call(
  rbind,
  agreement_results
)

rownames(agreement) <- NULL

stopifnot(
  nrow(agreement) == 930,
  all(is.finite(agreement$spearman_rho))
)

# ------------------------------------------------------------
# 5. Save all 930 correlation results
# ------------------------------------------------------------

write.csv(
  agreement,
  "results/WuImmune_method_agreement_all_62_modules.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 6. Summarize agreement by method pair
# ------------------------------------------------------------

pair_summary <- do.call(
  rbind,
  lapply(method_pairs, function(pair) {

    subset_data <- agreement[
      agreement$method_1 == pair[1] &
      agreement$method_2 == pair[2],
      ,
      drop = FALSE
    ]

    data.frame(
      method_1 = pair[1],
      method_2 = pair[2],
      n_modules = nrow(subset_data),
      mean_spearman = mean(
        subset_data$spearman_rho
      ),
      median_spearman = median(
        subset_data$spearman_rho
      ),
      min_spearman = min(
        subset_data$spearman_rho
      ),
      max_spearman = max(
        subset_data$spearman_rho
      ),
      stringsAsFactors = FALSE
    )
  })
)

pair_summary <- pair_summary[
  order(-pair_summary$mean_spearman),
  ,
  drop = FALSE
]

rownames(pair_summary) <- NULL

write.csv(
  pair_summary,
  "results/WuImmune_method_agreement_pair_summary.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 7. Build the 6 x 6 method agreement matrix
# ------------------------------------------------------------

agreement_matrix <- matrix(
  NA_real_,
  nrow = length(methods),
  ncol = length(methods),
  dimnames = list(methods, methods)
)

diag(agreement_matrix) <- 1

for (i in seq_len(nrow(pair_summary))) {

  m1 <- pair_summary$method_1[i]
  m2 <- pair_summary$method_2[i]

  value <- pair_summary$mean_spearman[i]

  agreement_matrix[m1, m2] <- value
  agreement_matrix[m2, m1] <- value
}

write.csv(
  agreement_matrix,
  "results/WuImmune_method_agreement_matrix.csv"
)

# ------------------------------------------------------------
# 8. Create method-agreement heatmap
# ------------------------------------------------------------

png(
  "figures/WuImmune_method_agreement_heatmap.png",
  width = 1200,
  height = 1000,
  res = 150
)

par(mar = c(9, 10, 4, 3))

image(
  x = seq_along(methods),
  y = seq_along(methods),
  z = agreement_matrix,
  col = colorRampPalette(
    c("#B2182B", "#F7F7F7", "#2166AC")
  )(101),
  zlim = c(-1, 1),
  axes = FALSE,
  xlab = "",
  ylab = "",
  main = "WuImmune: Mean Spearman agreement\nacross 62 modules"
)

axis(
  1,
  at = seq_along(methods),
  labels = methods,
  las = 2
)

axis(
  2,
  at = seq_along(methods),
  labels = methods,
  las = 2
)

for (i in seq_along(methods)) {
  for (j in seq_along(methods)) {

    text(
      x = i,
      y = j,
      labels = sprintf(
        "%.2f",
        agreement_matrix[i, j]
      ),
      cex = 0.85
    )
  }
}

box()
dev.off()

# ------------------------------------------------------------
# 9. Create method-pair agreement bar chart
# ------------------------------------------------------------

pair_labels <- paste(
  pair_summary$method_1,
  "vs",
  pair_summary$method_2
)

png(
  "figures/WuImmune_method_pair_agreement.png",
  width = 1500,
  height = 1100,
  res = 150
)

par(mar = c(5, 19, 4, 3))

barplot(
  rev(pair_summary$mean_spearman),
  names.arg = rev(pair_labels),
  horiz = TRUE,
  las = 1,
  xlim = c(-1, 1),
  col = "#4477AA",
  xlab = "Mean Spearman correlation across 62 modules",
  main = "WuImmune: Method-pair agreement"
)

abline(v = 0, col = "gray40")

dev.off()

# ------------------------------------------------------------
# 10. Identify modules with lowest overall agreement
# ------------------------------------------------------------

module_summary <- aggregate(
  spearman_rho ~ module,
  data = agreement,
  FUN = mean
)

names(module_summary)[2] <- "mean_pairwise_spearman"

module_summary <- module_summary[
  order(module_summary$mean_pairwise_spearman),
  ,
  drop = FALSE
]

rownames(module_summary) <- NULL

write.csv(
  module_summary,
  "results/WuImmune_module_agreement_summary.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 11. Print final results
# ------------------------------------------------------------

cat("\n========================================\n")
cat("METHOD AGREEMENT SUMMARY\n")
cat("========================================\n")

print(
  pair_summary,
  row.names = FALSE,
  digits = 4
)

cat("\nHighest-agreement pair:\n")
print(
  pair_summary[1, ],
  row.names = FALSE
)

cat("\nLowest-agreement pair:\n")
print(
  pair_summary[nrow(pair_summary), ],
  row.names = FALSE
)

cat("\nModules with lowest average agreement:\n")

print(
  head(module_summary, 10),
  row.names = FALSE,
  digits = 4
)

cat("\nFiles saved:\n")
cat("results/WuImmune_method_agreement_all_62_modules.csv\n")
cat("results/WuImmune_method_agreement_pair_summary.csv\n")
cat("results/WuImmune_method_agreement_matrix.csv\n")
cat("results/WuImmune_module_agreement_summary.csv\n")
cat("figures/WuImmune_method_agreement_heatmap.png\n")
cat("figures/WuImmune_method_pair_agreement.png\n")

cat("\nIMPORTANT:\n")
cat("High agreement means similar cell rankings.\n")
cat("It does not establish biological accuracy.\n")

cat("\nScript 22 completed successfully!\n")
