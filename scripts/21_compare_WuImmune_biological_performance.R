
# ============================================================
# Script 21: Biological comparison of six scoring methods
# Dataset: WuImmune
#
# Evaluation:
#   - Cell-ID alignment
#   - AUROC against existing cell annotations
#   - Method-level AUROC summaries
#   - AUROC heatmap and comparison plots
#
# No rescoring is performed.
# ============================================================

cat("\n========================================\n")
cat("WUIMMUNE BIOLOGICAL METHOD COMPARISON\n")
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
# 1. Load WuImmune reference annotations
# ------------------------------------------------------------

wu <- readRDS("data/WuImmune_reclustered.rds")

cell_ids <- colnames(wu)
metadata <- wu@meta.data

stopifnot(
  length(cell_ids) == nrow(metadata),
  !anyDuplicated(cell_ids),
  all(cell_ids %in% rownames(metadata))
)

metadata <- metadata[cell_ids, , drop = FALSE]

cat("\nReference cells:", length(cell_ids), "\n")

# ------------------------------------------------------------
# 2. Define candidate biological validation tasks
#
# Each task specifies:
# - a module
# - the annotation column
# - positive reference labels
# - negative reference labels
#
# These labels are dataset-derived reference annotations,
# not independent experimental ground truth.
# ------------------------------------------------------------

tasks <- list(
  list(
    module = "B_cell_identity_module",
    annotation = "celltype_minor",
    positive = c("B cells Memory", "B cells Naive"),
    negative = NULL
  ),
  list(
    module = "CD4_identity_module",
    annotation = "celltype_minor",
    positive = "T cells CD4+",
    negative = c("T cells CD8+")
  ),
  list(
    module = "CD8_identity_module",
    annotation = "celltype_minor",
    positive = "T cells CD8+",
    negative = c("T cells CD4+")
  ),
  list(
    module = "Memory_B_cell_module",
    annotation = "celltype_minor",
    positive = "B cells Memory",
    negative = "B cells Naive"
  ),
  list(
    module = "Naive_B_cell_module",
    annotation = "celltype_minor",
    positive = "B cells Naive",
    negative = "B cells Memory"
  ),
  list(
    module = "cDC1_identity_module",
    annotation = "celltype_subset",
    positive = "Myeloid_c3_cDC1_CLEC9A",
    negative = c(
      "Myeloid_c11_cDC2_CD1C",
      "Myeloid_c4_DCs_pDC_IRF7",
      "Myeloid_c0_DC_LAMP3"
    )
  ),
  list(
    module = "cDC2_identity_module",
    annotation = "celltype_subset",
    positive = "Myeloid_c11_cDC2_CD1C",
    negative = c(
      "Myeloid_c3_cDC1_CLEC9A",
      "Myeloid_c4_DCs_pDC_IRF7",
      "Myeloid_c0_DC_LAMP3"
    )
  ),
  list(
    module = "pDC_identity_module",
    annotation = "celltype_subset",
    positive = "Myeloid_c4_DCs_pDC_IRF7",
    negative = c(
      "Myeloid_c3_cDC1_CLEC9A",
      "Myeloid_c11_cDC2_CD1C",
      "Myeloid_c0_DC_LAMP3"
    )
  ),
  list(
    module = "LAMP3_DC_module",
    annotation = "celltype_subset",
    positive = "Myeloid_c0_DC_LAMP3",
    negative = c(
      "Myeloid_c3_cDC1_CLEC9A",
      "Myeloid_c11_cDC2_CD1C",
      "Myeloid_c4_DCs_pDC_IRF7"
    )
  ),
  list(
    module = "Treg_module",
    annotation = "celltype_subset",
    positive = "T_cells_c2_CD4+_T-regs_FOXP3",
    negative = c(
      "T_cells_c0_CD4+_CCR7",
      "T_cells_c1_CD4+_IL7R",
      "T_cells_c3_CD4+_Tfh_CXCL13"
    )
  )
)

# ------------------------------------------------------------
# 3. Validate task labels before calculating AUROC
# ------------------------------------------------------------

for (task in tasks) {

  stopifnot(task$annotation %in% colnames(metadata))

  available <- unique(
    as.character(metadata[[task$annotation]])
  )

  required <- c(task$positive, task$negative)

  missing_labels <- setdiff(required, available)

  if (length(missing_labels) > 0) {
    stop(
      "Missing reference labels for ",
      task$module, ": ",
      paste(missing_labels, collapse = ", ")
    )
  }
}

# ------------------------------------------------------------
# 4. AUROC calculation using ranks
#
# AUROC = probability that a randomly selected positive
# cell receives a higher score than a negative cell,
# with ties contributing 0.5.
# ------------------------------------------------------------

calculate_auc <- function(scores, labels) {

  stopifnot(
    length(scores) == length(labels),
    all(is.finite(scores)),
    all(!is.na(labels))
  )

  n_positive <- sum(labels)
  n_negative <- sum(!labels)

  if (n_positive == 0 || n_negative == 0) {
    return(NA_real_)
  }

  score_ranks <- rank(scores, ties.method = "average")

  auc <- (
    sum(score_ranks[labels]) -
      n_positive * (n_positive + 1) / 2
  ) / (n_positive * n_negative)

  as.numeric(auc)
}

# ------------------------------------------------------------
# 5. Read and align the six saved score tables
# ------------------------------------------------------------

all_results <- list()

for (method in methods) {

  cat("\n----------------------------------------\n")
  cat("Evaluating:", method, "\n")
  cat("----------------------------------------\n")

  path <- paste0(
    "results/WuImmune_",
    method,
    "_scores.rds"
  )

  if (!file.exists(path)) {
    stop("Missing score file: ", path)
  }

  scores <- readRDS(path)

  stopifnot(
    is.data.frame(scores),
    "cell_id" %in% colnames(scores),
    !anyNA(scores$cell_id),
    !anyDuplicated(scores$cell_id)
  )

  if (
    length(scores$cell_id) != length(cell_ids) ||
    !setequal(scores$cell_id, cell_ids)
  ) {
    stop("Cell IDs do not match WuImmune for ", method)
  }

  # Explicit alignment by cell ID.
  index <- match(cell_ids, scores$cell_id)

  if (anyNA(index)) {
    stop("Failed to align cell IDs for ", method)
  }

  scores <- scores[index, , drop = FALSE]

  stopifnot(identical(
    as.character(scores$cell_id),
    as.character(cell_ids)
  ))

  method_results <- list()

  for (i in seq_along(tasks)) {

    task <- tasks[[i]]

    score_column <- paste0(
      method,
      "__",
      task$module
    )

    if (!score_column %in% colnames(scores)) {
      stop("Missing score column: ", score_column)
    }

    annotations <- as.character(
      metadata[[task$annotation]]
    )

    positive_cells <- annotations %in% task$positive

    if (is.null(task$negative)) {

      # One-versus-all for broad B-cell identity.
      selected_cells <- !is.na(annotations)

    } else {

      # Specific biological contrast.
      selected_cells <- annotations %in% c(
        task$positive,
        task$negative
      )
    }

    labels <- positive_cells[selected_cells]

    values <- scores[[score_column]][selected_cells]

    if (!is.numeric(values)) {
      stop("Non-numeric scores in ", score_column)
    }

    if (any(!is.finite(values))) {
      stop("Non-finite scores in ", score_column)
    }

    auc <- calculate_auc(values, labels)

    method_results[[i]] <- data.frame(
      method = method,
      module = task$module,
      annotation = task$annotation,
      positive_labels = paste(
        task$positive,
        collapse = "; "
      ),
      negative_labels = if (is.null(task$negative)) {
        "All other annotated cells"
      } else {
        paste(task$negative, collapse = "; ")
      },
      n_positive = sum(labels),
      n_negative = sum(!labels),
      auroc = auc,
      stringsAsFactors = FALSE
    )

    cat(
      sprintf(
        "%-32s AUROC = %.4f (positive=%d, negative=%d)\n",
        task$module,
        auc,
        sum(labels),
        sum(!labels)
      )
    )
  }

  all_results[[method]] <- do.call(
    rbind,
    method_results
  )

  rm(scores)
  invisible(gc())
}

# ------------------------------------------------------------
# 6. Save detailed results
# ------------------------------------------------------------

results <- do.call(rbind, all_results)
rownames(results) <- NULL

detailed_path <- paste0(
  "results/",
  "WuImmune_biological_AUROC_comparison.csv"
)

write.csv(
  results,
  detailed_path,
  row.names = FALSE
)

# ------------------------------------------------------------
# 7. Method-level summaries
# ------------------------------------------------------------

summary_table <- do.call(
  rbind,
  lapply(methods, function(method) {

    values <- results$auroc[
      results$method == method
    ]

    data.frame(
      method = method,
      evaluated_modules = sum(!is.na(values)),
      mean_auroc = mean(values, na.rm = TRUE),
      median_auroc = median(values, na.rm = TRUE),
      min_auroc = min(values, na.rm = TRUE),
      max_auroc = max(values, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
)

summary_table <- summary_table[
  order(-summary_table$mean_auroc),
  ,
  drop = FALSE
]

rownames(summary_table) <- NULL

write.csv(
  summary_table,
  "results/WuImmune_biological_AUROC_summary.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 8. Build AUROC matrix for figures
# ------------------------------------------------------------

module_names <- vapply(
  tasks,
  function(task) task$module,
  character(1)
)

auc_matrix <- matrix(
  NA_real_,
  nrow = length(module_names),
  ncol = length(methods),
  dimnames = list(module_names, methods)
)

for (i in seq_len(nrow(results))) {

  auc_matrix[
    results$module[i],
    results$method[i]
  ] <- results$auroc[i]
}

write.csv(
  auc_matrix,
  "results/WuImmune_biological_AUROC_matrix.csv"
)

# ------------------------------------------------------------
# 9. AUROC heatmap
# ------------------------------------------------------------

png(
  "figures/WuImmune_biological_AUROC_heatmap.png",
  width = 1800,
  height = 1250,
  res = 160
)

par(
  mar = c(8, 17, 4, 3)
)

image(
  x = seq_len(ncol(auc_matrix)),
  y = seq_len(nrow(auc_matrix)),
  z = t(auc_matrix),
  col = colorRampPalette(
    c("#B2182B", "#F7F7F7", "#2166AC")
  )(101),
  zlim = c(0, 1),
  axes = FALSE,
  xlab = "",
  ylab = "",
  main = "WuImmune: Biological AUROC comparison"
)

axis(
  1,
  at = seq_along(methods),
  labels = methods,
  las = 2
)

axis(
  2,
  at = seq_along(module_names),
  labels = module_names,
  las = 2,
  cex.axis = 0.72
)

for (i in seq_len(nrow(auc_matrix))) {
  for (j in seq_len(ncol(auc_matrix))) {

    text(
      x = j,
      y = i,
      labels = sprintf("%.2f", auc_matrix[i, j]),
      cex = 0.68
    )
  }
}

box()
dev.off()

# ------------------------------------------------------------
# 10. Mean AUROC comparison
# ------------------------------------------------------------

png(
  "figures/WuImmune_biological_mean_AUROC.png",
  width = 1300,
  height = 900,
  res = 150
)

par(mar = c(6, 11, 4, 3))

barplot(
  rev(summary_table$mean_auroc),
  names.arg = rev(summary_table$method),
  horiz = TRUE,
  las = 1,
  xlim = c(0, 1),
  col = "#4477AA",
  xlab = "Mean AUROC across selected identity tasks",
  main = "WuImmune: Biological discrimination by method"
)

abline(v = 0.5, lty = 2, col = "gray40")

dev.off()

# ------------------------------------------------------------
# 11. Final summary
# ------------------------------------------------------------

cat("\n========================================\n")
cat("BIOLOGICAL COMPARISON SUMMARY\n")
cat("========================================\n")

print(
  summary_table,
  row.names = FALSE,
  digits = 4
)

cat("\nHighest mean AUROC in these tasks:",
    summary_table$method[1], "\n")

cat("\nFiles saved:\n")
cat("results/WuImmune_biological_AUROC_comparison.csv\n")
cat("results/WuImmune_biological_AUROC_summary.csv\n")
cat("results/WuImmune_biological_AUROC_matrix.csv\n")
cat("figures/WuImmune_biological_AUROC_heatmap.png\n")
cat("figures/WuImmune_biological_mean_AUROC.png\n")

cat("\nIMPORTANT:\n")
cat("These AUROCs measure agreement with existing annotations.\n")
cat("They do not establish independent biological ground truth.\n")
cat("The mean AUROC is descriptive, not a definitive ranking.\n")

cat("\nScript 21 completed successfully!\n")
