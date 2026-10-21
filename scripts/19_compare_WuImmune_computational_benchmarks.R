# ============================================================
# Script 19: Compare six WuImmune scoring methods
# Handles different benchmark CSV formats.
# Does NOT modify original benchmark files.
# ============================================================

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

methods <- c(
  "AddModuleScore", "UCell", "AUCell",
  "JASMINE", "VISION", "GSVA"
)

# Known external measurements from completed runs.
# These are used only when the original CSV lacks them.
external <- data.frame(
  method = methods,
  external_elapsed_seconds = c(
    25.18, 269.18, 170.73, 3266.37, 21.38, 9981.11
  ),
  max_resident_bytes = c(
    2018951168, 1990639616, 2051719168,
    2001649664, 2362179584, 2603712512
  ),
  peak_footprint_bytes = c(
    4329434432, 6882411968, 5427327168,
    5025691072, 7466237312, 8922401920
  ),
  stringsAsFactors = FALSE
)

get_value <- function(x, column, default = NA) {
  if (!column %in% names(x)) return(default)
  value <- x[[column]][1]
  if (length(value) == 0 || is.na(value)) return(default)
  value
}

read_method <- function(method) {

  path <- file.path(
    "results", paste0("benchmark_", method, ".csv")
  )

  if (!file.exists(path)) stop("Missing file: ", path)

  x <- suppressWarnings(
    read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  )

  if (nrow(x) != 1) {
    stop("Expected one benchmark row for ", method)
  }

  if (
    "method" %in% names(x) &&
    as.character(x$method[1]) != method
  ) {
    stop("Method mismatch for ", method)
  }

  ext <- external[external$method == method, ]

  resident <- get_value(
    x, "max_resident_bytes", ext$max_resident_bytes
  )

  footprint <- get_value(
    x, "peak_footprint_bytes", ext$peak_footprint_bytes
  )

  external_seconds <- get_value(
    x,
    "external_elapsed_seconds",
    ext$external_elapsed_seconds
  )

  # Prefer precise byte measurements over rounded CSV GiB.
  peak_ram <- round(as.numeric(resident) / 1024^3, 2)
  peak_footprint <- round(
    as.numeric(footprint) / 1024^3, 2
  )

  data.frame(
    method = method,
    package_version = as.character(
      get_value(x, "package_version", NA_character_)
    ),
    cells = as.numeric(get_value(x, "cells")),
    modules = as.numeric(get_value(x, "modules")),
    genes = as.numeric(get_value(x, "genes", 29733)),
    scoring_seconds = as.numeric(
      get_value(x, "scoring_seconds")
    ),
    total_execution_seconds = as.numeric(
      get_value(x, "total_execution_seconds")
    ),
    external_elapsed_seconds = as.numeric(
      external_seconds
    ),
    peak_ram_gib = peak_ram,
    peak_memory_footprint_gib = peak_footprint,
    max_resident_bytes = as.numeric(resident),
    peak_footprint_bytes = as.numeric(footprint),
    installation_difficulty = as.character(
      get_value(x, "installation_difficulty", NA_character_)
    ),
    missing_scores = as.numeric(
      get_value(x, "missing_scores")
    ),
    nonfinite_scores = as.numeric(
      get_value(x, "nonfinite_scores")
    ),
    status = as.character(
      get_value(x, "status", NA_character_)
    ),
    stringsAsFactors = FALSE
  )
}

cat("========================================\n")
cat("COMPUTATIONAL BENCHMARK COMPARISON\n")
cat("========================================\n")

comparison <- do.call(
  rbind,
  lapply(methods, read_method)
)

# Validate measurements that exist for every method.
stopifnot(
  nrow(comparison) == 6,
  all(comparison$cells == 50063),
  all(comparison$modules == 62),
  all(comparison$genes == 29733),
  all(is.finite(comparison$scoring_seconds)),
  all(comparison$scoring_seconds > 0),
  all(is.finite(comparison$peak_ram_gib)),
  all(is.finite(comparison$peak_memory_footprint_gib))
)

comparison$scoring_minutes <- round(
  comparison$scoring_seconds / 60, 2
)

comparison$external_elapsed_minutes <- round(
  comparison$external_elapsed_seconds / 60, 2
)

comparison$scores_generated <-
  comparison$cells * comparison$modules

comparison <- comparison[
  order(comparison$scoring_seconds),
]

rownames(comparison) <- NULL

# ------------------------------------------------------------
# Save full and summary tables
# ------------------------------------------------------------

write.csv(
  comparison,
  "results/WuImmune_computational_benchmark_comparison.csv",
  row.names = FALSE
)

summary_table <- comparison[, c(
  "method",
  "package_version",
  "scoring_seconds",
  "scoring_minutes",
  "external_elapsed_seconds",
  "peak_ram_gib",
  "peak_memory_footprint_gib",
  "installation_difficulty"
)]

write.csv(
  summary_table,
  "results/WuImmune_computational_summary_table.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# Plot helper
# ------------------------------------------------------------

make_plot <- function(values, filename, title, xlabel, color) {

  png(
    filename,
    width = 1800,
    height = 1100,
    res = 180
  )

  par(mar = c(5, 13, 4, 3))

  barplot(
    values,
    names.arg = comparison$method,
    horiz = TRUE,
    las = 1,
    col = color,
    border = NA,
    main = title,
    xlab = xlabel,
    xlim = c(0, max(values) * 1.12)
  )

  dev.off()
}

# ------------------------------------------------------------
# Generate three figures
# ------------------------------------------------------------

make_plot(
  comparison$scoring_minutes,
  "figures/WuImmune_scoring_runtime_comparison.png",
  "WuImmune: Scoring Runtime",
  "Scoring time (minutes)",
  "steelblue"
)

make_plot(
  comparison$peak_ram_gib,
  "figures/WuImmune_peak_RAM_comparison.png",
  "WuImmune: Peak Resident Memory",
  "Peak RAM (GiB)",
  "darkseagreen4"
)

make_plot(
  comparison$peak_memory_footprint_gib,
  "figures/WuImmune_peak_memory_footprint_comparison.png",
  "WuImmune: Peak Memory Footprint",
  "Peak memory footprint (GiB)",
  "darkorange"
)

# ------------------------------------------------------------
# Print final summary
# ------------------------------------------------------------

cat("\n========================================\n")
cat("COMPARISON SUMMARY\n")
cat("========================================\n")

print(summary_table, row.names = FALSE)

cat(
  "\nFastest method:",
  comparison$method[which.min(comparison$scoring_seconds)],
  "\n"
)

cat(
  "Lowest peak RAM:",
  comparison$method[which.min(comparison$peak_ram_gib)],
  "\n"
)

cat(
  "Lowest memory footprint:",
  comparison$method[
    which.min(comparison$peak_memory_footprint_gib)
  ],
  "\n"
)

cat("\nOutput files:\n")
cat("results/WuImmune_computational_benchmark_comparison.csv\n")
cat("results/WuImmune_computational_summary_table.csv\n")
cat("figures/WuImmune_scoring_runtime_comparison.png\n")
cat("figures/WuImmune_peak_RAM_comparison.png\n")
cat("figures/WuImmune_peak_memory_footprint_comparison.png\n")

cat("\nScript 19 completed successfully!\n")