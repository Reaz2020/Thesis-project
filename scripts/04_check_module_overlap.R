# ============================================================
# 04_check_module_overlap.R
#
# PROJECT:
# AI-Driven Gene Module Scoring for Cross-Dataset
# Single-Cell Immune Profiling
#
# DATASET:
# ICI_cRA_joint_small_for_module_pipeline.RData
#
# Seurat object:
# obj_small
#
# ORIGINAL MODULE LIBRARY:
# ICI_62_modules_annotation_references.xlsx
#
# INPUT USED BY THIS SCRIPT:
# ../data/module_list_clean.rds
#
# ------------------------------------------------------------
# BEFORE RUNNING THIS SCRIPT:
# ------------------------------------------------------------
#
# Run this script first:
#
#   Rscript 03_check_missing_module_genes.R
#
# Script 03:
#   1. Reads the latest module definitions from:
#      ICI_62_modules_annotation_references.xlsx
#
#   2. Loads:
#      ICI_cRA_joint_small_for_module_pipeline.RData
#
#   3. Checks every module gene against obj_small.
#
#   4. Removes genes that are not present in obj_small
#      from the WORKING COPY of the modules.
#
#   5. Saves the latest cleaned modules as:
#      ../data/module_list_clean.rds
#
# Therefore:
#
# If the Excel module definitions are changed,
# ALWAYS rerun script 03 before running this script.
#
# ------------------------------------------------------------
# PURPOSE OF THIS SCRIPT:
# ------------------------------------------------------------
#
# Find which of the 62 cleaned modules overlap with each other.
#
# DEFINITION OF OVERLAP:
# Two modules overlap if they contain at least one gene
# in common.
#
# IMPORTANT:
# This script only identifies overlapping module pairs.
# It does NOT modify, merge, remove, or create modules.
#
# ============================================================

# ------------------------------------------------------------
# 1. Load the latest cleaned module list
# ------------------------------------------------------------

module_list_clean <- readRDS(
  "../data/module_list_clean.rds"
)

cat("Dataset: ICI_cRA_joint_small_for_module_pipeline.RData\n")
cat("Cleaned module file: module_list_clean.rds\n")
cat("Number of cleaned modules:",
    length(module_list_clean), "\n")


# ------------------------------------------------------------
# 2. Compare every module with every other module
# ------------------------------------------------------------

# Get the names of all 62 modules
module_names <- names(module_list_clean)

# Create an empty table where we will store
# overlapping module pairs
overlap_results <- data.frame()

# Compare each module with every other module
for (i in 1:(length(module_names) - 1)) {

  for (j in (i + 1):length(module_names)) {

    # Names of the two modules being compared
    module_1 <- module_names[i]
    module_2 <- module_names[j]

    # Find genes that occur in BOTH modules
    shared_genes <- intersect(
      module_list_clean[[module_1]],
      module_list_clean[[module_2]]
    )

    # Only save the pair if at least one gene is shared
if (length(shared_genes) > 0) {

  overlap_results <- rbind(
    overlap_results,
    data.frame(
      Module_1 = module_1,
      Module_2 = module_2,

      # How many genes are shared
      Number_shared = length(shared_genes),

      # Names of the shared genes
      Shared_genes = paste(
        shared_genes,
        collapse = ", "
      )
    )
  )
}
  }
}


# ------------------------------------------------------------
# 3. Print the overlapping module pairs
# ------------------------------------------------------------

cat("\n========================================\n")
cat("OVERLAPPING MODULE PAIRS\n")
cat("========================================\n\n")

print(
  overlap_results,
  row.names = FALSE
)

cat("\nTotal overlapping module pairs:",
    nrow(overlap_results), "\n")    


    # ------------------------------------------------------------
# 4. Create a UNIQUE list of all overlapping genes
# ------------------------------------------------------------

# Take the Shared_genes column from all overlapping module pairs,
# split the comma-separated gene names,
# and combine them into one vector.

all_overlapping_genes <- unlist(
  strsplit(
    overlap_results$Shared_genes,
    ", "
  )
)

# Keep each gene only once
unique_overlapping_genes <- unique(
  all_overlapping_genes
)

# Sort alphabetically to make the list easier to read
unique_overlapping_genes <- sort(
  unique_overlapping_genes
)


# ------------------------------------------------------------
# 5. Print the unique overlapping genes
# ------------------------------------------------------------

cat("\n========================================\n")
cat("UNIQUE OVERLAPPING GENES\n")
cat("========================================\n\n")

print(unique_overlapping_genes)


# Print the total number of unique overlapping genes
cat("\n========================================\n")
cat("TOTAL UNIQUE OVERLAPPING GENES\n")
cat("========================================\n")

print(length(unique_overlapping_genes))