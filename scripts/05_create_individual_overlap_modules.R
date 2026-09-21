# ============================================================
# 05_create_individual_overlap_modules.R
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
# First run:
#
#   Rscript 03_check_missing_module_genes.R
#
# Script 03 creates the latest dataset-compatible module list:
#
#   ../data/module_list_clean.rds
#
# Script 04 can be run before this script to inspect the
# overlapping module pairs and shared genes:
#
#   Rscript 04_check_module_overlap.R
#
# However, this script will independently calculate the
# overlapping genes from module_list_clean.rds.
#
# ------------------------------------------------------------
# PURPOSE:
# ------------------------------------------------------------
#
# Create a NEW experimental module library in which:
#
# 1. Genes occurring in more than one original module are
#    identified automatically.
#
# 2. These overlapping genes are removed from ALL of the
#    original modules.
#
# 3. Each overlapping gene will later become its own
#    individual one-gene module.
#
# IMPORTANT:
# The original module_list_clean.rds will NOT be modified.
#
# ============================================================


# ------------------------------------------------------------
# 1. Load the latest cleaned module library
# ------------------------------------------------------------

module_list_clean <- readRDS(
  "../data/module_list_clean.rds"
)

cat("\nNumber of cleaned modules:",
    length(module_list_clean), "\n")


# ------------------------------------------------------------
# 2. Automatically identify overlapping genes
# ------------------------------------------------------------

# Combine genes from all 62 modules.
# A gene can appear multiple times here if it belongs
# to multiple modules.

all_module_genes <- unlist(
  module_list_clean,
  use.names = FALSE
)

# Count how many times each gene occurs across modules
gene_counts <- table(all_module_genes)

# Keep genes that occur more than once
overlapping_genes <- names(
  gene_counts[gene_counts > 1]
)

# Sort alphabetically for easier inspection
overlapping_genes <- sort(overlapping_genes)


# ------------------------------------------------------------
# 3. Check the result
# ------------------------------------------------------------

cat("\n========================================\n")
cat("UNIQUE OVERLAPPING GENES\n")
cat("========================================\n\n")

print(overlapping_genes)

cat("\n========================================\n")
cat("TOTAL UNIQUE OVERLAPPING GENES\n")
cat("========================================\n")

print(length(overlapping_genes))

# ------------------------------------------------------------
# 4. Remove overlapping genes from the original modules
# ------------------------------------------------------------

# For every module:
# keep only genes that are NOT in the overlapping_genes list.
#
# IMPORTANT:
# module_list_clean remains unchanged.
# We create a NEW object called module_list_nonoverlap.

module_list_nonoverlap <- lapply(
  module_list_clean,
  function(module_genes) {

    setdiff(
      module_genes,
      overlapping_genes
    )

  }
)


# ------------------------------------------------------------
# 5. Count how many genes remain in each original module
# ------------------------------------------------------------

remaining_gene_counts <- lengths(
  module_list_nonoverlap
)

cat("\n========================================\n")
cat("GENES REMAINING AFTER OVERLAP REMOVAL\n")
cat("========================================\n\n")

print(remaining_gene_counts)


# ------------------------------------------------------------
# 6. Check whether any modules became empty
# ------------------------------------------------------------

empty_modules <- names(
  module_list_nonoverlap
)[
  remaining_gene_counts == 0
]


cat("\n========================================\n")
cat("MODULES THAT BECAME EMPTY\n")
cat("========================================\n\n")

if (length(empty_modules) == 0) {

  cat("No modules became empty.\n")

} else {

  print(empty_modules)

}


cat("\n========================================\n")
cat("NUMBER OF EMPTY MODULES\n")
cat("========================================\n")

print(length(empty_modules))


# ------------------------------------------------------------
# 7. Inspect empty modules and their ORIGINAL genes
# ------------------------------------------------------------

cat("\n========================================\n")
cat("EMPTY MODULES AND THEIR ORIGINAL GENES\n")
cat("========================================\n\n")

for (module_name in empty_modules) {

  cat("Module:", module_name, "\n")

  cat(
    "Original genes:",
    paste(
      module_list_clean[[module_name]],
      collapse = ", "
    ),
    "\n\n"
  )
}

# ------------------------------------------------------------
# 8. Create before-vs-after audit table
# ------------------------------------------------------------

module_audit <- data.frame(
  Module = names(module_list_clean),

  # Number of genes in each module BEFORE removing overlap
  Original_genes = lengths(module_list_clean),

  # Number of overlapping genes found in each module
  Overlapping_genes = sapply(
    module_list_clean,
    function(module_genes) {
      length(
        intersect(
          module_genes,
          overlapping_genes
        )
      )
    }
  ),

  # Number of genes removed from each module
  Removed_genes = sapply(
    module_list_clean,
    function(module_genes) {
      length(
        intersect(
          module_genes,
          overlapping_genes
        )
      )
    }
  ),

  # Number of genes remaining AFTER overlap removal
  Genes_left = lengths(module_list_nonoverlap),

  row.names = NULL
)


# ------------------------------------------------------------
# 9. Print the audit table
# ------------------------------------------------------------

cat("\n========================================\n")
cat("MODULE OVERLAP REMOVAL AUDIT\n")
cat("========================================\n\n")

print(
  module_audit,
  row.names = FALSE
)


# ------------------------------------------------------------
# 10. Validate the audit table
# ------------------------------------------------------------

# For every module, this should be TRUE:
#
# Original genes - Removed genes = Genes left

audit_check <- (
  module_audit$Original_genes -
  module_audit$Removed_genes
) == module_audit$Genes_left


cat("\n========================================\n")
cat("AUDIT VALIDATION\n")
cat("========================================\n")

if (all(audit_check)) {

  cat("PASS: All 62 modules are mathematically consistent.\n")

} else {

  cat("WARNING: Some modules failed the validation.\n")

  print(
    module_audit[!audit_check, ]
  )
}