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
  "data/module_list_clean.rds"
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



# ============================================================
# STEP: Create one individual module for each overlapping gene
# ============================================================

# We already identified 102 genes that appeared in more than
# one of the original modules.
#
# Here, each overlapping gene becomes its OWN module.
#
# Example:
# AICDA  -> AICDA_individual_module = "AICDA"
# CTLA4  -> CTLA4_individual_module = "CTLA4"

individual_overlap_modules <- lapply(
  overlapping_genes,
  function(gene) {
    gene
  }
)

# Give each one-gene module a unique module name.
names(individual_overlap_modules) <- paste0(
  overlapping_genes,
  "_individual_module"
)


# ============================================================
# CHECK 1: How many individual modules were created?
# ============================================================

cat("\n========================================\n")
cat("INDIVIDUAL OVERLAP MODULES CREATED\n")
cat("========================================\n")

cat("Number of individual modules:",
    length(individual_overlap_modules), "\n")


# ============================================================
# CHECK 2: Does every module contain exactly one gene?
# ============================================================

individual_module_sizes <- lengths(individual_overlap_modules)

if (all(individual_module_sizes == 1)) {
  cat("PASS: Every individual module contains exactly one gene.\n")
} else {
  cat("WARNING: At least one individual module does not contain exactly one gene.\n")
}


# ============================================================
# CHECK 3: Print the individual modules
# ============================================================

cat("\n========================================\n")
cat("INDIVIDUAL MODULE LIST\n")
cat("========================================\n")

print(individual_overlap_modules)


# ============================================================
# SAVE THE 102 INDIVIDUAL OVERLAP MODULES
# ============================================================

# Save the list containing the 102 one-gene modules.
# This does NOT modify module_list_clean.rds.
saveRDS(
  individual_overlap_modules,
  file = "data/individual_overlap_modules.rds"
)

cat("\n========================================\n")
cat("INDIVIDUAL MODULES SAVED\n")
cat("========================================\n")
cat("Saved 102 individual modules to:\n")
cat("data/individual_overlap_modules.rds\n")


# ============================================================
# VERIFY THE SAVED 102 INDIVIDUAL MODULES
# ============================================================

# Read the file back from disk.
# This checks what was actually saved, not just what was in memory.
individual_modules_check <- readRDS(
  "data/individual_overlap_modules.rds"
)

cat("\n========================================\n")
cat("VERIFY SAVED INDIVIDUAL MODULES\n")
cat("========================================\n")

# Number of modules
cat("Total modules:",
    length(individual_modules_check), "\n")

# Check how many genes each module contains
cat("\nModule sizes:\n")
print(table(lengths(individual_modules_check)))

# Show all module names
cat("\nModule names:\n")
print(names(individual_modules_check))

# Show complete module -> gene mapping
cat("\nComplete module contents:\n")
print(individual_modules_check)


# ============================================================
# CHECK SCHEMA OF module_list_clean.rds
# ============================================================

module_list_schema_check <- readRDS(
  "data/module_list_clean.rds"
)

cat("\n========================================\n")
cat("SCHEMA OF module_list_clean.rds\n")
cat("========================================\n")

# What kind of R object is it?
cat("\nClass:\n")
print(class(module_list_schema_check))

# Overall internal structure
cat("\nStructure:\n")
str(module_list_schema_check)

# Number of modules
cat("\nNumber of modules:\n")
print(length(module_list_schema_check))

# Show first few module names
cat("\nFirst 10 module names:\n")
print(head(names(module_list_schema_check), 10))

# Show first few modules and their genes
cat("\nFirst 3 modules and contents:\n")
print(module_list_schema_check[1:3])