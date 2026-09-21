# ============================================================
# 03_check_missing_module_genes.R
#
# PURPOSE:
# Check whether genes from our 62 predefined modules
# are present in the scRNA-seq dataset.
#
# IMPORTANT:
# This script:
# 1. Identifies module genes that are missing from obj_small.
# 2. Removes those genes from a cleaned working copy.
# 3. Keeps the original Excel file unchanged.
# 4. Saves the cleaned modules as module_list_clean.rds.
# ============================================================


# ------------------------------------------------------------
# 1. Load required packages
# ------------------------------------------------------------

library(Seurat)
library(readxl)


# ------------------------------------------------------------
# 2. Load the Seurat object
# ------------------------------------------------------------

load("../data/ICI_cRA_joint_small_for_module_pipeline.RData")

# Check that the object loaded correctly
print(obj_small)


# ------------------------------------------------------------
# 3. Load the Excel file containing the 62 modules
# ------------------------------------------------------------

modules <- read_excel(
  "../data/ICI_62_modules_annotation_references.xlsx",
  sheet = 1
)

# Check number of modules
cat("\nNumber of modules:", nrow(modules), "\n")


# ------------------------------------------------------------
# 4. Get all genes that exist in obj_small
# ------------------------------------------------------------

dataset_genes <- rownames(obj_small)

cat("Number of genes in obj_small:",
    length(dataset_genes), "\n")


# ------------------------------------------------------------
# 5. Convert the Excel modules into an R list
# ------------------------------------------------------------

module_list <- setNames(
  
  lapply(modules$`Genes included`, function(x) {
    
    # Split genes at commas
    genes <- unlist(strsplit(x, ","))
    
    # Remove spaces before/after gene names
    genes <- trimws(genes)
    
    # Return the genes
    genes
  }),
  
  modules$`Module name`
)


# Check that we created 62 modules
cat("Number of modules in module_list:",
    length(module_list), "\n")


# ------------------------------------------------------------
# 6. Find genes that are NOT present in obj_small
# ------------------------------------------------------------

missing_genes <- lapply(
  
  module_list,
  
  function(module_genes) {
    
    # setdiff asks:
    #
    # Which genes are in the module
    # BUT NOT in dataset_genes?
    
    setdiff(module_genes, dataset_genes)
  }
)


# ------------------------------------------------------------
# 7. Keep only modules that have missing genes
# ------------------------------------------------------------

modules_with_missing_genes <-
  missing_genes[lengths(missing_genes) > 0]


# ------------------------------------------------------------
# 8. Print the results
# ------------------------------------------------------------

cat("\n========================================\n")
cat("MODULES WITH MISSING GENES\n")
cat("========================================\n\n")

print(modules_with_missing_genes)


# ------------------------------------------------------------
# 9. Count missing genes BEFORE cleaning
# ------------------------------------------------------------

missing_count <- lengths(missing_genes)

cat("\n========================================\n")
cat("TOTAL MISSING GENE ENTRIES\n")
cat("========================================\n")

print(sum(missing_count))


# ============================================================
# 10. REMOVE GENES THAT ARE NOT PRESENT IN THE DATASET
# ============================================================

# Create a NEW cleaned module list.
# The original module_list is NOT changed.

module_list_clean <- lapply(
  module_list,
  function(module_genes) {

    # Keep only genes that are present in obj_small
    intersect(module_genes, dataset_genes)

  }
)


# ============================================================
# 11. VERIFY THAT NO MISSING GENES REMAIN
# ============================================================

missing_after_cleaning <- lapply(
  module_list_clean,
  function(module_genes) {

    # Check whether any unavailable genes remain
    setdiff(module_genes, dataset_genes)

  }
)

# Keep only modules that still have missing genes
missing_after_cleaning <- missing_after_cleaning[
  lengths(missing_after_cleaning) > 0
]

cat("\n========================================\n")
cat("MISSING GENES AFTER CLEANING\n")
cat("========================================\n")

print(missing_after_cleaning)


# ============================================================
# 12. SAVE CLEANED MODULES
# ============================================================

saveRDS(
  module_list_clean,
  "../data/module_list_clean.rds"
)

cat("\n========================================\n")
cat("CLEANED MODULES SAVED SUCCESSFULLY\n")
cat("========================================\n")

cat("../data/module_list_clean.rds\n")