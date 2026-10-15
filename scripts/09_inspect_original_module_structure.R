
# ============================================================
# Script 09: Inspect original module structure
# No files are modified or overwritten.
# ============================================================

library(readxl)

# 1. Read original curated modules
module_table <- read_excel(
  "data/ICI_62_modules_annotation_references.xlsx"
)

cat("\n=== ORIGINAL EXCEL STRUCTURE ===\n")
print(dim(module_table))
print(names(module_table))
str(module_table)

cat("\n=== FIRST THREE MODULES ===\n")
print(module_table[1:3, ], width = Inf)

cat("\n=== ALL MODULE NAMES AND DESCRIPTIONS ===\n")
print(
  module_table[, c("Module name", "Description")],
  n = Inf,
  width = Inf
)

# 2. Inspect the WuImmune-compatible module list
modules <- readRDS(
  "data/WuImmune_62_modules_available.rds"
)

cat("\n=== CLEANED MODULE LIST STRUCTURE ===\n")
str(modules, max.level = 1)

cat("\n=== FIRST THREE MODULE GENE VECTORS ===\n")
print(modules[1:3])

# 3. Check names and structure
stopifnot(
  nrow(module_table) == 62,
  length(modules) == 62,
  !anyDuplicated(module_table[["Module name"]]),
  setequal(module_table[["Module name"]], names(modules)),
  all(vapply(modules, is.character, logical(1)))
)

cat("\nStructure checks passed.\n")
