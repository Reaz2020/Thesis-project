# ============================================================
# Script 20: Inspect WuImmune annotations and score files
# Preparation for biological comparison of six methods
# Does not modify the dataset or scoring results
# ============================================================

library(Seurat)

cat("========================================\n")
cat("WUIMMUNE BIOLOGICAL EVALUATION SETUP\n")
cat("========================================\n")

# 1. Load original dataset

wu <- readRDS("data/WuImmune_reclustered.rds")

cat("\nCells:", ncol(wu), "\n")
cat("Genes:", nrow(wu), "\n")

# 2. Inspect metadata columns

cat("\n========================================\n")
cat("AVAILABLE METADATA COLUMNS\n")
cat("========================================\n")

print(colnames(wu@meta.data))

# 3. Inspect existing Seurat cell identities

cat("\n========================================\n")
cat("CURRENT CELL IDENTITIES\n")
cat("========================================\n")

print(sort(table(Idents(wu)), decreasing = TRUE))

# 4. Identify potential annotation columns

annotation_columns <- grep(
  "celltype|cell_type|annotation|label|cluster|identity",
  colnames(wu@meta.data),
  ignore.case = TRUE,
  value = TRUE
)

cat("\n========================================\n")
cat("POTENTIAL ANNOTATION COLUMNS\n")
cat("========================================\n")

print(annotation_columns)

for (column in annotation_columns) {

  cat("\nColumn:", column, "\n")

  values <- wu@meta.data[[column]]

  cat("Number of unique values:",
      length(unique(values)), "\n")

  if (length(unique(values)) <= 50) {
    print(sort(table(values), decreasing = TRUE))
  }
}

# 5. Inspect module names

modules <- readRDS(
  "data/WuImmune_62_modules_available.rds"
)

cat("\n========================================\n")
cat("62 IMMUNE MODULE NAMES\n")
cat("========================================\n")

print(names(modules))

# 6. Locate saved score files

cat("\n========================================\n")
cat("AVAILABLE SCORE FILES\n")
cat("========================================\n")

score_files <- list.files(
  "results",
  pattern = "WuImmune.*scores.*\\.(csv|rds)$",
  full.names = TRUE
)

print(score_files)

cat("\n========================================\n")
cat("INSPECTION COMPLETED\n")
cat("========================================\n")
