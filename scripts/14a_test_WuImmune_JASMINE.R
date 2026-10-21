
# ============================================================
# Script 14A: Small JASMINE compatibility test
# ============================================================

library(Seurat)

source(
  "external/JASMINE/JASMINE_V1_11October2021.r"
)

# Load dataset and modules
wu <- readRDS("data/WuImmune_reclustered.rds")

modules <- readRDS(
  "data/WuImmune_62_modules_available.rds"
)

DefaultAssay(wu) <- "RNA"

# Obtain normalized expression matrix
expr <- SeuratObject::GetAssayData(
  wu,
  assay = "RNA",
  layer = "data"
)

# Select one module with at least two genes
valid_modules <- which(lengths(modules) >= 2)

stopifnot(length(valid_modules) > 0)

module_name <- names(modules)[valid_modules[1]]
module_genes <- modules[[valid_modules[1]]]

# Use only 20 cells for initial compatibility test
test_cells <- colnames(expr)[1:20]

test_expr <- as.matrix(
  expr[, test_cells, drop = FALSE]
)

cat("Testing module:", module_name, "\n")
cat("Module genes:", length(module_genes), "\n")
cat("Test cells:", ncol(test_expr), "\n")

# Run original JASMINE implementation
start_time <- Sys.time()

test_scores <- JASMINE(
  data = test_expr,
  genes = module_genes,
  method = "oddsratio"
)

elapsed <- as.numeric(
  difftime(Sys.time(), start_time, units = "secs")
)

cat("\nJASMINE output structure:\n")
str(test_scores)

cat("\nFirst few results:\n")
print(head(test_scores))

cat("\nTest execution time:", round(elapsed, 2), "seconds\n")

stopifnot(
  is.data.frame(test_scores),
  nrow(test_scores) == length(test_cells),
  "SampleID" %in% colnames(test_scores)
)

cat("\nJASMINE small test completed successfully!\n")

