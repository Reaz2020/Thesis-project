
# ============================================================
# Script 08: Check original 62 module genes against WuImmune
# ============================================================

library(Seurat)
library(readxl)

# 1. Load the reclustered WuImmune object
WuImmune <- readRDS("data/WuImmune_reclustered.rds")

# 2. Read the ORIGINAL module library
module_table <- read_excel(
  "data/ICI_62_modules_annotation_references.xlsx"
)

# Inspect the column names
print(colnames(module_table))

# 3. Convert the Excel gene lists into an R named list
# Each row contains a module name and its genes.

module_list <- lapply(
  module_table[["Genes included"]],
  function(x) {
    if (is.na(x) || !nzchar(trimws(x))) {
      return(character(0))
    }

    # Split genes separated by commas, semicolons,
    # or whitespace
    genes <- unlist(strsplit(x, "[,;[:space:]]+"))

    # Remove empty entries and duplicates within a module
    unique(genes[nzchar(genes)])
  }
)

names(module_list) <- module_table[["Module name"]]

# 4. Get all gene names available in WuImmune
wu_genes <- rownames(WuImmune[["RNA"]])

# 5. Check present and missing genes for every module
module_report <- do.call(
  rbind,
  lapply(names(module_list), function(module_name) {

    genes <- module_list[[module_name]]

    present <- intersect(genes, wu_genes)
    missing <- setdiff(genes, wu_genes)

    data.frame(
      Module = module_name,
      Total_genes = length(genes),
      Present_genes = length(present),
      Missing_genes = length(missing),
      Coverage_percent = if (length(genes) > 0) {
        round(100 * length(present) / length(genes), 1)
      } else {
        NA_real_
      },
      Missing_gene_names = paste(missing, collapse = ", "),
      stringsAsFactors = FALSE
    )
  })
)

# 6. Print the summary
print(module_report)

cat("\nNumber of modules:", length(module_list), "\n")
cat(
  "Modules with missing genes:",
  sum(module_report$Missing_genes > 0),
  "\n"
)
cat(
  "Modules with no available genes:",
  sum(module_report$Present_genes == 0),
  "\n"
)

# 7. Create a compatible gene list for WuImmune
# Keep only genes present in this dataset.
module_list_WuImmune <- lapply(
  module_list,
  function(genes) intersect(genes, wu_genes)
)

# 8. Save the results
dir.create("results", showWarnings = FALSE)

write.csv(
  module_report,
  "results/WuImmune_module_gene_coverage.csv",
  row.names = FALSE
)

saveRDS(
  module_list_WuImmune,
  "data/WuImmune_62_modules_available.rds"
)

cat("\nResults saved successfully.\n")
