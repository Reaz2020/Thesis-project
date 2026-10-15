
# ============================================================
# Script 09: Separate original 62 modules into two Excel files
# Does NOT modify the original workbook or gene lists
# ============================================================

library(readxl)
library(writexl)

# 1. Load original module library
original <- read_excel(
  "data/ICI_62_modules_annotation_references.xlsx"
)

# 2. Identity modules: proposed classification
identity_names <- c(
  "Age_associated_B_module",
  "B_cell_identity_module",
  "C1Q_module",
  "CD3_module",
  "CD4_identity_module",
  "CD8_identity_module",
  "cDC1_identity_module",
  "cDC2_identity_module",
  "cDC3_identity_module",
  "Innate_NK_like_T_module",
  "MAIT_like_T_module",
  "Memory_B_cell_module",
  "Naive_B_cell_module",
  "Naive_memory_module",
  "pDC_identity_module",
  "Tissue_resident_macrophage_module",
  "Treg_module"
)

# 3. All other modules are classified as functional state
all_names <- original[["Module name"]]
state_names <- setdiff(all_names, identity_names)

# 4. Validate the original library and classification
stopifnot(
  nrow(original) == 62,
  ncol(original) == 4,
  identical(
    names(original),
    c("Module name", "Genes included",
      "Description", "References")
  ),
  !anyDuplicated(all_names),
  all(identity_names %in% all_names),
  length(identity_names) + length(state_names) == 62
)

# 5. Separate rows, keeping original column order and content
identity_table <- original[
  original[["Module name"]] %in% identity_names,
  ,
  drop = FALSE
]

state_table <- original[
  original[["Module name"]] %in% state_names,
  ,
  drop = FALSE
]

# 6. Verify every original row is preserved exactly once
combined <- rbind(identity_table, state_table)

stopifnot(
  nrow(combined) == nrow(original),
  setequal(combined[["Module name"]], all_names),
  identical(
    combined[match(all_names, combined[["Module name"]]), ],
    original
  )
)

# 7. Save two separate Excel files
write_xlsx(
  identity_table,
  "data/identity_modules.xlsx"
)

write_xlsx(
  state_table,
  "data/functional_state_modules.xlsx"
)

# 8. Save a classification report for review
classification <- data.frame(
  Module = all_names,
  Category = ifelse(
    all_names %in% identity_names,
    "identity",
    "functional_state"
  )
)

dir.create("results", showWarnings = FALSE)

write.csv(
  classification,
  "results/original_module_classification.csv",
  row.names = FALSE
)

# 9. Print summary
cat("\nIdentity modules:", nrow(identity_table), "\n")
cat("Functional-state modules:", nrow(state_table), "\n")
cat("Total modules:", nrow(combined), "\n")
cat("Original schema preserved.\n")
cat("Original workbook unchanged.\n")
