############################################################
# 1. Load package for reading Excel files
############################################################

library(readxl)

############################################################
# 2. Set working directory to the data folder
############################################################

setwd("/Users/mdreazmorshed/Desktop/Thesis-project/data")

############################################################
# 3. Check that the Excel file is present
############################################################

list.files()

############################################################
# 4. See the sheet names inside the Excel file
############################################################

excel_sheets("ICI_62_modules_annotation_references.xlsx")

############################################################
# 5. Read the first sheet
############################################################

modules <- read_excel(
  "ICI_62_modules_annotation_references.xlsx",
  sheet = 1
)

############################################################
# 6. Inspect the structure of the file
############################################################

dim(modules)
colnames(modules)
head(modules)
str(modules)

############################################################
# 7. Print the first 20 rows
############################################################

# print(
#   head(modules, 20),
#   width = Inf
# )

############################################################
# 8. Print all 62 module names
#    This lets us see the complete module library provided
#    by Salim before doing any scoring.
############################################################

# print(modules$`Module name`)

############################################################
# 11. Print complete information for modules 1–10
#     We will study these modules carefully before assigning
#     biological categories.
############################################################

for (i in seq_len(nrow(modules))) {

  cat("\n")
  cat("====================================================\n")
  cat("MODULE", i, "\n")
  cat("====================================================\n")

  cat("\nModule name:\n")
  cat(modules$`Module name`[i], "\n")

  cat("\nGenes included:\n")
  cat(modules$`Genes included`[i], "\n")

  cat("\nDescription:\n")
  cat(modules$Description[i], "\n")

  cat("\nReferences:\n")
  cat(modules$References[i], "\n")
}

