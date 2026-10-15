# AI-Driven Gene Module Scoring for Single-Cell Immune Profiling

## 1. Project Overview

This repository contains the R scripts and analysis workflow for a Master's thesis in Bioinformatics at the University of Skövde.

**Thesis title:** Development and Evaluation of an Explainable Gene Module-Based Framework for Single-Cell Immune Profiling

**Students:** Md Reaz Morshed and Mahnoor Tariq  
**Supervisor:** Dr. Salim Ghannoum

The project investigates whether a curated library of **62 immune gene modules** can provide interpretable information about immune cell identity and functional state in single-cell RNA sequencing (scRNA-seq) datasets.

The project includes:

- Evaluation of gene module scoring methods.
- Development of an explainable machine learning framework.
- Comparison with established cell annotation methods.
- Validation across independent immune datasets.

## 2. Main Research Objectives

### Objective 1 — Gene Module Scoring

**Primary responsibility: Reaz**

Evaluate how well different scoring methods measure the activity of the 62 immune gene modules.

Methods specified in the project plan:

1. Seurat AddModuleScore
2. UCell
3. AUCell
4. singscore
5. GSVA

Evaluation will consider AUROC, biological specificity, and stability under changes in cell-type composition.

### Objective 2 — Explainable Machine Learning

**Primary responsibility: Mahnoor**

Use module scores as input features for machine learning models that predict immune cell identity and functional state.

Planned methods include random forest, XGBoost, regularized logistic regression, and SHAP.

### Objective 3 — Cross-Dataset Validation

**Shared responsibility**

Evaluate the framework on independent datasets and compare it with established annotation approaches.

## 3. Datasets and Files

### Main Development Dataset: WuImmune

WuImmune is an immune-cell subset of the Wu et al. (2021) breast cancer atlas.

- 50,063 immune cells
- 29,733 genes
- 26 patients
- 4 major annotated immune populations
- 12 minor annotated immune populations

**Original dataset:**

`data/WuImmune.RData`

**Reclustered working dataset:**

`data/WuImmune_reclustered.rds`

The original dataset was provided by the supervisor. It already contained immune cells, cell-type annotations, and PCA/UMAP reductions.

Reaz performed new PCA, clustering, and UMAP on the immune-cell subset and saved the updated object separately.

### Initial Exploratory Dataset

`data/ICI_cRA_joint_small_for_module_pipeline.RData`

This dataset was used for early exploration of Seurat objects and the gene module library.

**Important:** Gene availability results from this initial dataset must not be assumed to apply to WuImmune.

### Original Gene Module Library

`data/ICI_62_modules_annotation_references.xlsx`

This is the original curated library of 62 immune gene modules, including module names, genes, descriptions, and references.

**This is the main module library for the thesis scoring benchmark.**

## 4. Project Structure

```text
Thesis-project/
│
├── data/
│   ├── WuImmune.RData
│   ├── WuImmune_reclustered.rds
│   ├── ICI_cRA_joint_small_for_module_pipeline.RData
│   ├── ICI_62_modules_annotation_references.xlsx
│   ├── module_list_clean.rds
│   └── individual_overlap_modules.rds
│
├── scripts/
│   ├── 01_explore_obj_small.R
│   ├── 02_explore_62_modules.R
│   ├── 03_check_missing_module_genes.R
│   ├── 04_check_module_overlap.R
│   ├── 05_create_individual_overlap_modules.R
│   ├── 06_explore_WuImmune.R
│   └── 07_recluster_WuImmune.R
│
├── figures/
│   ├── WuImmune_elbow_plot.png
│   ├── WuImmune_umap_major.png
│   ├── WuImmune_umap_minor.png
│   └── WuImmune_umap_clusters.png
│
└── README.md
```

## 5. Completed Work — Initial Dataset and Module Exploration

These scripts document the work performed before starting the main WuImmune scoring benchmark.

### Script 01 — Explore Initial Seurat Object

**File:** `scripts/01_explore_obj_small.R`

**What I did:**

- Loaded the initial ICI/cRA Seurat object.
- Explored its dimensions and structure.
- Inspected gene-expression data and cell metadata.
- Became familiar with Seurat objects.

**Why this work matters:**

It helped establish the data-analysis workflow and understand the structure of single-cell datasets.

**When to use it again:**

When revisiting the initial dataset or documenting preliminary data exploration.

**Status:** COMPLETED — Keep for documentation.

---

### Script 02 — Explore the Original 62 Gene Modules

**File:** `scripts/02_explore_62_modules.R`

**What I did:**

- Loaded the original module library from Excel.
- Examined the 62 module names.
- Explored the genes included in each module.
- Reviewed the structure of the curated library.

**Why this work matters:**

This establishes which biological programs are represented by the module library.

**When to use it again:**

- When describing the original module library in the thesis.
- When reviewing module sizes and gene composition.
- When preparing tables of module names, descriptions, and references.
- When selecting modules for detailed evaluation.

**Status:** COMPLETED — Relevant throughout the project.

---

### Script 03 — Check Missing Module Genes in the Initial Dataset

**File:** `scripts/03_check_missing_module_genes.R`

**What I did:**

- Compared module gene names against genes available in the initial ICI/cRA Seurat object.
- Identified three missing genes:
  - ITGAE
  - CASP1
  - TRAV1-2
- Prepared a dataset-compatible module list.

**Output:**

`data/module_list_clean.rds`

**Important limitation:**

This cleaned module list was prepared for the initial ICI/cRA dataset.

It should **not automatically be used for WuImmune**, because gene availability may differ between datasets.

**When to use it again:**

- When analyzing the initial ICI/cRA dataset.
- When comparing missing-gene patterns across datasets.
- When documenting how module gene coverage changes between datasets.
- As a reference for developing the WuImmune missing-gene checking script.

**Status:** COMPLETED FOR INITIAL DATASET — Repeat the check separately for WuImmune.

---

### Script 04 — Investigate Overlapping Genes Between Modules

**File:** `scripts/04_check_module_overlap.R`

**What I did:**

- Compared the genes in all 62 original modules.
- Identified genes shared between different modules.
- Investigated redundancy within the module library.

**Results:**

- 116 module pairs share at least one gene.
- 102 unique genes occur in multiple modules.

**Why this work matters:**

Overlapping genes can make module scores correlated because the same genes contribute to multiple biological programs.

However, overlap is not automatically an error. A gene can legitimately participate in multiple immune functions.

**When to use it again:**

1. **During module library quality assessment:** Report the extent of redundancy in the original library.
2. **During scoring-method evaluation:** Investigate whether strongly correlated module scores are partly explained by shared genes.
3. **During machine learning interpretation:** Consider overlap when interpreting feature importance or SHAP values.
4. **During thesis writing:** Describe the limitations of overlapping gene modules.

**Important:** The project plan explicitly requires assessment of redundancy between modules.

**Status:** COMPLETED — Important supporting analysis for the main thesis.

---

### Script 05 — Experimental Non-Overlapping Module Analysis

**File:** `scripts/05_create_individual_overlap_modules.R`

**What I did:**

Investigated an alternative representation of the original module library.

The experiment involved:

1. Identifying the 102 genes shared between modules.
2. Removing shared genes from the original modules.
3. Keeping the remaining module-specific genes.
4. Creating separate one-gene modules for the shared genes.

**Results:**

- 62 original modules were examined.
- 102 overlapping genes were identified.
- 12 original modules became empty after shared-gene removal.
- 50 original modules remained nonempty.
- 102 separate one-gene modules were created.

**Saved output:**

`data/individual_overlap_modules.rds`

**Important distinction:**

This was an experimental modification of the module library, not the original 62-module scoring framework.

A potential alternative library could contain 50 nonempty modules plus 102 one-gene modules, giving 152 features. However, this combined library has not yet been established as the main analysis library.

**When could I use this work again?**

**Possible future use A — Sensitivity analysis**

After evaluating the original 62 modules, compare the results with a non-overlapping representation.

Research question:

"Do shared genes significantly influence module scoring results?"

**Possible future use B — Redundancy investigation**

Examine whether removing shared genes changes correlations between module scores.

Research question:

"Are highly correlated module scores driven by overlapping genes?"

**Possible future use C — Machine learning comparison**

If scientifically justified, compare model performance using original modules and a modified representation.

Research question:

"Does removing shared genes improve model interpretability or prediction performance?"

**Possible future use D — Limitations and future work**

Discuss the effect of module overlap even if the experimental library is not used for final benchmarking.

**Important limitations:**

- Single-gene features are not equivalent to multi-gene module scores.
- Removing overlapping genes can remove biologically important information.
- Some original modules become empty.
- Any comparison requires careful interpretation.
- The experimental representation is not part of the required main scoring benchmark in the current project plan.

**Decision:**

Keep this script and its outputs, but **do not use the modified modules for the main benchmark at this stage**.

First evaluate the original 62 modules. Revisit this experiment later only if needed and agreed with the supervisor.

**Status:** EXPLORATORY WORK COMPLETED — OPTIONAL FUTURE ANALYSIS.

## 6. Completed Work — WuImmune Dataset

### Script 06 — Explore WuImmune

**File:** `scripts/06_explore_WuImmune.R`

**What I did:**

- Loaded the WuImmune dataset provided by the supervisor.
- Checked the number of cells and genes.
- Examined major and minor immune cell annotations.
- Identified existing PCA and UMAP reductions.

**Results:**

- 50,063 cells
- 29,733 genes
- 4 major annotated cell populations
- 12 minor annotated cell populations

**Status:** COMPLETED.

---

### Script 07 — Recluster WuImmune

**File:** `scripts/07_recluster_WuImmune.R`

**What I did:**

- Normalized the RNA data.
- Selected 3,000 highly variable genes.
- Scaled the selected features.
- Recalculated PCA.
- Constructed a nearest-neighbor graph using 30 PCs.
- Performed clustering at resolution 0.5.
- Recalculated UMAP using 30 PCs.
- Preserved the original cell annotations.

**Results:**

- 27 computational clusters.
- Updated PCA and UMAP.
- All 50,063 cells retained.
- Original major and minor annotations preserved.

**Saved output:**

`data/WuImmune_reclustered.rds`

**Figures:**

- `figures/WuImmune_elbow_plot.png`
- `figures/WuImmune_umap_major.png`
- `figures/WuImmune_umap_minor.png`
- `figures/WuImmune_umap_clusters.png`

**Important distinction:**

The 27 clusters are computational groups, not 27 confirmed biological cell types.

The original cell-type annotations will be used as reference labels for evaluating module scores.

**Status:** COMPLETED — Main working Seurat object ready.

## 7. Which Module Library Should Be Used?

| Module representation | Purpose | When to use |
|---|---|---|
| Original 62 modules from Excel | Main thesis scoring benchmark | NOW |
| Dataset-compatible original modules | Same original definitions, with unavailable genes excluded for a specific dataset | During scoring |
| Previously cleaned ICI/cRA module list | Initial dataset analysis | Only for the matching dataset |
| Original overlap analysis | Redundancy assessment | During interpretation and thesis writing |
| Modified non-overlapping modules | Optional sensitivity experiment | Later, if justified |
| Individual shared-gene modules | Experimental investigation of shared genes | Later, if justified |

### Main decision

**Use the original 62 gene modules for the main scoring analysis.**

Do not replace them with the experimental non-overlapping representation without a clear scientific reason.

## 8. Current Work — What To Do Next

### Step 08 — Check Original Module Genes Against WuImmune

**Planned script:** `scripts/08_check_WuImmune_module_genes.R`

**Inputs:**

- `data/WuImmune_reclustered.rds`
- `data/ICI_62_modules_annotation_references.xlsx`

**Tasks:**

1. Load the reclustered WuImmune object.
2. Load the original 62 modules.
3. Compare module genes with WuImmune gene names.
4. Count present and missing genes for every module.
5. Save a report of missing genes.
6. Prepare WuImmune-compatible module lists while preserving the original library separately.
7. Identify modules with insufficient available genes for reliable scoring.

**Expected outputs:**

- A missing-gene summary for WuImmune.
- A WuImmune-compatible version of the original module library.
- Documentation of module gene coverage.

**Status:** NEXT TASK — NOT YET COMPLETED.

---

### Step 09 — Review Module Redundancy and Marker Coverage

**Previous work to reuse:**

- `04_check_module_overlap.R`
- Original 62-module library

**Tasks:**

- Use the completed overlap results to document redundancy.
- Check which important markers of annotated immune populations are absent from the module library.
- Record limitations of the module library.

**Status:** Overlap analysis completed; marker coverage assessment still pending.

---

### Step 10 — Calculate Module Scores

**Planned scoring methods:**

- AddModuleScore
- UCell
- AUCell
- singscore
- GSVA

**Tasks:**

- Apply all 62 original modules to WuImmune.
- Calculate one score per cell per module for each method.
- Save results separately for each scoring method.
- Document missing genes and method-specific requirements.

For 50,063 cells and 62 modules, each method can produce up to:

**50,063 × 62 = 3,103,906 cell-module scores.**

**Status:** PLANNED.

---

### Step 11 — Evaluate Module Scoring Methods

**Tasks:**

- Compare identity-module scores against annotated immune cell types.
- Compare functional-state module scores against annotated functional states.
- Evaluate AUROC and biological specificity.
- Examine score distributions across reference populations.
- Assess stability under altered cell-type proportions.
- Investigate possible circularity where reference state labels and module definitions share marker genes.

**Prioritized modules:**

- Eight cell-identity modules.
- Two functional-state modules: proliferation and interferon response.

**Status:** PLANNED.

---

### Step 12 — Optional Module Overlap Sensitivity Analysis

**Previous work to reuse:**

- `04_check_module_overlap.R`
- `05_create_individual_overlap_modules.R`
- `individual_overlap_modules.rds`

**Only consider this after the original 62-module benchmark.**

Possible tasks:

- Compare original and modified module representations.
- Examine changes in module-score correlations.
- Determine whether removing shared genes affects biological specificity.
- Evaluate whether the modified representation adds useful information.

**Status:** OPTIONAL — Not required before the main benchmark.

---

### Step 13 — Build the Module-Score Reference Atlas

Use the selected scoring method to generate module-score profiles for annotated WuImmune cells.

These scores will support the downstream machine learning framework.

**Status:** PLANNED.

---

### Step 14 — Machine Learning and Explainability

Planned methods include:

- Random forest
- XGBoost
- Regularized logistic regression
- SHAP

Separate prediction tasks will address immune cell identity and functional state.

**Status:** PLANNED — Primarily Mahnoor's responsibility.

---

### Step 15 — Cross-Dataset Validation

Evaluate the framework using independent datasets, including PBMCpedia, RA, and ICI-induced arthritis data.

Compare the framework against established annotation methods.

**Status:** PLANNED — Shared responsibility.

## 9. Work Progress Summary

| Task | Status |
|---|---|
| Explore initial Seurat object | COMPLETED |
| Explore original 62 modules | COMPLETED |
| Check missing genes in initial dataset | COMPLETED |
| Check overlap among 62 modules | COMPLETED |
| Experiment with non-overlapping modules | EXPLORATORY COMPLETED |
| Explore WuImmune | COMPLETED |
| Recluster WuImmune | COMPLETED |
| Check 62 modules against WuImmune | NEXT |
| Assess uncovered immune markers | PENDING |
| Run module scoring methods | PENDING |
| Evaluate scoring methods | PENDING |
| Optional overlap sensitivity analysis | OPTIONAL |
| Build module-score reference atlas | PENDING |
| Develop explainable ML framework | PENDING |
| Cross-dataset validation | PENDING |

## 10. Important Research Decisions

1. Preserve the original 62 gene module definitions.
2. Use `WuImmune_reclustered.rds` as the main WuImmune working object.
3. Do not reuse gene lists cleaned for another dataset without checking gene availability.
4. Keep all exploratory scripts for reproducibility.
5. Treat gene overlap as a possible source of redundancy, not automatically as an error.
6. Evaluate original modules before considering experimental modifications.
7. Distinguish computational clusters from biological cell-type annotations.
8. Use reference annotations to evaluate module scores.
9. Document scoring methods and results separately.
10. Keep confidential datasets and derived data within the authorized project environment.

## 11. Software

- R
- Seurat v5
- Visual Studio Code
- macOS

Additional package versions and dependencies will be documented as the analysis progresses.

## 12. Data Confidentiality

Some project datasets contain confidential research data.

Raw datasets, derived Seurat objects, and restricted annotations must not be distributed outside the authorized research team.

Data and outputs must be stored and shared only according to the supervisor's and University's approved data-management requirements.

## 13. Reference

Wu, S. Z., et al. (2021). A single-cell and spatially resolved atlas of human breast cancers. *Nature Genetics*, 53, 1334–1347.

https://doi.org/10.1038/s41588-021-00911-1

Further references are provided in the thesis project plan.

---

## Current Project Status

**Completed:** Initial module exploration, missing-gene checking for the initial dataset, module overlap analysis, experimental non-overlapping modules, WuImmune exploration, and WuImmune reclustering.

**Next:** Check the original 62 gene modules against the reclustered WuImmune dataset.

**Main analysis:** Benchmark the original 62 modules using the scoring methods specified in the project plan.

**Later optional work:** Revisit the non-overlapping module experiment to investigate redundancy or scoring sensitivity if scientifically justified.