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

**---**

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

**---**

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

**---**

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

**---**

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

**---**

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
\| Module representation | Purpose | When to use |

\|---|---|---|

\| Original 62 modules from Excel | Main thesis scoring benchmark | NOW |

\| Dataset-compatible original modules | Same original definitions, with unavailable genes excluded for a specific dataset | During scoring |

\| Previously cleaned ICI/cRA module list | Initial dataset analysis | Only for the matching dataset |

\| Original overlap analysis | Redundancy assessment | During interpretation and thesis writing |

\| Modified non-overlapping modules | Optional sensitivity experiment | Later, if justified |

\| Individual shared-gene modules | Experimental investigation of shared genes | Later, if justified |

### Main decision
**Use the original 62 gene modules for the main scoring analysis.**

Do not replace them with the experimental non-overlapping representation without a clear scientific reason.

## 8. WuImmune module preparation and scoring design

### Original library and dataset-specific availability

The authoritative source is `data/ICI_62_modules_annotation_references.xlsx`. The **actual gene-set input for the WuImmune benchmark** is `data/WuImmune_62_modules_available.rds`, a named list containing all **62 original modules**, with only genes absent from WuImmune excluded. Three unavailable gene entries were removed across two modules; module names and biologically meaningful overlap between modules were retained. The resulting module sizes range from **3 to 28 genes** (median **6.5**). The earlier `data/module_list_clean.rds` belongs to the initial ICI/cRA dataset and **must not be substituted** for the WuImmune scoring library.

The common working object is `data/WuImmune_reclustered.rds` (**50,063 cells × 29,733 genes**). Where applicable, scoring uses the **Seurat RNA `data` layer (normalized, log-transformed expression)**. Methods differ internally in ranking, filtering, weighting, normalization, and input expectations; therefore, the same input library does not imply identical mathematical assumptions.

The intended full score matrix per successful method has **50,063 × 62 = 3,103,906 scores**. All 62 modules must be retained, and cell IDs and module names must match the source. Missing, infinite, and constant scores are assessed before accepting an output. Method-specific limitations and exceptions must be documented rather than hidden.

### Scoring-first workflow

1. Check installation and interface; run a small compatibility test when needed.
2. Score all cells and all 62 original WuImmune-compatible modules without altering module definitions.
3. Validate score dimensions, cell/module identity, missing/non-finite values, and module variability.
4. Save method-specific scores and quality reports only after successful validation.
5. Record measured scoring time, process time, and macOS memory statistics. Installation time remains `NA` unless actually timed.
6. Evaluate biological specificity, AUROC, stability, and interpretability **after** completing the scoring phase.

**Benchmark terminology:** `scoring_seconds` measures the scoring operation inside R; `total_execution_seconds` measures the script's elapsed time; external `/usr/bin/time -l` elapsed time includes process startup/shutdown and may differ. On macOS, `maximum resident set size` and `peak memory footprint` are different metrics. Convert bytes to GiB using `bytes / 1024^3`. Never substitute a tool's estimated intermediate allocation for observed peak memory.

## 9. Scoring methods: actual progress as of 21 October 2026

| Method | Version / implementation | Status | Valid full scores | Measured scoring time | Practical note |
|---|---|---|---|---|---|
| Seurat AddModuleScore | Seurat built-in | **Completed** | 62 × 50,063 | Not separately timed | Simple integration with existing Seurat object; benchmark timing pending |
| UCell | 2.16.0 | **Completed** | 62 × 50,063 | **196.95 s** | Single-cell rank-based scoring |
| AUCell | 1.34.0 | **Completed** | 62 × 50,063 | **97.27 s** | 101 batches of 500 cells; 5% ranking threshold |
| JASMINE | Authors' R source, 11 Oct 2021 | **Completed, implementation validation pending** | 62 × 50,063 | **3,195.19 s (~53 min 15 s)** | Manual source preparation; optimized computation needs reference cross-check |
| VAM | 1.1.0 | **Paused / skipped** | No | — | `vamForSeurat` produced all-NA results; direct `vam` worked on one test module |
| SCSE (Single-Cell Signature Explorer) | Standalone implementation | **Skipped** | No | — | Additional standalone Go workflow and input/UMI requirements did not fit the current common pipeline |
| VISION | 3.0.2 | **Completed** | 62 × 50,063 | **10.09 s** | Lowered signature detection threshold to retain `Th2_module` |
| GSVA | 2.6.12 | **Full run in progress** | Not yet confirmed | Not yet measured | Gaussian sparse GSVA on all cells; prolonged ECDF computation |
| ANS | Not established | **Deferred** | No | — | Implementation and requirements to review later |
| singscore | Not established | **Not yet run** | No | — | In original project plan; no verified full benchmark |

**Interpretation:** JASMINE's approximately 58-minute recollection should be corrected using the measured log: **3,195.19 seconds for scoring (53.25 minutes)** and **3,266.37 seconds for the complete R process (54.44 minutes)**. Neither measurement is exactly 58 minutes. The reason VAM was paused was **invalid scoring output**, not merely warning messages; its installation/build warnings and the observed NA failure should be documented separately. SCSE was skipped because of its additional implementation/input requirements.

### 9.1 AddModuleScore — Script 11

- Script: `scripts/11_score_WuImmune_AddModuleScore.R`
- Uses Seurat's built-in `AddModuleScore()`; no separate scoring package.
- Successful output: `results/WuImmune_AddModuleScore_scores.csv`, `results/WuImmune_AddModuleScore_scores.rds`, `data/WuImmune_AddModuleScore_scored.rds`.
- Full scoring completed; **a directly comparable timed runtime and peak-memory benchmark have not been established**. Do not invent these measurements.

### 9.2 UCell — Script 12

- Script: `scripts/12_score_WuImmune_UCell.R`
- Install: `BiocManager::install("UCell", ask = FALSE, update = FALSE)`.
- Version: **2.16.0**. Scored all **62 × 50,063** combinations with **zero NA**.
- Measured scoring: **196.95 s**; complete R process: **269.18 s**.
- macOS maximum resident set size: **1,990,639,616 bytes** (~1.85 GiB); peak footprint: **6,882,411,968 bytes** (~6.41 GiB).
- Dedicated single-cell rank-based method; note the method's ranking and implementation settings in the scoring script.

### 9.3 AUCell — Script 13

- Script: `scripts/13_score_WuImmune_AUCell.R`
- Install: `BiocManager::install("AUCell", ask = FALSE, update = FALSE)`.
- Version: **1.34.0**. All **62 × 50,063** scores produced, **zero NA**.
- **101 batches** of up to **500 cells**, `auc_max_rank = 1487` (~5% of 29,733 genes).
- Ranking: **29.19 s**; AUC calculation: **60.69 s**; measured scoring operation: **97.27 s**; total execution: **170.73 s**.
- Maximum resident set size: **2,051,719,168 bytes** (~1.911 GiB); peak footprint: **5,427,327,168 bytes** (~5.055 GiB).
- Benchmark: `results/benchmark_AUCell.csv`. Installation time was not measured.

### 9.4 JASMINE — Script 14

- Script: `scripts/14_score_WuImmune_JASMINE.R`; compatibility test: `scripts/14a_test_WuImmune_JASMINE.R`.
- Authors' source: `https://github.com/NNoureen/JASMINE`.
- Source file stored at `external/JASMINE/JASMINE_V1_11October2021.r`.
- Download example:

  ```bash
  mkdir -p external/JASMINE
  curl -L https://raw.githubusercontent.com/NNoureen/JASMINE/main/JASMINE_V1_11October2021.r -o external/JASMINE/JASMINE_V1_11October2021.r
  ```

- Source-loading issue: an unguarded example call on line 25 (`Result = JASMINE(data,genes,method =c('oddsratio','likelihood'))`) was commented out so sourcing did not execute an example. This does **not** itself change the scoring formula.
- Required manual source loading and dependencies (including `stringr` and `GSA`). Compared with package-manager installations, setup was more demanding.
- An optimized rank/odds-ratio batch workflow with a safeguard for constant-value normalization completed **62 × 50,063** scores with **zero NA**. **Numerical agreement with the original authors' function on a controlled subset still requires verification** before using the optimized scores as fully validated reference-method results.
- Measured scoring: **3,195.19 s = 53 min 15.19 s**. Total process: **3,266.37 s = 54 min 26.37 s**.
- Maximum resident set size: **2,001,649,664 bytes** (~1.86 GiB); peak footprint: **5,025,691,072 bytes** (~4.68 GiB).
- Avoid independently normalizing batches if the intended method's normalization depends on the whole population of cells.

### 9.5 VAM — Script 15 (paused)

- Version **1.1.0**; attempted using `vamForSeurat` on the WuImmune object.
- The Seurat wrapper produced **3,103,906 NA scores** (all 62 modules × 50,063 cells), including when `gamma = FALSE`.
- A direct `VAM::vam` call for one module produced finite scores, suggesting a wrapper or workflow issue, but the cause has **not been confirmed**.
- Diagnostic output: `results/WuImmune_VAM_diagnostics.csv` and the run log.
- **Decision:** skip/hold VAM for the current scoring comparison; do not report it as a successful benchmark and do not interpret its all-NA matrix as biological results. Installation or runtime warnings are distinct from the confirmed score-validation failure.

### 9.6 SCSE — skipped because of requirements

Single-Cell Signature Explorer (SCSE) was reviewed but **not executed** in the common R/Seurat scoring benchmark. Its standalone Go implementation and different UMI/input workflow would require additional preparation and reduce direct workflow comparability. This is a **practical exclusion**, not evidence that SCSE performs poorly. No runtime, memory, or biological result was measured.

### 9.7 VISION — Script 16

- Script: `scripts/16_score_WuImmune_VISION.R`.
- Installed version **3.0.2** using:

  ```r
  remotes::install_github("YosefLab/VISION", dependencies = NA, upgrade = "never")
  ```

- GitHub/source compilation produced warnings during installation, but the package installed and the 100-cell/one-module test succeeded with finite scores.
- Built **62 weighted signatures** with unit gene weights; `min_signature_genes = 3`, `sig_norm_method = NULL`, `sig_gene_importance = TRUE`.
- Initially `sig_gene_threshold = 0.01` retained only **61** modules, excluding `Th2_module` (`IL4`, `IL5`, `IL13`, `GATA3`). The full successful run used **`sig_gene_threshold = 0`**, retaining all **62** modules. The initial VISION object reported its own 0.1% gene-detection filter, so this change should **not** be described as disabling every internal filter.
- Full results: **50,063 cells × 62 modules**, **zero missing**, **zero non-finite**, **zero constant-score modules**; range **−0.4310203 to 10.2924**. Negative values are valid.
- Scoring: **10.09 s**; R script elapsed: **19 s**; external wall time: **21.38 s**.
- Maximum resident set size: **2,362,179,584 bytes** (~2.20 GiB); peak footprint: **7,466,237,312 bytes** (~6.95 GiB).
- Outputs: `results/WuImmune_VISION_scores.csv`, `.rds`, `results/WuImmune_VISION_quality_report.csv`, `results/benchmark_VISION.csv`, `results/WuImmune_VISION_run_log.txt`.
- **Reproducibility correction:** both Script 16's benchmark-writing code and the saved benchmark CSV must show **`sig_gene_threshold = 0`**, not the older `0.01`. The benchmark was initially printed with the stale `0.01` value even after the successful full run; verify the correction has actually been saved.
- The common Seurat log-normalized RNA layer was used, while VISION's internal processing has its own conventions. This input choice should be noted in the methods section and, if necessary, revisited in a sensitivity analysis.

### 9.8 GSVA — Script 17 (running, not completed)

- Version **2.6.12**, installed from Bioconductor with:

  ```r
  BiocManager::install("GSVA", ask = FALSE, update = FALSE)
  ```

- `GSVA::gsvaParam()` is available. The test on **100 cells** using `kcdf = "Gaussian"`, `sparse = TRUE`, `minSize = 1` produced **58 modules × 100 cells**, zero missing and zero non-finite scores, in **0.62 s**. Four signatures were filtered in the small test.
- A detection-based diagnostic found at least two variably detected genes for **all 62 modules across the full dataset**. This was encouraging but did not guarantee full GSVA retention.
- Full script: `scripts/17_score_WuImmune_GSVA.R`.
- Full-run settings: normalized RNA `data` layer; **29,733 genes, 50,063 cells, 62 modules**; `kcdf = "Gaussian"`, `sparse = TRUE`, `minSize = 1`, `maxSize = Inf`.
- At the last observed output, GSVA had reached **row-wise Gaussian-kernel ECDF estimation** and reported splitting a large intermediate calculation into two chunks, including a **7.200 GiB** estimate. It reported filtering **4,026 constant rows** and **1,731 constant-nonzero rows**.
- A process check at **1 h 36 min 28 s** showed the actual R process (PID 8000) using approximately **98.4% CPU**, consistent with ongoing computation rather than a confirmed freeze. This PID is historical and will change on reruns.
- **Status: still running at last observation.** No final GSVA scores, elapsed benchmark, memory peak, or completed status should be claimed until the script finishes and all **62 × 50,063** scores pass validation.
- Do not assume the full run will be as fast as the 100-cell test; GSVA's ECDF calculations depend on the full cell collection. Do not split the full analysis into independent scoring batches without assessing the effect on the GSVA estimates.

### 9.9 Deferred methods

- **ANS:** implementation and special requirements to be reviewed later; not installed, scored, or benchmarked in the confirmed record.
- **singscore:** listed in the original project plan; no confirmed full WuImmune run yet.
- **Other proposed methods:** assess scientific relevance, compatibility, installation difficulty, and resource requirements before adding them to the comparison.

## 10. Benchmark summary and interpretation

The values below are **measured results from completed scoring runs**, not estimates. The methods differ algorithmically, and runtime comparisons do not establish biological accuracy.

| Method | Scoring seconds | R total seconds | Max resident (GiB) | Peak footprint (GiB) | Notes |
|---|---:|---:|---:|---:|---|
| AddModuleScore | — | — | — | — | Scores completed; timed benchmark pending |
| UCell | 196.95 | 269.18 | 1.85 | 6.41 | Completed |
| AUCell | 97.27 | 170.73 | 1.911 | 5.055 | Completed |
| JASMINE | 3,195.19 | 3,266.37 | 1.86 | 4.68 | Optimized implementation; original-function validation pending |
| VISION | 10.09 | 19.00 | 2.20 | 6.95 | Completed; `sig_gene_threshold = 0` |
| VAM | — | — | — | — | Invalid all-NA wrapper output; skipped |
| SCSE | — | — | — | — | Not run; requirements incompatible with current workflow |
| GSVA | — | — | — | — | Full run in progress at last check |

For complete provenance, consult each `results/benchmark_<METHOD>.csv` and its matching execution log. **Installation seconds are not known** for these methods unless separately measured. Installation difficulty is a qualitative assessment of this project's experience, not a measured performance outcome.

Example macOS execution command:

```bash
/usr/bin/time -l Rscript scripts/17_score_WuImmune_GSVA.R 2>&1 | tee results/WuImmune_GSVA_run_log.txt
```

**Important:** `tee` saves the displayed output but can mask a failing pipeline's exit code in some shells unless `pipefail` is enabled. Verify the final summary and expected output files before marking a method completed.

## 11. Updated project structure (selected files)

```text
Thesis-project/
├── data/
│   ├── WuImmune.RData
│   ├── WuImmune_reclustered.rds
│   ├── ICI_62_modules_annotation_references.xlsx
│   ├── WuImmune_62_modules_available.rds
│   ├── module_list_clean.rds                 # Initial dataset only
│   └── individual_overlap_modules.rds       # Experimental; not benchmark input
├── scripts/
│   ├── 01_explore_obj_small.R
│   ├── 02_explore_62_modules.R
│   ├── 03_check_missing_module_genes.R
│   ├── 04_check_module_overlap.R
│   ├── 05_create_individual_overlap_modules.R
│   ├── 06_explore_WuImmune.R
│   ├── 07_recluster_WuImmune.R
│   ├── 11_score_WuImmune_AddModuleScore.R
│   ├── 12_score_WuImmune_UCell.R
│   ├── 13_score_WuImmune_AUCell.R
│   ├── 14_score_WuImmune_JASMINE.R
│   ├── 14a_test_WuImmune_JASMINE.R
│   ├── 15_score_WuImmune_VAM.R              # Attempted; paused
│   ├── 16_score_WuImmune_VISION.R
│   └── 17_score_WuImmune_GSVA.R            # Full run in progress
├── external/
│   └── JASMINE/
│       └── JASMINE_V1_11October2021.r
├── results/
│   ├── WuImmune_AddModuleScore_scores.csv
│   ├── WuImmune_UCell_scores.csv
│   ├── WuImmune_AUCell_scores.csv
│   ├── WuImmune_JASMINE_scores.csv
│   ├── WuImmune_VISION_scores.csv
│   ├── WuImmune_VISION_quality_report.csv
│   ├── WuImmune_VAM_diagnostics.csv
│   ├── benchmark_AUCell.csv
│   ├── benchmark_VISION.csv
│   └── WuImmune_GSVA_run_log.txt            # In progress
├── figures/
│   ├── WuImmune_elbow_plot.png
│   ├── WuImmune_umap_major.png
│   ├── WuImmune_umap_minor.png
│   └── WuImmune_umap_clusters.png
└── README.md
```

The tree is a **selected inventory**, not an assertion that every optional filename exists. In particular, GSVA score and benchmark files are created only if full validation succeeds. Other preliminary scripts, figures, and reports may also be present.

## 12. Work progress and next actions

| Task | Status |
|---|---|
| Explore initial Seurat object and module library | Completed |
| Check missing genes in initial dataset | Completed |
| Assess gene overlap and experimental non-overlapping representation | Completed exploratory work |
| Explore and recluster WuImmune | Completed |
| Prepare original 62 modules for WuImmune | Completed |
| Score AddModuleScore, UCell, AUCell | Completed |
| Score JASMINE | Completed; optimized-vs-original validation pending |
| Score VISION | Completed; verify benchmark threshold correction |
| Diagnose VAM | Attempted; paused because of invalid scores |
| Assess SCSE | Skipped because of implementation/input requirements |
| Score GSVA | Full run in progress at last observation |
| Evaluate singscore and ANS | Pending / deferred |
| Review marker coverage and module redundancy | Overlap assessed; marker coverage review pending |
| Evaluate biological specificity, AUROC, stability | Pending until scoring comparisons are ready |
| Build module-score reference atlas | Planned |
| Explainable ML (random forest, XGBoost, regularized logistic regression, SHAP) | Planned; primarily Mahnoor |
| Cross-dataset validation (PBMCpedia, RA, ICI-induced arthritis) | Planned; shared |

**Immediate next action:** wait for the current GSVA process to finish or report an error; review the final log and memory data. If GSVA cannot complete with available resources, document that as a computational limitation rather than recording an invented benchmark. Separately verify JASMINE's optimized scores against the original function and confirm VISION's saved threshold field is `0`.

## 13. Software, reproducibility, and limitations

- Environment: **macOS on Apple Silicon**, **R 4.6.1**, **Seurat v5**, Visual Studio Code, Bioconductor packages as documented above.
- Run all commands from the repository root (`Thesis-project`), using **relative paths**.
- Preserve the original 62 module definitions and record dataset-specific gene availability.
- Preserve the original cell annotations separately from computational clustering. WuImmune reclustering produced **27 clusters**, not 27 experimentally confirmed cell types.
- Do not treat a method as biologically validated merely because its scores are finite. Biological specificity, AUROC, and stability remain separate analyses.
- Document potential circularity when gene markers used in a module overlap with markers used to define reference labels.
- Treat overlapping module genes as a possible source of correlated scores and reduced feature independence, not automatically an error.
- Preserve method-specific score scales; do not assume all methods produce scores on [0, 1].
- Record installation and usability limitations separately from runtime and memory results.
- Keep sensitive datasets, raw expression objects, and derived results within the authorized research environment; do not publish restricted data or patient information.

## 14. Reference

Wu, S. Z., et al. (2021). A single-cell and spatially resolved atlas of human breast cancers. *Nature Genetics*, 53, 1334–1347. https://doi.org/10.1038/s41588-021-00911-1

Additional biological and module-specific references are recorded in `data/ICI_62_modules_annotation_references.xlsx` and the thesis project plan.

---

**README status:** Updated **21 October 2026** using the completed method logs and the most recent GSVA progress reported. Update the GSVA section and benchmark table after a validated final run; do not change its status merely because the process is consuming CPU.

## Gene Module Scoring and Method Comparison

Six gene-module scoring methods were evaluated on the WuImmune single-cell RNA-seq dataset (**50,063 cells and 62 immune gene modules**):

- **AddModuleScore**
- **UCell**
- **AUCell**
- **JASMINE** (custom optimized implementation)
- **VISION**
- **GSVA**

### Evaluation Approach

The methods were compared using three criteria:

1. **Computational performance:** Scoring runtime, peak memory usage, and practical installation considerations were evaluated.
2. **Biological performance:** AUROC was calculated for 10 selected immune identity modules against existing WuImmune cell-type annotations to assess how well each method distinguished relevant immune-cell populations.
3. **Method agreement:** Spearman rank correlations were calculated across all 62 modules and 50,063 cells for all 15 method pairs (930 correlations), assessing similarity in cell rankings without requiring cell-type labels.

### Key Findings

| Method | Mean AUROC (10 modules) | Scoring time |
|---|---:|---:|
| VISION | **0.8926** | **10.09 s** |
| AddModuleScore | 0.8918 | 18.30 s |
| UCell | 0.8831 | 196.95 s |
| AUCell | 0.8824 | 97.27 s |
| GSVA | 0.8764 | 9,965.57 s |
| JASMINE | 0.8059 | 3,195.19 s |

The highest method agreement was observed between **AddModuleScore and VISION** (mean Spearman ρ = 0.9628), followed by **UCell and AUCell** (ρ = 0.9369).

### Preliminary Method Selection

**VISION and AddModuleScore** emerged as the strongest candidates based on biological discrimination and computational efficiency. **UCell** was retained as a complementary rank-based method for further evaluation.

These findings support method shortlisting rather than establishing a universally best scoring algorithm. AUROC measures agreement with existing annotations, and high inter-method correlation does not necessarily indicate biological accuracy. The custom JASMINE implementation requires further validation.

The selected methods and their saved module scores can be forwarded for downstream analysis and validation.