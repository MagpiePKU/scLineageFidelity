# scLineageFidelity

An R package for lineage fidelity analysis and Maximum Likelihood Estimation (MLE) cell type annotation of single-cell RNA sequencing data.

scLineageFidelity is developed by Yi Zhang's lab in Institute of Genetics and Molecular Medicine, Chinese Institutes for Medical Research, Beijing, China. 

## Table of Contents

- [Introduction](#introduction)
- [Installation](#installation)
- [Dependencies](#dependencies)
- [Main Features](#main-features)
- [Notes](#notes)
- [Citation](#citation)
- [License](#license)
- [Contact](#contact)
- [Usage Examples](#usage-examples)
- [Data Format Requirements](#data-format-requirements)
- [Reference Data](#reference-data)
- [Function Documentation](#function-documentation)
  - [sc_lineage_fidelity.R](#sc_lineage_fidelityr)
  - [Maximal_Likelihood_Estimation_Annotator.r](#maximal_likelihood_estimation_annotatorr)

## Introduction

This package provides two main analytical approaches:

1. **Correlation-based method** (`sc_lineage_fidelity`): Evaluates lineage fidelity by calculating correlation coefficients between single cells and reference datasets
2. **Maximum Likelihood Estimation method** (`Predict_MLE`): Performs more precise cell type annotation and lineage fidelity assessment using a negative binomial distribution model

## Installation

After downloading the code files to your local machine, load them in R:

```r
source("sc_lineage_fidelity.R")
source("Maximal_Likelihood_Estimation_Annotator.r")
```

## Dependencies

### Required Packages
- `WGCNA` - For correlation calculations
- `dplyr` - Data manipulation
- `reshape2` - Data reshaping
- `parallel` - Parallel computing
- `Seurat` - Single-cell data analysis
- `tidyr` - Data tidying
- `ggplot2` - Plotting

### Optional Packages (for parallel computing)
- `pbmcapply` - Parallel computing with progress bar

### Optional Packages (for tree analysis)
- `ape` - Phylogenetic analysis
- `stringr` - String processing

Install dependencies:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("WGCNA", "Seurat"))
install.packages(c("dplyr", "reshape2", "parallel", "tidyr", "ggplot2", 
                   "pbmcapply", "ape", "stringr"))
```

## Main Features

1. **Lineage Fidelity Assessment**: Evaluates the extent to which single cells maintain their original lineage characteristics during differentiation or transformation
2. **Cell Type Annotation**: Performs precise cell type annotation based on reference datasets
3. **Transcriptional Span Calculation**: Quantifies transcriptional heterogeneity of cells on a reference cell type tree
4. **Confidence Interval Analysis**: Provides cell type matching results with 95% confidence intervals

## Notes

1. **Memory Usage**: For large datasets, it is recommended to use the `chunk_size` parameter for batch processing
2. **Parallel Computing**: Adjust the `ncores` parameter according to system resources
3. **Reference Data Normalization**: Ensure reference data is correctly normalized, otherwise results may be inaccurate
4. **Gene Name Matching**: Ensure query data and reference data use the same gene naming convention
5. **Iterative Optimization**: The MLE method automatically performs iterative optimization to improve annotation accuracy, but this will increase computation time

## Citation

If you use this tool, please cite the relevant paper:

```
Xiao Y., Jin W., Chen F., Qian K., Ju L., Zhang Y. 
Evolution of cancer metastases via lineage trans-differentiation; 
in revision (2025)
```

## License

MIT License 

## Contact

For questions or suggestions, please submit an Issue.



## Usage Examples

### Example 1: Calculate Lineage Fidelity Using Correlation Method

```r
library(Seurat)

# Load data
seurat_obj <- readRDS("your_seurat_object.rds")

# Get expression data
count_data <- GetAssayData(seurat_obj, slot = 'data')

# Calculate lineage fidelity
result <- sc_lineage_fidelity(
  scdata = count_data,
  species = 'human',
  ncores = 10
)

# Add results to Seurat object
seurat_obj$lineage_fidelity <- result$lineage_fidelity$lineage_fidelity[
  match(rownames(seurat_obj@meta.data), result$lineage_fidelity$Cell)
]
```

### Example 2: Cell Type Annotation Using MLE Method

```r
library(Seurat)
library(celldex)

# Load reference data (using celldex package)
hpca <- celldex::HumanPrimaryCellAtlasData()
ref_mtx_count <- round(exp(assays(hpca)[[1]]) - 1)
ref_mtx <- t(t(ref_mtx_count) / colSums(ref_mtx_count))

# Load query data
seurat_obj <- readRDS("your_seurat_object.rds")

# Preprocessing
preprocess_result <- Preprocess_for_MLE(
  input_seurat_obj = seurat_obj,
  reference_data_matrix = ref_mtx,
  ncores = 3
)

# MLE prediction
MLE_result_list <- Predict_MLE(
  input_seurat_obj = seurat_obj,
  reference_data_matrix = ref_mtx,
  preprocess_result = preprocess_result,
  top_correlation_percentile_for_MLE = 0.95,
  ncores = 3,
  MLE_hit_probability_threshold = 0.99,
  MLE_hit_ci_95_matches_threshold = 5
)

# Calculate transcriptional span
MLE_span <- Compute_Transcriptional_Span_by_MLE_95CI(
  ref_mtx = ref_mtx,
  MLE_result_list = MLE_result_list,
  ncore = 3
)

# Merge results
MLE_span_annotated <- left_join(
  MLE_span,
  MLE_result_list$MLE_results
)

# Add results to Seurat object
seurat_obj$MLE_cell_type <- MLE_result_list$MLE_results$ref_cell_type[
  match(rownames(seurat_obj@meta.data), MLE_result_list$MLE_results$cell_id)
]
seurat_obj$transcriptional_span <- MLE_span$sum_dist_MLE_95CI[
  match(rownames(seurat_obj@meta.data), MLE_span$cell_id)
]
```

### Example 3: Complete Workflow (from Raw Data to Visualization)

```r
library(Seurat)
library(dplyr)
library(ggplot2)

# 1. Load and preprocess data
seurat_obj <- readRDS("your_seurat_object.rds")
ref_mtx <- readRDS("reference_data.rds")  # Normalized reference matrix

# 2. MLE annotation
preprocess_result <- Preprocess_for_MLE(
  input_seurat_obj = seurat_obj,
  reference_data_matrix = ref_mtx,
  ncores = 3
)

MLE_result_list <- Predict_MLE(
  input_seurat_obj = seurat_obj,
  reference_data_matrix = ref_mtx,
  preprocess_result = preprocess_result,
  ncores = 3
)

# 3. Calculate transcriptional span
MLE_span <- Compute_Transcriptional_Span_by_MLE_95CI(
  ref_mtx = ref_mtx,
  MLE_result_list = MLE_result_list,
  ncore = 3
)

# 4. Merge results
results <- left_join(
  MLE_result_list$MLE_results,
  MLE_span,
  by = "cell_id"
)

# 5. Visualization
ggplot(results, aes(x = ref_cell_type, y = sum_dist_MLE_95CI + 1)) +
  geom_boxplot() +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Cell Type", y = "Transcriptional Span (log scale)")
```

## Data Format Requirements

### Input Data Format

1. **Seurat Object** (`input_seurat_obj`):
   - Must contain an RNA assay
   - Count data can be obtained using `GetAssayData(seurat_obj, slot = 'count', assay = 'RNA')`

2. **Reference Data Matrix** (`reference_data_matrix`):
   - Genes as rows, reference cells/cell types as columns
   - **Must be normalized**: Normalized by column sums (`t(t(ref_mtx_count) / colSums(ref_mtx_count))`)
   - If input is a count matrix, normalize it first
   - If input is a log1p matrix, convert back to count matrix first, then normalize

3. **Single-cell Expression Matrix** (`scdata`):
   - Genes as rows, cells as columns
   - Can be normalized expression matrix or raw count matrix

### Output Data Format

- **Correlation method**: Returns a list containing correlation matrix, top matching results, and lineage fidelity
- **MLE method**: Returns a list containing MLE annotation results, preprocessing information, and transcriptional span

## Reference Data

### Default Reference Data Paths

The code has hardcoded default reference data paths (modify according to your actual situation):

- **Human**: `/gpfs/output/ECS_Research/data/scRNA/HPA_scRNA_ref/ref_expr_human.HPA.rds`
- **Mouse**: `/gpfs/output/ECS_Research/data/scRNA/HPA_scRNA_ref/scMCA_mouse_ref_expr.rds`

**We provide the human and mouse references in `scLineageFidelity_ref_data`**

### Using Custom Reference Data

1. **From celldex package**:
```r
library(celldex)
hpca <- HumanPrimaryCellAtlasData()
ref_mtx_count <- round(exp(assays(hpca)[[1]]) - 1)
ref_mtx <- t(t(ref_mtx_count) / colSums(ref_mtx_count))
```

2. **Load from file**:
```r
ref_mtx <- readRDS("your_reference_data.rds")
# Ensure it is normalized
```

3. **Prepare your own reference data**:
   - Reference data should be a normalized expression matrix
   - Column names should contain cell type information
   - Row names should be gene names (matching query data)






## Function Documentation

### sc_lineage_fidelity.R

#### `sc_lineage_fidelity()`

Calculates lineage fidelity based on correlation analysis.

**Function Signature:**
```r
sc_lineage_fidelity(
  scdata, 
  ref_data = NULL,
  numbers_plot = 3,
  ncores = 10,
  top_n_candidate_for_fidelity = 10,
  species = 'human'
)
```

**Parameters:**

- `scdata` (required): Single-cell expression matrix with genes as rows and cells as columns
- `ref_data` (optional): Reference dataset. If NULL, default reference data will be automatically loaded based on the `species` parameter
  - `species = 'human'`: Loads human HPA reference data
  - `species = 'mouse'`: Loads mouse scMCA reference data
- `numbers_plot` (default=3): Number of top correlations to retain for each cell
- `ncores` (default=10): Number of cores for parallel computing
- `top_n_candidate_for_fidelity` (default=10): Number of candidate cell types used for fidelity calculation
- `species` (default='human'): Species type, either 'human' or 'mouse'

**Return Value:**

Returns a list containing the following elements:

- `cors_matrix`: Complete correlation matrix (reference cell types × query cells)
- `top_cors`: Number of top correlations retained per cell
- `scMCA`: Best matching cell type for each cell (based on highest correlation)
- `scMCA_probility`: Top correlation results for each cell (data frame format with Cell, Cell type, and Score columns)
- `lineage_fidelity`: Lineage fidelity results (data frame format with Cell, lineage_fidelity, and lineage_fidelity_var columns)

**Usage Example:**

```r
# Get expression data from Seurat object
count_data <- GetAssayData(seurat_obj, slot = 'data')

# Calculate lineage fidelity
result <- sc_lineage_fidelity(
  scdata = count_data,
  ref_data = NULL,  # Use default reference data
  numbers_plot = 3,
  ncores = 10,
  top_n_candidate_for_fidelity = 10,
  species = 'human'
)

# View results
head(result$lineage_fidelity)
head(result$scMCA_probility)
```

---

### Maximal_Likelihood_Estimation_Annotator.r

#### `Preprocess_for_MLE()`

Performs data preprocessing for MLE analysis, including feature gene selection and correlation calculation.

**Function Signature:**
```r
Preprocess_for_MLE(
  input_seurat_obj,
  reference_data_matrix,
  nfeatures_for_variable_in_obj = 3000,
  nfeatures_for_variable_in_ref = 40,
  chunk_size = 5000,
  ncores = 1
)
```

**Parameters:**

- `input_seurat_obj` (required): Seurat object that must contain an RNA assay
- `reference_data_matrix` (required): Reference expression matrix, must be normalized (normalized by column sums)
- `nfeatures_for_variable_in_obj` (default=3000): Number of highly variable genes to select from query data
- `nfeatures_for_variable_in_ref` (default=40): Number of highly variable genes to select per cell type from reference data
- `chunk_size` (default=5000): Number of cells per batch when cell count exceeds 5000
- `ncores` (default=1): Number of cores for parallel computing

**Return Value:**

Returns a list containing the following elements:

- `common_vars`: List of common highly variable genes
- `cor_mtx`: Correlation matrix (query cells × reference cell types)
- `count_data`: Raw count matrix
- `ref_mtx`: Reference expression matrix

**Usage Example:**

```r
# Prepare reference data (needs to be normalized)
ref_mtx <- t(t(ref_mtx_count) / colSums(ref_mtx_count))

# Preprocessing
preprocess_result <- Preprocess_for_MLE(
  input_seurat_obj = seurat_obj,
  reference_data_matrix = ref_mtx,
  nfeatures_for_variable_in_obj = 3000,
  nfeatures_for_variable_in_ref = 40,
  chunk_size = 5000,
  ncores = 3
)
```

---

#### `Predict_MLE()`

Performs cell type annotation using Maximum Likelihood Estimation.

**Function Signature:**
```r
Predict_MLE(
  input_seurat_obj,
  reference_data_matrix,
  preprocess_result = NULL,
  top_correlation_percentile_for_MLE = 0.95,
  ncores = 3,
  MLE_hit_probability_threshold = 0.99,
  MLE_hit_ci_95_matches_threshold = 5,
  compute_cells_limit = NULL,
  dispersion_min = 1e-5
)
```

**Parameters:**

- `input_seurat_obj` (required): Seurat object
- `reference_data_matrix` (required): Reference expression matrix (normalized)
- `preprocess_result` (optional): Preprocessing result. If NULL, the function will internally call `Preprocess_for_MLE`
- `top_correlation_percentile_for_MLE` (default=0.95): Correlation percentile threshold for MLE analysis (between 0 and 1)
- `ncores` (default=3): Number of cores for parallel computing
- `MLE_hit_probability_threshold` (default=0.99): MLE probability threshold; iterative optimization will be performed if below this value
- `MLE_hit_ci_95_matches_threshold` (default=5): 95% confidence interval match count threshold; iterative optimization will be performed if above this value
- `compute_cells_limit` (optional): Limit the number of cells to compute (for testing purposes)
- `dispersion_min` (default=1e-5): Minimum dispersion parameter for negative binomial distribution

**Return Value:**

Returns a list containing the following elements:

- `MLE_results`: Data frame containing annotation results for each cell:
  - `best_prob`: Probability of best match
  - `cell_id`: Cell ID
  - `ref_cell_type`: Best cell type predicted by MLE
  - `top_cor_celltype`: Cell type with highest correlation
  - `ci_95_matches`: Number of matches within 95% confidence interval
  - `log_likelihood`: Log-likelihood value
  - `best_cell_types`: Best candidate cell types (semicolon-separated)
  - `compare_set_size`: Number of reference cell types used for comparison
  - `iter_num`: Number of iterations
- `common_vars`: List of common highly variable genes
- `cor_mtx`: Correlation matrix
- `count_data`: Raw count matrix
- `ref_mtx`: Reference expression matrix

**Usage Example:**

```r
# Method 1: Using preprocessing result
preprocess_result <- Preprocess_for_MLE(
  input_seurat_obj = seurat_obj,
  reference_data_matrix = ref_mtx
)

MLE_result_list <- Predict_MLE(
  input_seurat_obj = seurat_obj,
  reference_data_matrix = ref_mtx,
  preprocess_result = preprocess_result,
  top_correlation_percentile_for_MLE = 0.95,
  ncores = 3,
  MLE_hit_probability_threshold = 0.99,
  MLE_hit_ci_95_matches_threshold = 5
)

# Method 2: Without providing preprocessing result (function will automatically preprocess)
MLE_result_list <- Predict_MLE(
  input_seurat_obj = seurat_obj,
  reference_data_matrix = ref_mtx,
  ncores = 3
)

# View results
head(MLE_result_list$MLE_results)
```

---

#### `Compute_Transcriptional_Span_by_MLE_95CI()`

Calculates transcriptional span based on MLE 95% confidence intervals, used to quantify transcriptional heterogeneity of cells.

**Function Signature:**
```r
Compute_Transcriptional_Span_by_MLE_95CI(
  ref_mtx,
  MLE_result_list,
  clustering_method = 'average',
  correlation_method_for_ref_mtx = 'spearman',
  ncore = 3
)
```

**Parameters:**

- `ref_mtx` (required): Reference expression matrix (should be the same matrix used for MLE analysis)
- `MLE_result_list` (required): Return result from `Predict_MLE()`
- `clustering_method` (default='average'): Hierarchical clustering method; default is recommended
- `correlation_method_for_ref_mtx` (default='spearman'): Correlation calculation method for reference matrix; default is recommended
- `ncore` (default=3): Number of cores for parallel computing

**Return Value:**

Returns a data frame containing:

- `cell_id`: Cell ID
- `sum_dist_MLE_95CI`: Transcriptional span on the reference cell type tree (sum of distances)

**Usage Example:**

```r
# Calculate transcriptional span
MLE_span <- Compute_Transcriptional_Span_by_MLE_95CI(
  ref_mtx = ref_mtx,
  MLE_result_list = MLE_result_list,
  ncore = 3
)

# Merge results
MLE_span_annotated <- left_join(
  MLE_span, 
  MLE_result_list$MLE_results
)
```

---

#### `Compute_Individual_Cell_Transcriptional_Span_on_Tree()`

Calculates transcriptional span of a single cell on the reference tree (helper function).

**Function Signature:**
```r
Compute_Individual_Cell_Transcriptional_Span_on_Tree(
  tree,
  node_names,
  cophenetic_mat
)
```

**Parameters:**

- `tree`: hclust object, hierarchical clustering tree of reference cell types
- `node_names`: Vector of node names for which to calculate span
- `cophenetic_mat`: Cophenetic distance matrix built from the hclust object

**Return Value:**

Returns a numeric value representing the sum of distances between nodes.

---

#### `MLE_NB()` (Internal Function)

Calculates Maximum Likelihood Estimation using negative binomial distribution model.

**Function Signature:**
```r
MLE_NB(
  obs,
  ref,
  count_obs_size,
  top_cell_type_for_fine_tune = NULL,
  dispersion_min = 1e-5,
  cell_id = NULL,
  top_correlated_ref_cell_names = NULL
)
```

**Parameters:**

- `obs`: Observed gene expression vector
- `ref`: Reference expression matrix
- `count_obs_size`: Total count of observed cell
- `top_cell_type_for_fine_tune`: Number of top cell types for fine-tuning
- `dispersion_min`: Minimum dispersion
- `cell_id`: Cell ID
- `top_correlated_ref_cell_names`: Names of most correlated reference cell types

---

#### `Select_Dynamic_Marker_Genes()` (Internal Function)

Dynamically selects marker genes for MLE iterative optimization.

**Function Signature:**
```r
Select_Dynamic_Marker_Genes(
  obs,
  ref,
  top_candidates,
  fine_tune_marker_select_num = 10
)
```

**Parameters:**

- `obs`: Observed gene expression vector
- `ref`: Reference expression matrix
- `top_candidates`: Candidate cell type names
- `fine_tune_marker_select_num`: Number of marker genes to select

---

#### `Get_Cell_Type_by_MLE()` (Internal Function)

Performs MLE annotation for a single cell (with iterative optimization).

**Function Signature:**
```r
Get_Cell_Type_by_MLE(
  cell_id,
  cell_expression_counts_full,
  ref_expr,
  count_obs,
  top_correlated_ref_cell_names,
  common_vars,
  iter_num_max = 3,
  fine_tune_marker_select_num = 10,
  top_cell_type_for_fine_tune = 5,
  probability_threshold = 0.99,
  ci_95_matches_threshold = 5,
  dispersion_min = 1e-5
)
```
