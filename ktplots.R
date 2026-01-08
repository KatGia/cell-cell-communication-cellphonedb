## Visualization with `ktplots` (R)
# Once the Python analysis is complete, you can visualize the results in R using the `ktplots` library. This allows you to create **Heatmaps**, Dot Plots, and **Chord Diagrams** of cell-cell interactions.

### 1. Requirements & Setup
library(ktplots)
library(Seurat)
library(SingleCellExperiment)

# Load the CellphoneDB output files
# Note: Check the exact filename (timestamp) in your results folder
pvals_path <- "results/analisi_2026_2/statistical_analysis_pvalues_...txt"
means_path <- "results/analisi_2026_2/statistical_analysis_means_...txt"
decon_path <- "results/analisi_2026_2/statistical_analysis_deconvoluted_...txt"

pvals_Cate <- read.delim(pvals_path, check.names = FALSE)
means_Cate <- read.delim(means_path, check.names = FALSE)
decon_Cate <- read.delim(decon_path, check.names = FALSE)


### 2. Global Interaction Heatmap
# The heatmap provides a quick overview of the **number of significant interactions** between different cell types.
# Basic heatmap
plot_cpdb_heatmap(pvals = pvals_Cate, cellheight = 10, cellwidth = 10)

# Custom heatmap
plot_cpdb_heatmap(
    pvals = pvals_Cate, 
    title = "Significant Interactions Count",
    figsize = c(6, 6)
)


### 3. Detailed Dot Plots
# To visualize specific ligand-receptor pairs, you need to convert your Seurat object to a SingleCellExperiment object.
# 1. Convert Seurat to SingleCellExperiment
sce <- as.SingleCellExperiment(pbmc_normalized)

# 2. Match Metadata names with CellphoneDB results
# (Remove spaces and + to ensure names match exactly)
sce$scType_CPDB <- gsub(" ", "_", sce$scType_classification)
sce$scType_CPDB <- gsub("\\+", "", sce$scType_CPDB)

# 3. Clean column names (Python often replaces '-' with '.')
colnames(means_Cate) <- gsub("\\.", "_", colnames(means_Cate))
colnames(pvals_Cate) <- gsub("\\.", "_", colnames(pvals_Cate))

# 4. Plot specific cell-type interactions
ktplots::plot_cpdb(
    scdata = sce,
    cell_type1 = "Naive_B_cells",
    cell_type2 = "Classical_Monocytes",
    celltype_key = "scType_CPDB",
    means = means_Cate,
    pvals = pvals_Cate,
    standard_scale = TRUE
)


### 4. Advanced Interaction Grouping
# You can also filter interactions by Gene Families (e.g., Coinhibitory, Costimulatory) or visualize specific pairs using plot_cpdb2.

ktplots::plot_cpdb(
    scdata = sce,
    cell_type1 = "Naive_B_cells",
    cell_type2 = "Classical_Monocytes",
    celltype_key = "scType_CPDB",
    means = means_Cate,
    pvals = pvals_Cate,
    gene_family = c("Coinhibitory", "Costimulatory"),
    standard_scale = TRUE,
    cluster_rows = FALSE
)



### Common Troubleshooting
# * Column Mismatch: If the plot is empty, check if `unique(sce$scType_CPDB)` exactly matches the column headers in your `pvals_Cate` file.
# * Gene Names: CellphoneDB uses Gene Symbols. If your data uses Ensembl IDs, the plot will not find any matches.
# * Data Scale: Using `standard_scale = TRUE` is recommended to make the dot sizes comparable across different interactions.
