# ==============================================================================
# CellPhoneDB v5 Input Preparation & Downstream Visualization (Seurat -> Python)
# Author: Caterina Giachino
# Description: Exports normalized counts (.h5ad), cell metadata, cluster DEGs,
#              and inferred active transcription factors (DoRothEA + VIPER)
#              from a Seurat object for CellPhoneDB v5 cell-cell communication
#              analysis, followed by dot plot visualization via ktplots.
# ==============================================================================

# 1. Libraries -----------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(zellkonverter)
  library(decoupleR)
  library(dplyr)
  library(ktplots)
})

# 2. Configuration & Paths -----------------------------------------------------
input_seurat_path <- "data/pbmc_normalized.rds" # Update with your relative path
output_dir        <- "results/cellphonedb_input"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Define the metadata column containing cell type annotations
# Note: Use annotated labels (e.g., "SingleR.labels" or "scType_classification")
# if plotting specific populations like "Myeloid" downstream.
celltype_col <- "SingleR.labels" 

# Load Seurat object
pbmc_normalized <- readRDS(input_seurat_path) # Or load() if .RData/.rDa
Idents(pbmc_normalized) <- celltype_col

# 3. Export Counts Matrix to AnnData (.h5ad) -----------------------------------
# Convert Seurat object to SingleCellExperiment
sce <- as.SingleCellExperiment(pbmc_normalized)

# Remove dimensional reductions to prevent HDF5/AnnData conversion conflicts
reducedDims(sce) <- NULL

# Write .h5ad using log-normalized counts ('logcounts') as primary matrix X
# (CellPhoneDB v5 expects normalized expression data)
writeH5AD(
  sce,
  file   = file.path(output_dir, "pbmc_normalized_counts.h5ad"),
  X_name = "logcounts"
)

# 4. Export Cell Metadata (META_FILE) ------------------------------------------
# Required format: 2 tab-separated columns ('Cell' and 'celltype')
cpdb_meta <- pbmc_normalized@meta.data %>%
  select(celltype = all_of(celltype_col)) %>%
  tibble::rownames_to_column(var = "Cell")

write.table(
  cpdb_meta,
  file      = file.path(output_dir, "pbmc_meta.txt"),
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

# 5. Compute & Export Differentially Expressed Genes (DEGS_FILE) ---------------
# Required for CellPhoneDB v5 DEG-based interaction filtering
pbmc_markers <- FindAllMarkers(
  pbmc_normalized,
  only.pos        = TRUE,
  min.pct         = 0.1,
  logfc.threshold = 0.25,
  test.use        = "wilcox"
)

# Format columns to match CellPhoneDB specifications (cluster in col 1, gene in col 2)
cpdb_degs <- pbmc_markers %>%
  filter(p_val_adj < 0.05) %>%
  select(cluster, gene, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)

write.table(
  cpdb_degs,
  file      = file.path(output_dir, "pbmc_degs.txt"),
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

# 6. Infer & Export Active Transcription Factors (ACTIVE_TFS_FILE) -------------
# Aggregate log-normalized expression at the cell-type level (pseudobulk)
pseudobulk_expr <- AggregateExpression(
  pbmc_normalized,
  assays   = "RNA",
  layer    = "data",
  group.by = celltype_col,
  fun      = "mean"
)$RNA

# Retrieve curated human TF-target regulons (DoRothEA confidence levels A, B, C)
dorothea_net <- decoupleR::get_dorothea(
  organism = "human",
  levels   = c("A", "B", "C")
)

# Run VIPER to estimate TF Normalized Enrichment Scores (NES) per cell type
viper_scores <- decoupleR::run_viper(
  mat     = as.matrix(pseudobulk_expr),
  net     = dorothea_net,
  .source = "source",
  .target = "target",
  .mor    = "mor",
  minsize = 5
)

# Filter significantly active/repressed TFs (|NES| > 1.5) and format for CPDB
cpdb_active_tfs <- viper_scores %>%
  filter(abs(score) > 1.5) %>%
  select(celltype = condition, active_TF = source)

write.table(
  cpdb_active_tfs,
  file      = file.path(output_dir, "pbmc_active_tfs.txt"),
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

message("CellPhoneDB v5 input files successfully generated in: ", output_dir)

# ==============================================================================
# 7. Downstream Visualization (Run after executing CellPhoneDB v5 in Python)
# ==============================================================================
cpdb_out_dir <- "results/cellphonedb_output"

# Load CellPhoneDB statistical/DEGs output tables
means_df   <- read.delim(file.path(cpdb_out_dir, "statistical_analysis_means.txt"), check.names = FALSE)
pvalues_df <- read.delim(file.path(cpdb_out_dir, "statistical_analysis_pvalues.txt"), check.names = FALSE)

# Generate Dot Plot for significant interactions involving Myeloid cells
pdf(file.path(cpdb_out_dir, "cpdb_myeloid_interactions_dotplot.pdf"), width = 10, height = 8)
plot_cpdb(
  scdata                = pbmc_normalized,
  cell_type1            = "Myeloid",
  cell_type2            = ".", # All other cell types
  celltype_key          = celltype_col,
  means                 = means_df,
  pvals                 = pvalues_df,
  keep_significant_only = TRUE,
  title                 = "Significant Cell-Cell Interactions: Myeloid vs. All"
)
dev.off()
