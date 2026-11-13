# ==============================================================================
# Reference-Based Cell Type Annotation with SingleR
# Author: Caterina Giachino
# Description: Annotates single-cell RNA-seq PBMC populations (including rare
#              cell subtypes) by comparing individual cell expression profiles
#              against the Human Primary Cell Atlas (celldex) reference,
#              validates prediction confidence, and maps labels onto UMAP.
# References:
#   - https://bioconductor.org/packages/release/bioc/vignettes/SingleR/inst/doc/SingleR.html
#   - https://swbioinf.github.io/scRNAseqInR_Doco/singler.html
# ==============================================================================

# 1. Libraries -----------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(SingleR)
  library(celldex)
  library(ggplot2)
})

# 2. Configuration & Paths -----------------------------------------------------
input_seurat_path  <- "data/pbmc_normalized.rds" # Update with your relative path
output_dir         <- "results/singler_annotation"
output_seurat_path <- "data/pbmc_annotated.rds"
min_cells_per_type <- 10 # Threshold to group ultra-rare calls into 'Other'

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_seurat_path), recursive = TRUE, showWarnings = FALSE)

# Load normalized Seurat object and convert to SingleCellExperiment
pbmc_normalized <- readRDS(input_seurat_path) # Or load() if using .RData/.rDa
sce <- as.SingleCellExperiment(pbmc_normalized)

# 3. Load Reference Dataset (celldex) ------------------------------------------
# HumanPrimaryCellAtlasData provides both broad (label.main) and granular
# (label.fine) transcriptomic profiles across 36+ human cell types.
# Other immune-focused alternatives: MonacoImmuneData(), BlueprintEncodeData()
ref_set <- celldex::HumanPrimaryCellAtlasData()

# 4. Run SingleR Cell-Level Classification -------------------------------------
# Cell-level annotation preserves rare populations (e.g., pDCs) that may not
# form independent variance-driven clusters in Seurat.
pred_main <- SingleR::SingleR(
  test      = sce,
  ref       = ref_set,
  labels    = ref_set$label.main,
  de.method = "classic" # Options: "classic", "wilcox", "t"
)

# Inspect label frequencies (raw vs. pruned low-confidence calls)
print(table(Raw = pred_main$labels))
print(table(Pruned = pred_main$pruned.labels, useNA = "ifany"))

# 5. Diagnostic Confidence Plots -----------------------------------------------
# Score heatmap: visualizes correlation scores across all reference labels per cell
pdf(file.path(output_dir, "singler_score_heatmap.pdf"), width = 10, height = 7)
plotScoreHeatmap(pred_main)
dev.off()

# Delta distribution: checks separation between assigned label and runner-up scores
pdf(file.path(output_dir, "singler_delta_distribution.pdf"), width = 10, height = 7)
plotDeltaDistribution(pred_main, ncol = 4)
dev.off()

# 6. Transfer Annotations to Seurat Metadata & Visualize -----------------------
# Retain labels exceeding the minimum cell count threshold; group rarer calls as 'Other'
lbls_keep <- table(pred_main$labels) > min_cells_per_type

pbmc_normalized$SingleR.labels <- ifelse(
  lbls_keep[pred_main$labels],
  pred_main$labels,
  "Other"
)

# Visualize SingleR annotations on UMAP
umap_singler <- DimPlot(
  pbmc_normalized,
  reduction = "umap",
  group.by  = "SingleR.labels",
  label     = TRUE,
  repel     = TRUE
) +
  ggtitle("PBMC Cell Type Annotation (SingleR - HPCA)") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  filename = file.path(output_dir, "umap_singler_annotations.pdf"),
  plot     = umap_singler,
  width    = 9,
  height   = 6
)

# 7. Save Annotated Object for Downstream Analyses (e.g., CellPhoneDB) ---------
saveRDS(pbmc_normalized, file = output_seurat_path)
saveRDS(pred_main, file = file.path(output_dir, "singler_predictions_main.rds"))

message("SingleR annotation complete. Annotated Seurat object saved to: ", output_seurat_path)
