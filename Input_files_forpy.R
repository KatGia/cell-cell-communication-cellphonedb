## R Tutorial: Preparing Input Files from Seurat
# To run the analysis in Python, you first need to export your Seurat Object into two specific files: the Counts Matrix and the Metadata Table.

### 1. Requirements
library(Seurat)
library(dplyr)
library(tibble)


### 2. Generate the Metadata File
# The metadata file must link each cell barcode to its annotated cell type.
# Crucial: Avoid spaces and special characters (like +) in cell type names.


# 1. Extract metadata
metadata <- pbmc_normalized@meta.data

# 2. Create the data frame (Barcode and Cell Type)
meta_df <- data.frame(
  Cell = rownames(metadata),
  cell_type = metadata$scType_classification # Use your cluster/annotation column
)

# 3. Clean names: Replace spaces and '+' with underscores
meta_df$cell_type <- gsub(" ", "_", meta_df$cell_type)
meta_df$cell_type <- gsub("\\+", "", meta_df$cell_type)

# 4. Save as .txt (Tab-separated, NO headers, NO quotes)
write.table(meta_df, 
            file = "cellphonedb_meta_2026_2.txt", 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)



### 3. Generate the Counts Matrix
# CellphoneDB requires **Log-Normalized** data (not raw counts and not z-scores).


# 1. Get Log-Normalized data
counts_matrix <- GetAssayData(pbmc_normalized, assay = "RNA", slot = "data")

# 2. Check: Max value should be low (e.g., between 5 and 15)
# If it's > 500, you are using raw counts (WRONG!)
print(max(counts_matrix))

# 3. Convert to data frame and clean gene names (remove .1, .2 suffixes)
counts_df <- as.data.frame(as.matrix(counts_matrix)) %>%
  rownames_to_column("Gene") %>%
  mutate(Gene = sub("\\..*$", "", Gene))

# 4. Handle Duplicate Genes (Aggregate by Mean)
# Required because CellphoneDB needs a unique value per gene/cell
counts_final <- counts_df %>%
  group_by(Gene) %>%
  summarise(across(everything(), mean)) %>%
  ungroup()

# 5. Prepare for export (Genes as rownames)
counts_matrix_ready <- as.data.frame(counts_final[, -1])
rownames(counts_matrix_ready) <- counts_final$Gene

# 6. Save with the "empty top-left cell" format
# col.names = NA is fundamental for proper alignment in Python/CellphoneDB
write.table(counts_matrix_ready, 
            file = "cellphonedb_counts_2026_final.txt", 
            sep = "\t", 
            quote = FALSE, 
            row.names = TRUE, 
            col.names = NA) 

# NOTE: 
###  Quality Check before Python
# Before running the Python script, verify these points:
# * Gene Names: Ensure they are HGNC symbols (e.g., CD3D, not ENSG000..).
# * Alignment: The number of columns in the Counts file (excluding the gene column) must exactly match the number of rows in the Meta file.
# * Special Characters: Check that your cell types are named like CD8_T_cell and not CD8+ T cell.

### Pro Tip: Using .h5ad
# For very large datasets, exporting to .txt can create huge files (several GBs). In those cases, it is recommended to use the anndata library in R or SeuratDisk to convert your object directly to .h5ad format, which is much faster and lighter for CellphoneDB to read.
