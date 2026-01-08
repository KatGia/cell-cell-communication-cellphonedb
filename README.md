# <img width="80" height="100" alt="image" src="https://github.com/user-attachments/assets/29743ad8-0f59-4e03-afc5-b4e88471b046" />  Analysis

**CellphoneDB** is a publicly available repository of curated receptors, ligands, and their interactions, coupled with a statistical framework to infer cell-cell communication from single-cell RNA sequencing (**scRNA-seq**) data.

---

## 🔬 Methodology Overview

| Method | Name | Brief Description | When to Use |
| :--- | :--- | :--- | :--- |
| **Method 1** | **Simple Analysis** | Calculates mean expression levels only. | For a quick look at present molecules without statistics. |
| **Method 2** | **Statistical Analysis** | Uses random shuffling to assign $P$-values. | The gold standard for discovering specific/significant interactions. |
| **Method 3** | **DEGs Analysis** | Uses a list of "important" genes (DEGs) to filter interactions. | When comparing conditions (e.g., Healthy vs. Disease). |

---

## 🛠 Detailed Comparison

| Feature | **METHOD 1: Simple** | **METHOD 2: Statistical** | **METHOD 3: DEGs** |
| :--- | :--- | :--- | :--- |
| **Main Goal** | Descriptive calculation of means. | Identify specific interactions via significance tests. | Identify interactions involving genes of interest. |
| **Calculation Logic** | Mean of Ligand-Receptor (L-R) components. | Random permutations (shuffling) of cell labels. | Filter based on user-provided DEG list. |
| **Complex Handling** | Uses the value of the lowest expressed member. | Same as Method 1, but used for null distribution. | At least one member must be in the DEG list. |
| **Expression Threshold**| Default **10%** (gene must be in >10% of cluster). | Default **10%**. If not met, interaction is ignored. | Mandatory **10%** for all interaction members. |
| **Significance** | **No** statistical test performed. | **Yes**. Based on real mean vs null distribution. | **No**. Relevance depends on your DEG list. |
| **Required Input** | Counts Matrix + Metadata. | Counts Matrix + Metadata. | Counts + Metadata + DEG File. |
| **Outputs** | `means.csv`, `deconvoluted.csv` | + `pvalues.csv`, `significant_means.csv` | `relevant_interactions.txt`, `significant_means.csv` |

---

## 🤔 Which method should I choose?

* **Choose Method 1** if you only need a fast overview of "what is there" without worrying about statistical specificity.
* **Choose Method 2** for standard publications. It is the most rigorous way to claim an interaction between **Type A** and **Type B** is significantly higher than in the rest of the tissue.
    > 💡 **Note:** For massive datasets, use **subsampling** (e.g., geometric sketching) to avoid long computation times during the 1,000 permutations.
* **Choose Method 3** if you have already performed a differential analysis (**Seurat/Scanpy**) and found interesting genes (e.g., up-regulated in a tumor). This "forces" CellphoneDB to focus only on those specific markers.

### 📊 Visualization with `ktplots`
**Method 2** and **Method 3** are required for visualization.

**Why doesn't Method 1 work?**
Method 1 does not generate significance (`pvalues.txt`) or relevance files. Since `ktplots` requires these values to filter and plot **Dot Plots** or **Heatmaps**, visualization is not possible for Method 1.



---

## 📁 Input Files

### 1. Mandatory Files
* **`cpdb_file_path` (`cellphonedb.zip`):** The core database containing curated biological information.
    * *Check for compatibility between the database and software version.*
* **`meta_file_path` (`metadata.tsv`):** A two-column file (**Cell Barcode** vs. **Cell Type**).
    * *Avoid spaces or special characters; use underscores (`_`).*
* **`counts_file_path` (`counts.h5ad`):** The data matrix.
    * Must be **normalized** (e.g., Log-Normalize) but **NOT scaled** (z-score).
    * **`.h5ad` (AnnData)** is highly recommended for better memory management.

### 2. Optional Files
* **`microenvs_file_path`:** Restricts interactions based on physical proximity.
* **`active_tf_path`:** Integrates gene regulation data (Transcription Factors).

---

## 📚 Resources & Documentation

For a deeper dive into the statistical theory and file preparation:
* [**CellphoneDB Results Documentation**](https://github.com/ventolab/CellphoneDB/blob/master/docs/RESULTS-DOCUMENTATION.md) – Deep dive into statistical methods.
* [**Input File Preparation Guide**](https://cellphonedb.readthedocs.io/en/latest/RESULTS-DOCUMENTATION.html#input-files) – Useful tips for creating files using **R**.

For plots: 
* [**In R library(ktplots)**](https://zktuong.github.io/ktplots/articles/vignette.html).
* [**In python (end of the web site)**](https://github.com/ventolab/CellphoneDB/blob/master/docs/RESULTS-DOCUMENTATION.md).

---
Now we'll use [**Seurat data**](https://satijalab.org/seurat/articles/pbmc3k_tutorial) (a dataset of Peripheral Blood Mononuclear Cells - PBMC) to check their communication. 
