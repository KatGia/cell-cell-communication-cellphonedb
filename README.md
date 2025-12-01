# cell-cell-communication-CellPhoneDB


This repository provides a simple workflow for performing **cell–cell communication analysis** using **CellPhoneDB** in **Python**, followed by **graphical visualization** of the results in **R**.

CellPhoneDB is a widely used tool that infers ligand–receptor interactions from single-cell RNA-seq data.
Here, the Python implementation is used to run the statistical analysis, while R is used for generating dot plots, heatmaps, and other visualizations.

---

## 📂 Workflow Overview

The analysis is divided into two phases:

1.  **Python (`src/run_analysis.py`):** Executes the statistical analysis, handles data reading (e.g., `.h5ad`), and generates raw result files (`.txt`).
2.  **R:** Loads the generated `.txt` files and creates publication-ready visualizations (Dot Plots, Heatmaps).

---

## 🧪 Requirements

* **Python ≥ 3.8** (Recommended for stability packages)
* **R ≥ 4.0**
* **CellPhoneDB** (Python package)

### 📌 R Packages for Visualization (To be added)

* `ktplots` (R)
* `ggplot2` (or `ktplotspy` in Python)

---

## 🚀 Getting Started: Environment Setup

We recommend using **Conda** for environment management.

### 1. Create and Activate the Environment

Create an environment named `cpdb_env` with the required Python version:

```bash
# 1. Ensure no other conda environment is active
conda deactivate 
# (You may need to repeat 'conda deactivate' until prompt is clean)

# 2. Create the new environment
conda create -n cpdb_env python=3.8

# 3. Activate the new environment
conda activate cpdb_env

```
### 2. Install CellPhoneDB Package
With the environment (cpdb_env) active, install the CellPhoneDB Python package:

```bash 
pip install cellphonedb
```
### 3. Install CellPhoneDB
With the environment (cpdb_env) active, install the CellPhoneDB Python package:

```bash
pip install cellphonedb
```

### 4. Download the Database
Download the required version of the database (v5.0.0 is used in the script) directly into the data/ folder:
```bash
# Ensure you are running this from the project's root directory:
python -c "from cellphonedb.utils import db_utils; db_utils.download_database(target_dir='data', cpdb_version='v5.0.0')"
```

(The script src/run_analysis.py expects the database file to be named cellphonedb.zip or the explicit version name used in the script).

---

## 📘 License

This project is released under the MIT License.

---

## 🙌 Acknowledgments

CellPhoneDB: [https://www.cellphonedb.org](https://www.cellphonedb.org)
Original publication by Efremova et al.

--- 

## ℹ️ Project Note

This repository is a work in progress and aims to help students transition smoothly between Python and R when working with CellPhoneDB for cell–cell communication analysis.
It is intended to facilitate the use of CellPhoneDB for students. CellPhoneDB is developed by Efremova et al., 2020. This project does not modify the original software.
