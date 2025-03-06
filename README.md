# Single-Cell RNA-Seq Analysis with Scanpy 🧬

## 📌 Overview
This project is a **self-learning exploration of single-cell RNA sequencing (scRNA-seq) analysis** using **Scanpy**, a Python-based toolkit. The dataset used is **PBMC 3k (Peripheral Blood Mononuclear Cells) from 10X Genomics**. The goal is to preprocess the data, perform dimensionality reduction, cluster the cells, and identify marker genes.

## 📂 Dataset
- **Source:** 10X Genomics PBMC 3k dataset
- **Format:** `.h5ad` (AnnData format used in Scanpy)
- **Description:** Contains single-cell RNA-seq data from **3,000 PBMCs**

## 🛠 Workflow
1. **Load and preprocess data**  
   - Filter low-quality cells and genes
   - Normalize and log-transform the data
   - Identify highly variable genes

2. **Dimensionality reduction**  
   - Perform **Principal Component Analysis (PCA)**
   - Compute **UMAP embeddings** for visualization

3. **Clustering**  
   - Construct a nearest-neighbor graph
   - Apply the **Leiden clustering algorithm**

4. **Marker gene identification**  
   - Perform **differential expression analysis** to find unique cell-type markers

5. **Visualization**  
   - Generate UMAP plots colored by cluster labels
   - Plot top marker genes for each cluster

## 📌 Requirements
Make sure you have the required Python packages installed:
```bash
pip install scanpy matplotlib seaborn pandas numpy
```

## 🚀 Running the Analysis
Execute the Jupyter Notebook step by step or run the script:
```python
python SingleCellAnalysis.py
```

## 📈 Results
- **UMAP Plot:** Shows clusters of cells based on gene expression profiles
- **Marker Gene Plot:** Identifies top genes that define each cluster
- **Processed Dataset:** Saved as `pbmc3k_processed.h5ad`

## References
- Scanpy Documentation: [https://scanpy.readthedocs.io](https://scanpy.readthedocs.io)
- 10X Genomics PBMC Dataset: [https://www.10xgenomics.com/resources/datasets](https://www.10xgenomics.com/resources/datasets)


