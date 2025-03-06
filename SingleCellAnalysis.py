#!/usr/bin/env python
# coding: utf-8

# In[2]:


pip install scanpy pandas numpy seaborn matplotlib


# In[2]:


import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Step 1: Load the PBMC 3k dataset from 10X Genomics
adata = sc.datasets.pbmc3k()


# In[3]:


# Step 2: Preprocessing
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, inplace=True)


# In[4]:


# Step 3: Normalization and Log Transformation
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]


# In[5]:


# Step 4: PCA and UMAP for Dimensionality Reduction
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)


# In[6]:


pip install igraph


# In[7]:


pip install leidenalg


# In[7]:


# Step 5: Clustering
sc.tl.leiden(adata, resolution=0.5, flavor="igraph", directed=False, n_iterations=2)


# In[8]:


# Step 6: Find marker genes
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
markers = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(5)


# In[9]:


# Step 7: Visualization
sc.pl.umap(adata, color=['leiden'], save='_clusters.png')
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, save='_marker_genes.png')


# In[10]:


# Step 8: Save processed data
adata.write('pbmc3k_processed.h5ad')


# In[11]:


# Display first few marker genes
print("Top marker genes:")
print(markers)

