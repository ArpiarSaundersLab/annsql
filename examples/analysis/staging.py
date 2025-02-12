import scanpy as sc
import pandas as pd
import numpy as np
from AnnSQL import AnnSQL
from AnnSQL.MakeDb import MakeDb
from statsmodels.stats.multitest import multipletests
import time
import gc
import duckdb
import polars as pl
import seaborn as sns
import matplotlib.pyplot as plt
from numpy.linalg import eig
import warnings
from memory_profiler import memory_usage
warnings.filterwarnings('ignore')


# filepath = "../data/splatter/data_5000.h5ad"
# adata = sc.read_h5ad(filepath)
# MakeDb(adata=adata, db_name="data_5000", db_path="../db")


asql = AnnSQL(db="../db/data_5000.asql")
asql.save_raw()
asql.show_tables()
asql.calculate_total_counts(chunk_size=500,print_progress=True)
sns.violinplot(x="total_counts", data=asql.query("SELECT total_counts FROM obs"))
asql.filter_by_cell_counts(min_cell_count=2000, max_cell_count=15000)
sns.violinplot(x="total_counts", data=asql.query("SELECT total_counts FROM obs"))
asql.calculate_gene_counts(chunk_size=500, print_progress=True)
sns.violinplot(data=asql.query("SELECT gene_counts FROM var"))
asql.filter_by_gene_counts(min_gene_counts=100, max_gene_counts=10000)
sns.violinplot(data=asql.query("SELECT gene_counts FROM var"))
asql.calculate_variable_genes(chunk_size=500, print_progress=True, save_var_names=True, save_top_variable_genes=2000)

asql.raw_to_X()
asql.query("SELECT * FROM var")

asql.calculate_pca(n_pcs=25, top_variable_genes=2000, chunk_size=500, print_progress=True, zero_center=False)
pca_scores = asql.return_pca_scores_matrix()

pca_scores

#plot the PCA
fig, ax = plt.subplots()
sns.scatterplot(x=pca_scores[0], y=pca_scores[1],size=0.5, ax=ax)
ax.set_xlabel("PC1")
ax.set_ylabel("PC2")
plt.show()