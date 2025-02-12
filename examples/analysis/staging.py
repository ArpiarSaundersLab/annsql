import scanpy as sc
import pandas as pd
import numpy as np
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
from AnnSQL import AnnSQL
from AnnSQL.MakeDb import MakeDb

#make the db
filepath = "../data/splatter/data_5000.h5ad"
adata = sc.read_h5ad(filepath)
MakeDb(adata=adata, db_name="data_5000", db_path="../db")

#open the annsql database
asql = AnnSQL(db="../db/data_5000.asql")

#save a raw matrix layer
asql.save_raw()
#asql.raw_to_X() #can be now used if there's an oopsie

#view all tables (layers)
asql.show_tables()

#calculate total counts
asql.calculate_total_counts(chunk_size=500,print_progress=True)
sns.violinplot(x="total_counts", data=asql.query("SELECT total_counts FROM obs"))

#filter by total umi counts
asql.filter_by_cell_counts(min_cell_count=2000, max_cell_count=15000)
sns.violinplot(x="total_counts", data=asql.query("SELECT total_counts FROM obs"))

#calculate gene counts
asql.calculate_gene_counts(chunk_size=500, print_progress=True)
sns.violinplot(data=asql.query("SELECT gene_counts FROM var"))

#filter by gene counts
asql.filter_by_gene_counts(min_gene_counts=100, max_gene_counts=10000)
sns.violinplot(data=asql.query("SELECT gene_counts FROM var"))

#normalize the data
asql.expression_normalize(total_counts_per_cell=10000, chunk_size=500, print_progress=True)

#log transform the data (Ln)
asql.expression_log(chunk_size=500, print_progress=True)

#calculate the variable genes
asql.calculate_variable_genes(chunk_size=500, print_progress=True, save_var_names=True, save_top_variable_genes=2000)

#determine principal components based on the top 2000 variable genes
asql.calculate_pca(n_pcs=25, top_variable_genes=2000, chunk_size=500, print_progress=True, zero_center=False)

#convert the long form table pcs to a matrix for viewing (memory intensive)
pca_scores = asql.return_pca_scores_matrix()
pca_scores

#plot the PCA
fig, ax = plt.subplots()
sns.scatterplot(x=pca_scores[0], y=pca_scores[1],size=0.5, ax=ax)
ax.set_xlabel("PC1")
ax.set_ylabel("PC2")
plt.show()

#take a look at the main tables and the data
asql.show_tables()
asql.query("SELECT * FROM obs LIMIT 10")
asql.query("SELECT * FROM var LIMIT 10")
asql.query("SELECT * FROM X LIMIT 10")

