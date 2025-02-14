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
import os
warnings.filterwarnings('ignore')
from AnnSQL import AnnSQL
from AnnSQL.MakeDb import MakeDb


# #make the db
# filepath = "../data/splatter/data_5000_test.h5ad"
# adata = sc.read_h5ad(filepath)
# os.remove("../db/data_5000.asql")
# MakeDb(adata=adata, db_name="data_5000", db_path="../db")

#open the annsql database
asql = AnnSQL(db="../db/data_5000.asql")

#save a raw matrix layer
asql.save_raw()
#asql.raw_to_X() #can be now used if there's an oopsie

#view all tables (layers)
asql.show_tables()

#calculate total counts
asql.calculate_total_counts(chunk_size=750,print_progress=True)
sns.violinplot(x="total_counts", data=asql.query("SELECT total_counts FROM obs"))

#calculate gene counts
asql.calculate_gene_counts(chunk_size=750, print_progress=True)
sns.violinplot(data=asql.query("SELECT gene_counts FROM var"))

#filter by total umi counts
asql.filter_by_cell_counts(min_cell_count=2000, max_cell_count=15000)
sns.violinplot(x="total_counts", data=asql.query("SELECT total_counts FROM obs"))

#filter by gene counts
asql.filter_by_gene_counts(min_gene_counts=100, max_gene_counts=10000)
sns.violinplot(data=asql.query("SELECT gene_counts FROM var"))

#normalize the data
asql.expression_normalize(total_counts_per_cell=10000, chunk_size=750, print_progress=True)

#log transform the data (Ln)
asql.expression_log(chunk_size=750, print_progress=True, log_type="LN")

#calculate the variable genes
asql.calculate_variable_genes(chunk_size=750, print_progress=True, save_var_names=True, save_top_variable_genes=2000)

#determine principal components based on the top 2000 variable genes
start_time = time.time()
asql.calculate_pca(n_pcs=50, top_variable_genes=2000, chunk_size=750, print_progress=True, zero_center=False, max_cells_memory_threshold=1000)
print("Time to calculate PCA: ", time.time()-start_time)

#convert the long form table pcs to a matrix for viewing (memory intensive)
pca_scores = asql.return_pca_scores_matrix()

#plot the PCA
sns.scatterplot(x=pca_scores[0], y=pca_scores[1],size=0.5)
plt.xlabel("PCA1")
plt.ylabel("PCA2")



#TODO
#variance explained
#calculate the nearest neighbors
#leiden clustering
#umap
#plot the umap
#update the obs table with annotations
#differential expression
#plot the differential expression as volcano
