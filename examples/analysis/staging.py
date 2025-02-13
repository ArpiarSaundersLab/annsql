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


#make the db
filepath = "../data/splatter/data_1000.h5ad"
adata = sc.read_h5ad(filepath)
os.remove("../db/data_1000.asql")
MakeDb(adata=adata, db_name="data_1000", db_path="../db")

#open the annsql database
asql = AnnSQL(db="../db/data_1000.asql")

#save a raw matrix layer
asql.save_raw()
#asql.raw_to_X() #can be now used if there's an oopsie

#view all tables (layers)
asql.show_tables()

#calculate total counts
asql.calculate_total_counts(chunk_size=500,print_progress=True)
sns.violinplot(x="total_counts", data=asql.query("SELECT total_counts FROM obs"))

#calculate gene counts
asql.calculate_gene_counts(chunk_size=500, print_progress=True)
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
asql.expression_log(chunk_size=500, print_progress=True, log_type="LN")

#calculate the variable genes
asql.calculate_variable_genes(chunk_size=500, print_progress=True, save_var_names=True, save_top_variable_genes=2000)

#determine principal components based on the top 2000 variable genes
asql.calculate_pca(n_pcs=50, top_variable_genes=2000, chunk_size=500, print_progress=True, zero_center=False)

#convert the long form table pcs to a matrix for viewing (memory intensive)
pca_scores = asql.return_pca_scores_matrix()

#take a look at the main tables and the data
asql.show_tables()
asql.query("SELECT * FROM obs")
asql.query("SELECT * FROM var")
asql.query("SELECT * FROM X")

#plot the PCA
sns.scatterplot(x=pca_scores[0], y=pca_scores[1],size=0.5)


##########################################################################################
#compare the PCA to the scanpy PCA
##########################################################################################
filepath = "../data/splatter/data_1000_test.h5ad"
adata = sc.read_h5ad(filepath)

#qc the data
sc.pp.calculate_qc_metrics(adata, inplace=True)

#filter the data
sc.pp.filter_cells(adata, min_counts=2000)
sc.pp.filter_cells(adata, max_counts=15000)
sc.pp.filter_genes(adata, min_counts=100)
sc.pp.filter_genes(adata, max_counts=10000)

#normalize the data
sc.pp.normalize_total(adata, target_sum=10000)
sc.pp.log1p(adata)

#get all the genes that are highly variable from annsql
genes = asql.query("SELECT gene_names FROM var ")

#filter by the same genes in the genes	table
adata = adata[:,genes["gene_names"]]

#calculate the PCA
sc.tl.pca(adata, n_comps=50, zero_center=False, svd_solver=None)

#compare the PCA
sns.scatterplot(x=adata.obsm["X_pca"][:,0], y=adata.obsm["X_pca"][:,1],size=0.5)
