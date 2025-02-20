import scanpy as sc
import pandas as pd
import numpy as np
import time
import duckdb
import seaborn as sns
import matplotlib.pyplot as plt
from numpy.linalg import eig
import warnings
from memory_profiler import memory_usage
import os
warnings.filterwarnings('ignore')
from AnnSQL import AnnSQL
from AnnSQL.MakeDb import MakeDb
from scipy import stats


# #make the db
# filepath = "../data/splatter/data_1000_test.h5ad"
# adata = sc.read_h5ad(filepath)
# os.remove("../db/data_1000.asql")
# MakeDb(adata=adata, db_name="data_1000", db_path="../db")

#open the annsql database
asql = AnnSQL(db="../db/data_1000.asql")

#save a raw matrix layer
asql.save_raw()
#asql.raw_to_X() #can be now used if there's an oopsie

#view all tables (layers)
asql.show_tables()

#calculate total counts
asql.calculate_total_counts(chunk_size=750,print_progress=True)
asql.plot_total_counts()

#calculate gene counts
asql.calculate_gene_counts(chunk_size=750, print_progress=True)
asql.plot_gene_counts()

#filter by total umi counts
asql.filter_by_cell_counts(min_cell_count=2000, max_cell_count=15000)
asql.plot_total_counts()

#filter by gene counts
asql.filter_by_gene_counts(min_gene_counts=100, max_gene_counts=10000)
asql.plot_gene_counts()

#normalize the data
asql.expression_normalize(total_counts_per_cell=10000, chunk_size=750, print_progress=True)

#log transform the data (Ln)
asql.expression_log(chunk_size=750, print_progress=True, log_type="LN")

#calculate the variable genes
asql.calculate_variable_genes(chunk_size=750, print_progress=True, save_var_names=False)

#take a look at the highly variable
asql.plot_highly_variable_genes(top_variable_genes=2000)

#save the highly variable genes to X that seem correct from the plot above
asql.save_highly_variable_genes(top_variable_genes=2000)

#determine principal components based on the top 2000 variable genes
asql.calculate_pca(n_pcs=50, top_variable_genes=2000, chunk_size=750, print_progress=True, zero_center=False, max_cells_memory_threshold=1000)

#show the variance explained
asql.pca_variance_explained()

#plot the PCA
asql.plot_pca(PcX=1, PcY=2)

#calculate UMAP
asql.calculate_umap()

#plot the UMAP
asql.plot_umap()

#calculate leiden clusters
asql.calculate_leiden_clusters(resolution=1.0, n_neighbors=30)

#plot the leiden clusters
asql.plot_umap(color_by="leiden_clusters", legend_location="on data")

#calculate marker genes
asql.calculate_marker_genes(obs_key="leiden_clusters")

#plot the marker genes of 1 vs all
asql.plot_differential_expression(pvalue_threshold=0.05, logfc_threshold=0.5, group1="1", group2="ALL")

#plot marker genes for each cluster
asql.plot_marker_genes(obs_key="leiden_clusters", columns=3)

#return a list of marker genes
asql.get_marker_genes(obs_key="leiden_clusters", group="0")

#take a look at a marker gene
asql.plot_umap(color_by="gene_4601")

#create cell type annotations
cell_types = {"0": "Microglia","1":"Inhibitory","2":"OPCs","3":"Astrocytes","4":"Excitatory"}

#add the cell types to the obs table
asql.add_observations(obs_key="cell_type", obs_values=cell_types, match_on="leiden_clusters")

#plot the cell types
asql.plot_umap(color_by="cell_type")

#explicit calculate differential expression between groups
asql.calculate_differential_expression(obs_key="cell_type", group1_value="Excitatory", group2_value="Inhibitory", drop_table=False)

#plot the differential expression between groups
asql.plot_differential_expression(pvalue_threshold=0.05, logfc_threshold=0.5, group1="Excitatory", group2="Inhibitory")

