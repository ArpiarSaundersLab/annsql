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
asql.plot_umap(color_by="leiden_clusters")

#create cell type annotations
cell_types = {"0": "Microglia","1":"Inhibitory","2":"OPCs","3":"Astrocytes","4":"Excitatory"}	

#add the cell types to the obs table
asql.add_observations(obs_key="cell_types", obs_values=cell_types, match_on="leiden_clusters")

#plot the cell types
asql.plot_umap(color_by="cell_types")

#diff exp testing
obs_key = "cell_types"
group1_value = "Excitatory"
group2_value = "Inhibitory"
genes = asql.query("SELECT gene_names FROM var ORDER BY variance DESC LIMIT 2000;")["gene_names"].values.tolist()

start_time = time.time()
asql.close_db()
asql.open_db()
asql.conn.create_function("t_cdf", AnnSQL.t_cdf, [float, float], float)
for gene in genes:
	results = asql.conn.execute(f"""
	WITH stats AS (
		SELECT 
			obs.cell_types AS group_label,
			AVG(X.{gene})  AS mean_1,
			COUNT(X.{gene}) AS count_1,
			var_samp(X.{gene}) AS variance_1
		FROM X
		INNER JOIN obs ON X.cell_id = obs.cell_id 
		WHERE {obs_key} IN ('{group1_value}', '{group2_value}')
		GROUP BY obs.cell_types
	),
	calc AS (
		SELECT
			(s1.mean_1 - s2.mean_1) / SQRT(s1.variance_1/s1.count_1 + s2.variance_1/s2.count_1) AS tstat,
			LOG(s1.mean_1 / s2.mean_1) / LOG(2) AS logfc,
			POWER(s1.variance_1/s1.count_1 + s2.variance_1/s2.count_1, 2) /
			(
				POWER(s1.variance_1/s1.count_1, 2)/(s1.count_1 - 1) +
				POWER(s2.variance_1/s2.count_1, 2)/(s2.count_1 - 1)
			) AS df,
			2.0 * (1.0 - t_cdf(ABS(tstat), df)) AS pval
		FROM stats s1
		CROSS JOIN stats s2
		WHERE s1.group_label = '{group1_value}'
		AND s2.group_label = '{group2_value}'
	)
	SELECT * FROM calc ORDER BY pval DESC;
	""").df()
	#print(results)
asql.close_db()
end_time = time.time()
print(f"Time taken: {end_time - start_time}")

#check the ttest in python
group1 = asql.query(f"""
SELECT X.{gene} FROM X
INNER JOIN obs ON X.cell_id = obs.cell_id
WHERE obs.cell_types = '{group1_value}';
""")
group2 = asql.query(f"""
SELECT X.{gene} FROM X
INNER JOIN obs ON X.cell_id = obs.cell_id
WHERE obs.cell_types = '{group2_value}';
""")
from scipy.stats import ttest_ind
ttest_ind(group1, group2, equal_var=False)




#TODO
#differential expression
#plot the differential expression as volcano
