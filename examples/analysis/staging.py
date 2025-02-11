#import libraries
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
warnings.filterwarnings('ignore')

#parameters
n_pcs = 50
zero_center = False
top_variable_genes = 100
chunks = 25

#load the annsql object
asql = AnnSQL(db="../db/random/data_250000_test.asql")

#create table with structure (cell_id, gene, value)
asql.query_raw("DROP TABLE IF EXISTS X_standard; CREATE TABLE X_standard (cell_id STRING, gene STRING, value DOUBLE);")

#get all the gene names
genes = asql.query("SELECT gene_names FROM var")
genes = genes['gene_names'].tolist()
genes = genes[:top_variable_genes]

#build the query
if zero_center == True: 
	query_parts = [
		f"SELECT cell_id, '{gene}' AS gene, ( {gene} - AVG({gene}) OVER () ) / STDDEV({gene}) OVER () AS value FROM X"
		for gene in genes
	]
else:
	query_parts = [
		f"SELECT cell_id, '{gene}' AS gene, {gene} - AVG({gene}) OVER () AS value FROM X"
		for gene in genes
	]

#iterate query in a loop and print chunks
for i in range(0, len(query_parts), chunks):
	print("Chunk", i)
	full_query = "INSERT INTO X_standard " + " UNION ALL ".join(query_parts[i:i+chunks])
	asql.query_raw(full_query)

#create a new table for covariance
asql.query_raw("DROP TABLE IF EXISTS X_covariance; CREATE TABLE X_covariance (gene1 STRING, gene2 STRING, value DOUBLE);")

#convert the list to a string
genes_str = ', '.join([f"'{gene}'" for gene in genes])

#insert the covariance values
asql.query_raw(f"""
INSERT INTO X_covariance
SELECT 
    x.gene AS gene_1,
    y.gene AS gene_2,
    covar_samp(x.value, y.value) AS covariance
FROM X_standard x
JOIN X_standard y 
    ON x.cell_id = y.cell_id  -- Match cells
WHERE x.gene IN ({', '.join([f"'{gene}'" for gene in genes])}) 
  AND y.gene IN ({', '.join([f"'{gene}'" for gene in genes])}) 
   /* AND x.gene <= y.gene */
GROUP BY x.gene, y.gene;
""")

#take a look at the covariance
df = asql.query("SELECT * FROM X_covariance ORDER BY gene1, gene2")
df['gene1'] = df['gene1'].str.extract('(\d+)').astype(int) #sorts 
df['gene2'] = df['gene2'].str.extract('(\d+)').astype(int)
df = df.sort_values(by=['gene1', 'gene2'])

#pivots and create square covariance matrix. Smaller matrix that should be only genes x genes
cov_matrix = df.pivot(index="gene1", columns="gene2", values="value").fillna(0)

#convert to np
cov_matrix_np = cov_matrix.to_numpy()

#NOTE: Standardizing produces a matrix that mostly match perfectly to the covariance matrix produced by numpy. 
#A small amout are within 1e-4. (1.001002 vs 1.001001), but the fast majority of the values are the same.

#compute the eigenvectors and eigenvalues
eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix_np)
sorted_idx = np.argsort(eigenvalues)[::-1] #sort in descending order
eigenvalues_sorted = eigenvalues[sorted_idx]
eigenvectors_sorted = eigenvectors[:, sorted_idx]

#df for loadings: one row per gene for each of the top n_pcs PCs.
pc_loading_rows = []
for gene_idx, gene in enumerate(genes):
    for pc in range(n_pcs):
        pc_loading_rows.append({
            "gene": gene,
            "pc": pc,
            "loading": eigenvectors_sorted[gene_idx, pc]
        })

#save the pc loadings to a dataframe
pc_loadings_df = pd.DataFrame(pc_loading_rows)

#insert loadings to PC_loadings table 
asql.query_raw("DROP TABLE IF EXISTS PC_loadings;")
asql.query_raw("CREATE TABLE PC_loadings (gene STRING, pc INT, loading DOUBLE);")
values_list = []
for idx, row in pc_loadings_df.iterrows():
    gene_val = row["gene"].replace("'", "''")
    values_list.append(f"('{gene_val}', {int(row['pc'])}, {float(row['loading'])})")

values_str = ", ".join(values_list)
insert_sql = f"INSERT INTO PC_loadings (gene, pc, loading) VALUES {values_str};"
asql.query_raw(insert_sql)

#dot product entirely in SQL. Weird, but it works!
dot_product_query = f"""
SELECT 
    X.cell_id,
    P.pc,
    SUM(X.value * P.loading) AS pc_score
FROM X_standard X
JOIN PC_loadings P ON X.gene = P.gene
WHERE X.gene IN ({', '.join([f"'{gene}'" for gene in genes])})
GROUP BY X.cell_id, P.pc
ORDER BY X.cell_id, P.pc;
"""

#insert the dot product scores
asql.query_raw("DROP TABLE IF EXISTS PC_scores; CREATE TABLE PC_scores (cell_id STRING, pc INT, pc_score DOUBLE);")
asql.query_raw(f"INSERT INTO PC_scores {dot_product_query}")

# Execute the query to get the PC scores.
pc_scores_df = asql.query("SELECT * FROM PC_scores")
pc_scores_matrix_df = pc_scores_df.pivot(index="cell_id", columns="pc", values="pc_score").fillna(0)


#open anndata object and compare pca to annsql
adata = sc.read_h5ad("../data/random/data_1000_test.h5ad")
adata = adata[:, genes[:top_variable_genes]] #filter to top_variable_genes
sc.pp.pca(adata, n_comps=n_pcs)
pd.DataFrame(adata.obsm["X_pca"])
pd.DataFrame(pc_scores_matrix_df)