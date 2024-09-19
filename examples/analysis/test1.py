#import libraries
import scanpy as sc
from AnnSQL import AnnSQL
import time
import pandas as pd
from MakeDb import MakeDb

# #read in the data from scanpy
adata = sc.datasets.pbmc68k_reduced()

#call the AnnSQL class
adata_sql = AnnSQL(db="../db/pbmc68k_reduced.asql")

#return all x
adata_sql.query("SELECT SUM(COLUMNS(*)) FROM (SELECT * EXCLUDE (cell_id) FROM X)")
adata_sql.query("SELECT SUM(cast(columns(*exclude(cell_id,)) as double)) FROM X")
adata_sql.query("SUMMARIZE X")

#get a lst of all the genes
genes = adata_sql.query("SELECT * FROM var_names")["gene"].tolist()
column_sum = " + ".join([f'"{gene}"' for gene in genes])
query = f"SELECT ({column_sum}) AS total_sum FROM X;"
adata_sql.query(query)

time_start = time.time()
adata_sql.query("SELECT corr(ITGB2,SSU72) as correlation FROM adata WHERE bulk_labels = 'Dendritic' AND (ITGB2 > 0 OR SSU72 >0)")
time_end = time.time()
print(f"Time taken: {time_end - time_start}")

time_start = time.time()
results = adata[
    ((adata[:, "ITGB2"].X > 0) | (adata[:, "SSU72"].X > 0)).flatten() & 
    (adata.obs["bulk_labels"] == 'Dendritic'),["ITGB2","SSU72"]
].to_df()
results.corr()
time_end = time.time()
print(f"Time taken: {time_end - time_start}")