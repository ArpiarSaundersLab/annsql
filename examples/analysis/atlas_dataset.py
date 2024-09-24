#import libraries
import scanpy as sc
from MakeDb import MakeDb
from AnnSQL import AnnSQL
import time
import pandas as pd

# #same as above, but using scanpy's read function
#adata = sc.read_h5ad("/home/kenny/Documents/OHSU/Projects/TAP/data/celltypist_models/chunked_approach/Macosko_Mouse_Atlas_Single_Nuclei.Use_Backed.h5ad", backed="r")

# #build a database to query later
#MakeDb(adata=adata, 
# 		db_name="Macosko_Mouse_Atlas_4mil", 
# 		db_path="../db/", 
# 		create_basic_indexes=True,
# 		layers=["X", "obs","var"])

# #call the AnnSQL class
adata_sql = AnnSQL(db="../db/Macosko_Mouse_Atlas.asql")

###################################################
# attempt to normalize the data
###################################################
#get the gene names
gene_names = adata_sql.query(f"Describe X")['column_name'][1:].values

#removes total_counts from the array & sets to 0
if "total_counts" in gene_names:
	adata_sql.update_query(f"UPDATE X SET total_counts = 0;")
	gene_names = gene_names[:-1] 
else:
	adata_sql.query(f"ALTER TABLE X ADD COLUMN total_counts FLOAT DEFAULT 0;")

#iterates gene_names in chunks
start_time = time.time()
chunk_size = 990 #Ddb limits arguments to 1k
for i in range(0, len(gene_names), chunk_size):
	chunk = gene_names[i:i+chunk_size]
	chunk = " + ".join(chunk) + " + total_counts"
	adata_sql.update_query(f"UPDATE X SET total_counts = ({chunk});")
end_time = time.time()
print("Total Counts Time: ", end_time-start_time)

#necessary to break dependency to change column types
adata_sql.query(f"DROP INDEX IF EXISTS idx_obs_cell_id;")
adata_sql.query(f"DROP INDEX IF EXISTS idx_X_cell_id;")

#normalize to 10k and log2
start_time = time.time()
for gene in gene_names:
	gene_start_time = time.time()
	adata_sql.update_query(f"UPDATE X SET {gene} = ROUND(LOG2((({gene} / total_counts) * 1e4) + 1e-4),5);")
	gene_end_time = time.time()
	print(f"{gene}: ", gene_end_time-gene_start_time)
end_time = time.time()
print("Normalize & Log2 Time: ", end_time-start_time)


# #Highly Variable Genes
# """
# SELECT 'gene_1' AS gene, VARIANCE(gene_1) AS variance FROM X
# UNION
# SELECT 'gene_2', VARIANCE(gene_2) FROM X
# UNION
# SELECT 'gene_3', VARIANCE(gene_3) FROM X
# ORDER BY variance DESC
# LIMIT 2000;
# """
