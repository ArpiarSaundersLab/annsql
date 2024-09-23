#import libraries
import scanpy as sc
from MakeDb import MakeDb
from AnnSQL import AnnSQL
import time
import pandas as pd

# #same as above, but using scanpy's read function
#adata = sc.read_h5ad("/home/kenny/Documents/OHSU/Projects/TAP/data/celltypist_models/chunked_approach/Macosko_Mouse_Atlas_Single_Nuclei.Use_Backed.h5ad", backed="r")

# #build a database to query later
#MakeDb(adata=adata, db_name="Macosko_Mouse_Atlas_4mil", db_path="../db/", layers=["X", "obs","var"])

# #call the AnnSQL class
adata_sql = AnnSQL(db="../db/Macosko_Mouse_Atlas.asql")

#total counts per gene | Runtime: X seconds
start_time = time.time()
print(adata_sql.query("SELECT SUM(COLUMNS(*)) FROM (SELECT * EXCLUDE (cell_id) FROM X)"))
end_time = time.time()
print("Time taken: ", end_time-start_time)


# start_time = time.time()
# print(adata_sql.query("SELECT (ENSMUSG00000051951/1000*10000) FROM X"))
# end_time = time.time()
# print("Time taken: ", end_time-start_time)

#try normalizing the data 
gene_names = adata_sql.query("DESCRIBE X")["column_name"][0:990].values
gene_names = " + ".join(gene_names[1:])
start_time = time.time()
adata_sql.query(f"SELECT cell_id, ({gene_names}) as total_counts FROM X;")
end_time = time.time()
print("Time taken: ", end_time-start_time)



# adata_sql.query("SELECT cell_id, (ENSMUSG00000051951+ENSMUSG00000025900+ENSMUSG00000095041) FROM X")
# adata_sql.query("""
# SELECT 
#   cell_id, 
#   ROUND(ENSMUSG00000051951 / (ENSMUSG00000051951 + ENSMUSG00000025900 + ENSMUSG00000095041 + 1e-10)) AS total_counts, 
#   LOG2(CASE 
#          WHEN total_counts = 0 THEN 1e-10 
#          ELSE total_counts
#        END) AS log2_total_counts 
# FROM X
# """)


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
