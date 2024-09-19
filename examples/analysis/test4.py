#import libraries
import scanpy as sc
from MakeDb import MakeDb
from AnnSQL import AnnSQL
import time
import pandas as pd

# #same as above, but using scanpy's read function
adata = sc.read_h5ad("/home/kenny/Documents/OHSU/Projects/TAP/data/celltypist_models/chunked_approach/Macosko_Mouse_Atlas_Single_Nuclei.Use_Backed.h5ad", backed="r")

# #build a database to query later
MakeDb(adata=adata, db_name="Macosko_Mouse_Atlas", db_path="../db/", layers=["X", "obs"])

# # #call the AnnSQL class
# adata_sql = AnnSQL(db="../db/Macosko_Mouse_Atlas.asql")

# #query a table
# start_time = time.time()
# adata_sql.query("SELECT ENSMUSG00000051951 FROM X WHERE ENSMUSG00000051951 > 0 LIMIT 5")
# end_time = time.time()
# print("Time taken: ", end_time-start_time)

# #do the same but with pandas
# start_time = time.time()
# pd.DataFrame(adata[adata[:,"ENSMUSG00000051951"].X > 0,"ENSMUSG00000051951"].X[:5], columns=["ENSMUSG00000051951"])
# end_time = time.time()
# print("Time taken: ", end_time-start_time)


# adata_sql.query("SELECT (ENSMUSG00000051951/1000*10000) FROM X")
# #adata_sql.update_query("UPDATE X SET ENSMUSG00000051951 = (ENSMUSG00000051951/1000*10000)")


#Highly Variable Genes
"""
SELECT 'gene_1' AS gene, VARIANCE(gene_1) AS variance FROM X
UNION
SELECT 'gene_2', VARIANCE(gene_2) FROM X
UNION
SELECT 'gene_3', VARIANCE(gene_3) FROM X
ORDER BY variance DESC
LIMIT 2000;
"""
