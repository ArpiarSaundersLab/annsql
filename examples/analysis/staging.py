#import libraries
import scanpy as sc
import pandas as pd
from AnnSQL import AnnSQL
from AnnSQL.MakeDb import MakeDb
import time


# #open the database
# asql = AnnSQL(db="../db/Macosko_Mouse_Atlas_RAW.asql")

# #query the database
# asql.query("SELECT ENSMUSG00000070880 FROM X WHERE ENSMUSG00000070880 > 0 LIMIT 5")

# #anndata test
# adata = sc.read_h5ad("/home/kenny/Documents/OHSU/Projects/TAP/data/celltypist_models/chunked_approach/Macosko_Mouse_Atlas_Single_Nuclei.Use_Backed.h5ad", backed="r")

# #filter the ENSMUSG00000070880 the same way as the query but in scanpy
# start_time = time.time()
# adata[adata[:, "ENSMUSG00000070880"].X > 0, "ENSMUSG00000070880"]
# print("Filtering runtime no df:", time.time() - start_time)

# #from runtime analysis with df
# start_time = time.time()
# adata[adata[:, "ENSMUSG00000070880"].X > 0, "ENSMUSG00000070880"].to_df()
# print("Filtering runtime df:", time.time() - start_time)

# #from runtime analysis with df from pandas
# start_time = time.time()
# pd.DataFrame(adata[adata[:, "ENSMUSG00000070880"].X > 0, "ENSMUSG00000070880"], columns=["ENSMUSG00000070880"])
# print("Filtering runtime raw pandas df:", time.time() - start_time)