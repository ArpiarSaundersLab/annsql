#import libraries
import scanpy as sc
from AnnSQL import AnnSQL
import time
import pandas as pd
from MakeDb import MakeDb
import os

#remove the file if it exists
if os.path.exists("../db/pbmc3k_raw.asql"):
	os.remove("../db/pbmc3k_raw.asql")

# load the data/pbmc68k_raw.h5ad file
adata = sc.read("../data/pbmc3k_raw.h5ad", backed="r+")

#make the database using backed mode
MakeDb(adata, db_path="../db/", db_name="pbmc3k_raw")

#call the AnnSQL class
adata_sql = AnnSQL(db="../db/pbmc3k_raw.asql")

#return all x
adata_sql.query("SELECT * FROM adata")

#show the adata index using the adata object
adata.obs





