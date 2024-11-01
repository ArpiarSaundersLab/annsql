from AnnSQL import AnnSQL
from AnnSQL.MakeDb import MakeDb
import scanpy as sc
import os
from memory_profiler import memory_usage

#load the scanpy object
#adata = sc.datasets.pbmc3k()

#load in backed mode for testing
adata = sc.read_h5ad("../data/pbmc3k_raw.h5ad", backed="r")

#this delete command is for testing purposes only. Remove this line in production 
if os.path.exists("../db/pbmc3k.asql"):
	os.remove("../db/pbmc3k.asql")
if os.path.exists("../db/pbmc3k.asql.wal"):
	os.remove("../db/pbmc3k.asql.wal")

def wrapper():
	MakeDb(adata=adata, 
			db_name="pbmc3k", 
			db_path="../db", 
			convenience_view=False,
			chunk_size=1000,
			make_buffer_file=False)

mem_usage = memory_usage(wrapper)
print(f"Memory usage: {max(mem_usage) - min(mem_usage)} MB")
