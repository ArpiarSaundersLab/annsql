#import libraries
import scanpy as sc
from MakeDb import MakeDb

#read in the data from scanpy
adata = sc.datasets.pbmc68k_reduced()

#build a database to query later
MakeDb(adata=adata, db_name="pbmc68k_reduced", db_path="db/")
