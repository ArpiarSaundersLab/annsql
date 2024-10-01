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
