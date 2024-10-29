#import libraries
import scanpy as sc
from AnnSQL import AnnSQL
from MakeDb import MakeDb
import time
import os

#delete the database if it exists
if os.path.exists("../db/Test.asql"):
	os.remove("../db/Test.asql")

#load the dataset
adata = sc.read("/home/kenny/Documents/OHSU/Projects/MouseMarmoset/adata_mouse_allages_allregions_forBroadSCPortal.h5ad", backed='r')

#take the first x cells of the adata
adata = adata[:10000].to_memory()

#make a annsql db from the subset
MakeDb(adata=adata,db_name="Test",db_path="../db/")

# #open the database
asql = AnnSQL(db="../db/Test.asql")
asql.show_tables()
