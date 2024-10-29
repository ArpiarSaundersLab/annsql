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
adata = sc.read_h5ad("Macosko_Mouse_Atlas_Single_Nuclei.Use_Backed_test.h5ad", backed="r+")

#make a annsql db from the subset
MakeDb(adata=adata,db_name="Test",db_path="../db/")

# #open the database
asql = AnnSQL(db="../db/Test.asql")
asql.show_tables()