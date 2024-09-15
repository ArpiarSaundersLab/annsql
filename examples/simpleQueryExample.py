#import libraries
from AnnSQL import AnnSQL
import scanpy as sc

#read in the data from scanpy
adata = sc.datasets.pbmc68k_reduced()

#call the AnnSQL class
adata_sql = AnnSQL(adata=adata)

#show all of the tables
adata_sql.show_tables()

#query the X table
adata_sql.query("SELECT * FROM X LIMIT 5")