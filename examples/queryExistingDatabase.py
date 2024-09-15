#import libraries
import scanpy as sc
from AnnSQL import AnnSQL

#call the AnnSQL class
adata_sql = AnnSQL(db="db/pbmc68k_reduced.asql")

#show all of the tables
adata_sql.show_tables()

#query a table
adata_sql.query("SELECT * FROM X")