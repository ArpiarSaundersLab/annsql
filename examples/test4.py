#import libraries
import scanpy as sc
from MakeDb import MakeDb
from AnnSQL import AnnSQL

#same as above, but using scanpy's read function
adata = sc.read_h5ad("/home/kenny/Documents/OHSU/Projects/TAP/data/celltypist_models/chunked_approach/Macosko_Mouse_Atlas_Single_Nuclei.Use_Backed.h5ad", backed="r")

#build a database to query later
MakeDb(adata=adata, db_name="Macosko_Mouse_Atlas", db_path="db/", layers=["X", "obs"])


# #call the AnnSQL class
# adata_sql = AnnSQL(db="db/Macosko_Mouse_Atlas.asql")

# #show all of the tables
# adata_sql.show_tables()

# #query a table
# adata_sql.update_query("UPDATE X SET ENSMUSG00000051951 = (ENSMUSG00000051951/1000*10000)")
