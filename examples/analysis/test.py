import scanpy as sc
from MakeDb import MakeDb
from AnnSQL import AnnSQL

#read sample data
adata = sc.datasets.pbmc68k_reduced()

#build the AnnSQL database
MakeDb(adata=adata, db_name="pbmc3k_reduced", db_path="..db/")

#open the AnnSQL database
adata_sql = AnnSQL(db="..db/pbmc3k_reduced.asql")

#query the expression table
adata_sql.query("SELECT * FROM var_names")