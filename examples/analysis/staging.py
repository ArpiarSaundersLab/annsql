#import libraries
import scanpy as sc
from AnnSQL import AnnSQL
from MakeDb import MakeDb
import time

#load the atlas dataset (4.4 million cells)
adata = sc.read("/home/kenny/Documents/OHSU/Projects/TAP/data/celltypist_models/chunked_approach/Macosko_Mouse_Atlas_Single_Nuclei.Use_Backed.h5ad", backed='r')

#take the first 5 cells of the adata
subset = adata[:5].to_memory()

#make a annsql db from the subset
MakeDb(adata=subset,db_name="Macosko_Mouse_Atlas_Test",db_path="../db/")

#open the database
adata_sql = AnnSQL(db="../db/Macosko_Mouse_Atlas_Test.asql")

#calculate the gene counts
adata_sql.calculate_gene_counts(chunk_size=900, gene_field="gene_names_orig")

adata_sql.query("select * from var order by gene_counts DESC")

