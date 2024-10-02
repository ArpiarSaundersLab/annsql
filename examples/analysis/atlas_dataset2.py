#import libraries
import scanpy as sc
from AnnSQL import AnnSQL
import time

# #load the atlas dataset (4.4 million cells)
# adata = sc.read_h5ad("/home/kenny/Documents/OHSU/Projects/TAP/data/celltypist_models/chunked_approach/Macosko_Mouse_Atlas_Single_Nuclei.Use_Backed.h5ad", backed="r")

# #build a database to query later
# MakeDb(adata=adata, 
# 		db_name="Macosko_Mouse_Atlas_4mil", 
# 		db_path="../db/", 
# 		create_basic_indexes=True,
# 		layers=["X", "obs","var"])


#open the database
adata_sql = AnnSQL(db="../db/Macosko_Mouse_Atlas.asql")

#calculate total counts
start_time = time.time()
adata_sql.calculate_total_counts(chunk_size=950)
print("--- %s seconds calculate_total_counts---" % (time.time() - start_time))

#calculates the total gene counts
start_time = time.time()
adata_sql.calculate_gene_counts(chunk_size=950)
print("--- %s seconds calculate_gene_counts---" % (time.time() - start_time))

#normalize the expression
start_time = time.time()
adata_sql.expression_normalize(chunk_size=300, print_progress=True) 
print("--- %s seconds expression_normalize---" % (time.time() - start_time))



#calculate the log2 expression
start_time = time.time()
adata_sql.expression_log(log_type="LOG2", chunk_size=950)
print("--- %s seconds expression_log---" % (time.time() - start_time))

#calculates gene variance (Experimental: sample variance per gene, includes Bessel's bias correction)
start_time = time.time()
adata_sql.calculate_variable_genes(chunk_size=950) 
print("--- %s seconds calculate_variable_genes---" % (time.time() - start_time))
