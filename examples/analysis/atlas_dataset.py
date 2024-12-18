#import libraries
import scanpy as sc
from AnnSQL import AnnSQL
from AnnSQL.MakeDb import MakeDb
import time

#load the atlas dataset (4.4 million cells) in either memory or backed mode
adata = sc.read_h5ad("Macosko_Mouse_Atlas_Single_Nuclei.Use_Backed.h5ad", backed="r+")
#adata_mem = sc.read_h5ad("Macosko_Mouse_Atlas_Single_Nuclei.Use_Backed_test.h5ad")

# #filter gene expressed > 0
# start_time = time.time()
# adata_mem[adata_mem[:, "ENSMUSG00000070880"].X > 0, "ENSMUSG00000070880"]
# print("--- %s seconds filter---" % (time.time() - start_time))

# #calculate metrics
# start_time = time.time()
# sc.pp.calculate_qc_metrics(adata_mem, inplace=True)
# print("--- %s seconds calculate_qc_metrics---" % (time.time() - start_time))

# #normalize the expression
# start_time = time.time()
# sc.pp.normalize_total(adata_mem, target_sum=1e4)
# print("--- %s seconds normalize_total---" % (time.time() - start_time))

# #log the expression
# start_time = time.time()
# sc.pp.log1p(adata_mem)
# print("--- %s seconds log1p---" % (time.time() - start_time))


#build a database to query later
MakeDb(adata=adata, 
		db_name="Macosko_Mouse_Atlas_Normalized_Ln", 
		db_path="../db/", 
		create_basic_indexes=True,
		layers=["X", "obs","var"])

#open the database
adata_sql = AnnSQL(db="../db/Macosko_Mouse_Atlas_Normalized_Ln.asql")

#calculate total counts
start_time = time.time()
adata_sql.calculate_total_counts(chunk_size=300)
print("--- %s seconds calculate_total_counts---" % (time.time() - start_time))

#normalize the expression
start_time = time.time()
adata_sql.expression_normalize(total_counts_per_cell=1e4, chunk_size=250, print_progress=True) 
print("--- %s seconds expression_normalize---" % (time.time() - start_time))

#calculate the log2 expression
start_time = time.time()
adata_sql.expression_log(log_type="LN", chunk_size=300, print_progress=True)
print("--- %s seconds expression_log---" % (time.time() - start_time))

#calculates gene variance (Experimental: sample variance per gene, includes Bessel's bias correction)
start_time = time.time()
adata_sql.calculate_variable_genes(chunk_size=800) 
print("--- %s seconds calculate_variable_genes---" % (time.time() - start_time))

#calculates the total gene counts
start_time = time.time()
adata_sql.calculate_gene_counts(chunk_size=800)
print("--- %s seconds calculate_gene_counts---" % (time.time() - start_time))