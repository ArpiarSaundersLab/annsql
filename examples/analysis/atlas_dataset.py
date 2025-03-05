#import libraries
import scanpy as sc
from memory_profiler import memory_usage
import time
from AnnSQL import AnnSQL
from AnnSQL.MakeDb import MakeDb
import pandas as pd

####################################################################################################
# AnnData/Scanpy runtime & memory analysis
####################################################################################################
#load the atlas dataset (4.4 million cells) in either memory or backed mode
#adata = sc.read_h5ad("Macosko_Mouse_Atlas_Single_Nuclei.Use_Backed.h5ad", backed="r+")
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

####################################################################################################

####################################################################################################

# #build a database to query later
# MakeDb(adata=adata, 
# 		db_name="Macosko_Mouse_Atlas_Proceesed", 
# 		db_path="../db/", 
# 		create_basic_indexes=True,
# 		layers=["X", "obs","var"])

#open the database
asql = AnnSQL(db="../db/Macosko_Mouse_Atlas_Processed.asql")
asql.show_tables()
asql.query("SELECT * FROM X_standard_wide LIMIT 1")


# Get top 2000 genes
genes_df = asql.query("SELECT gene_names_orig FROM var ORDER BY variance DESC")
genes = genes_df["gene_names_orig"].tolist()[:2000]

chunk_size = 200  
covariance_matrix = pd.DataFrame(index=genes, columns=genes, dtype=float)

for i, gene_1 in enumerate(genes[:2]): 
    start_time = time.time()
    gene_2_list = genes[i:]  #upper triangle
    for j in range(0, len(gene_2_list), chunk_size):
        gene_2_chunk = gene_2_list[j:j + chunk_size] 
        genes_clause = ", ".join([f"covar_samp({gene_1}, {gene_2}) AS cov_{gene_1}_{gene_2}" for gene_2 in gene_2_chunk])
        result = asql.query(f"SELECT {genes_clause} FROM X_standard_wide")
        for k, gene_2 in enumerate(gene_2_chunk):
            covariance_matrix.at[gene_1, gene_2] = result.iloc[0, k]

    print("--- %s seconds covariance ---" % (time.time() - start_time))


def test_filter():
	#simple filter
	asql.query("SELECT ENSMUSG00000070880 FROM X WHERE ENSMUSG00000070880 > 0")

def calculate_total_counts_memory_wrapper():
	#calculate total counts
	asql.calculate_total_counts(chunk_size=950, print_progress=True)

def calculate_gene_counts_memory_wrapper():
	#calculates the total gene counts
	asql.calculate_gene_counts(chunk_size=950, print_progress=True, gene_field="gene_names_orig")

def expression_normalize_memory_wrapper():
	#normalize the expression
	asql.expression_normalize(total_counts_per_cell=10000, chunk_size=250, print_progress=True) 

def expression_log_memory_wrapper():
	#calculate the log2 expression
	asql.expression_log(log_type="LN", chunk_size=250, print_progress=True)

def calculate_variable_genes_memory_wrapper():
	#calculates gene variance (Experimental: sample variance per gene, includes Bessel's bias correction)
	asql.calculate_variable_genes(chunk_size=250, print_progress=True, gene_field="gene_names_orig") 

def save_highly_variable_genes_memory_wrapper():
	#filter the highly variable genes
	asql.save_highly_variable_genes(top_variable_genes=2000, gene_field="gene_names_orig")

def calculate_pca_memory_wrapper():
	#pca
	asql.calculate_pca(n_pcs=50, top_variable_genes=2000, chunk_size=20, zero_center=False, print_progress=True, gene_field="gene_names_orig")

def calculate_umap_memory_wrapper():
	#umap
	asql.calculate_umap()

def calculate_leiden_clusters_memory_wrapper():
	#leiden
	asql.calculate_leiden_clusters(resolution=1, n_neighbors=5)

def calculate_differential_expression_memory_wrapper():
	#de of cluster 0 vs cluster 1
	asql.calculate_differential_expression(obs_key="leiden_clusters", group1_value="0", group2_value="1")


#create a df to store the results
# df = pd.DataFrame(columns=["function", "max_memory","runtime"])
# df.to_csv("../results/atlas_runtime_memory.csv", index=False)

# #simple filter
# start_time = time.time()
# result = memory_usage(test_filter)
# runtime = time.time() - start_time
# max_memory = max(result) - min(result)
# print("test_filter", max_memory, runtime)
# df = pd.read_csv("../results/atlas_runtime_memory.csv")
# df = pd.concat([df, pd.DataFrame([["test_filter", max_memory, runtime]], columns=["function", "max_memory","runtime"])])
# df.to_csv("../results/atlas_runtime_memory.csv", index=False)

# #run total counts 
# start_time = time.time()
# result = memory_usage(calculate_total_counts_memory_wrapper)
# runtime = time.time() - start_time
# max_memory = max(result) - min(result)
# print("calculate_total_counts_memory_wrapper", max_memory, runtime)
# df = pd.read_csv("../results/atlas_runtime_memory.csv")
# df = pd.concat([df, pd.DataFrame([["calculate_total_counts", max_memory, runtime]], columns=["function", "max_memory","runtime"])])
# df.to_csv("../results/atlas_runtime_memory.csv", index=False)

# #run calculate_gene_counts NOTE: only gene counts (excluded mean when profiling)
# start_time = time.time()
# result = memory_usage(calculate_gene_counts_memory_wrapper)
# runtime = time.time() - start_time
# max_memory = max(result) - min(result)
# print("calculate_gene_counts_memory_wrapper", max_memory, runtime)
# df = pd.read_csv("../results/atlas_runtime_memory.csv")
# df = pd.concat([df, pd.DataFrame([["calculate_gene_counts", max_memory, runtime]], columns=["function", "max_memory","runtime"])])
# df.to_csv("../results/atlas_runtime_memory.csv", index=False)

# #run expression_normalize
# start_time = time.time()
# result = memory_usage(expression_normalize_memory_wrapper)
# runtime = time.time() - start_time
# max_memory = max(result) - min(result)
# print("expression_normalize_memory_wrapper", max_memory, runtime)
# df = pd.read_csv("../results/atlas_runtime_memory.csv")
# df = pd.concat([df, pd.DataFrame([["expression_normalize", max_memory, runtime]], columns=["function", "max_memory","runtime"])])
# df.to_csv("../results/atlas_runtime_memory.csv", index=False)

# #run expression_log
# start_time = time.time()
# result = memory_usage(expression_log_memory_wrapper)
# runtime = time.time() - start_time
# max_memory = max(result) - min(result)
# print("expression_log_memory_wrapper", max_memory, runtime)
# df = pd.read_csv("../results/atlas_runtime_memory.csv")
# df = pd.concat([df, pd.DataFrame([["expression_log", max_memory, runtime]], columns=["function", "max_memory","runtime"])])
# df.to_csv("../results/atlas_runtime_memory.csv", index=False)

# #run calculate_variable_genes
# #NOTE fast, but log failed. Around 4 minutes
# start_time = time.time()
# result = memory_usage(calculate_variable_genes_memory_wrapper)
# runtime = time.time() - start_time
# max_memory = max(result) - min(result)
# print("calculate_variable_genes_memory_wrapper", max_memory, runtime)
# df = pd.read_csv("../results/atlas_runtime_memory.csv")
# df = pd.concat([df, pd.DataFrame([["calculate_variable_genes", max_memory, runtime]], columns=["function", "max_memory","runtime"])])

# #run save_highly_variable_genes
# start_time = time.time()
# result = memory_usage(save_highly_variable_genes_memory_wrapper)
# runtime = time.time() - start_time
# max_memory = max(result) - min(result)
# print("save_highly_variable_genes_memory_wrapper", max_memory, runtime)
# df = pd.read_csv("../results/atlas_runtime_memory.csv")
# df = pd.concat([df, pd.DataFrame([["save_highly_variable_genes", max_memory, runtime]], columns=["function", "max_memory","runtime"])])


#TODO Start here
#run calculate_pca
start_time = time.time()
result = memory_usage(calculate_pca_memory_wrapper)
runtime = time.time() - start_time
max_memory = max(result) - min(result)
print("calculate_pca_memory_wrapper", max_memory, runtime)
df = pd.read_csv("../results/atlas_runtime_memory.csv")
df = pd.concat([df, pd.DataFrame([["calculate_pca", max_memory, runtime]], columns=["function", "max_memory","runtime"])])

#run calculate_umap
start_time = time.time()
result = memory_usage(calculate_umap_memory_wrapper)
runtime = time.time() - start_time
max_memory = max(result) - min(result)
print("calculate_umap_memory_wrapper", max_memory, runtime)
df = pd.read_csv("../results/atlas_runtime_memory.csv")
df = pd.concat([df, pd.DataFrame([["calculate_umap", max_memory, runtime]], columns=["function", "max_memory","runtime"])])
df.to_csv("../results/atlas_runtime_memory.csv", index=False)

#run calculate_leiden_clusters
start_time = time.time()
result = memory_usage(calculate_leiden_clusters_memory_wrapper)
runtime = time.time() - start_time
max_memory = max(result) - min(result)
print("calculate_leiden_clusters_memory_wrapper", max_memory, runtime)
df = pd.read_csv("../results/atlas_runtime_memory.csv")
df = pd.concat([df, pd.DataFrame([["calculate_leiden_clusters", max_memory, runtime]], columns=["function", "max_memory","runtime"])])
df.to_csv("../results/atlas_runtime_memory.csv", index=False)

#run calculate_differential_expression
start_time = time.time()
result = memory_usage(calculate_differential_expression_memory_wrapper)
runtime = time.time() - start_time
max_memory = max(result) - min(result)
print("calculate_differential_expression_memory_wrapper", max_memory, runtime)
df = pd.read_csv("../results/atlas_runtime_memory.csv")
df = pd.concat([df, pd.DataFrame([["calculate_differential_expression", max_memory, runtime]], columns=["function", "max_memory","runtime"])])
df.to_csv("../results/atlas_runtime_memory.csv", index=False)

