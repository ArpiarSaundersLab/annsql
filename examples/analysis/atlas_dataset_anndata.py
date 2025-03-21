#import libraries
import scanpy as sc
from memory_profiler import memory_usage
import time
from AnnSQL import AnnSQL
from AnnSQL.MakeDb import MakeDb
import pandas as pd
import gc

####################################################################################################
# Run the AnnSQL runtime and memory of all major procedures
####################################################################################################

def filter_annsql_memory_backed_wrapper():
	#filter gene expressed > 0
	adata_backed[adata_backed[:, "ENSMUSG00000070880"].X > 0, "ENSMUSG00000070880"]

def filter_annsql_inmem_memory_wrapper():
	#filter gene expressed > 0
	adata_mem[adata_mem[:, "ENSMUSG00000070880"].X > 0, "ENSMUSG00000070880"]

def calculate_qc_metrics_memory_wrapper():
	#calculate metrics
	sc.pp.calculate_qc_metrics(adata_mem, inplace=True)

def calculate_qc_metrics_inmem_memory_wrapper():
	#calculate metrics
	sc.pp.calculate_qc_metrics(adata_mem, inplace=True)

def normalize_total_memory_wrapper():
	#normalize the expression
	sc.pp.normalize_total(adata_mem, target_sum=1e4)

####################################################################################################
# AnnData/Scanpy runtime and memory of all major procedures
# NOTE: all other procedures for adata failed except filtering on laptop and HPC
####################################################################################################

# #load the atlas dataset (4.4 million cells) in memory
# adata_backed = sc.read_h5ad("/home/exacloud/gscratch/SaundersLab/24-09_CellTypist/chunked_approach/Macosko_Mouse_Atlas_Single_Nuclei.Use_Backed.h5ad", backed="r")

# #simple filter in backed
# start_time = time.time()
# result = memory_usage(filter_annsql_memory_backed_wrapper)
# runtime = time.time() - start_time
# max_memory = max(result) - min(result)
# print("filter_annsql_memory_backed_wrapper", max_memory, runtime)
# df = pd.read_csv("atlas_profile_hpc.csv")
# df = pd.concat([df, pd.DataFrame([["filter_annsql_memory_backed_wrapper", max_memory, runtime]], columns=["function", "max_memory","runtime"])])
# df.to_csv("atlas_profile_hpc.csv", index=False)


#load the atlas dataset (4.4 million cells) in backed mode
adata_mem = sc.read_h5ad("/home/exacloud/gscratch/SaundersLab/24-09_CellTypist/chunked_approach/Macosko_Mouse_Atlas_Single_Nuclei.Use_Backed.h5ad")

#simple filter memory mode
start_time = time.time()
result = memory_usage(filter_annsql_inmem_memory_wrapper)
runtime = time.time() - start_time
max_memory = max(result) - min(result)
print("filter_annsql_memory_inmem_wrapper", max_memory, runtime)
df = pd.read_csv("atlas_profile_hpc.csv")
df = pd.concat([df, pd.DataFrame([["filter_annsql_memory_inmem_wrapper", max_memory, runtime]], columns=["function", "max_memory","runtime"])])
df.to_csv("atlas_profile_hpc.csv", index=False)

#simple filter memory mode
start_time = time.time()
result = memory_usage(filter_annsql_inmem_memory_wrapper)
runtime = time.time() - start_time
max_memory = max(result) - min(result)
print("filter_annsql_memory_inmem_wrapper", max_memory, runtime)
df = pd.read_csv("atlas_profile_hpc.csv")
df = pd.concat([df, pd.DataFrame([["filter_annsql_memory_inmem_wrapper", max_memory, runtime]], columns=["function", "max_memory","runtime"])])
df.to_csv("atlas_profile_hpc.csv", index=False)

#simple filter memory mode
start_time = time.time()
result = memory_usage(filter_annsql_inmem_memory_wrapper)
runtime = time.time() - start_time
max_memory = max(result) - min(result)
print("filter_annsql_memory_inmem_wrapper", max_memory, runtime)
df = pd.read_csv("atlas_profile_hpc.csv")
df = pd.concat([df, pd.DataFrame([["filter_annsql_memory_inmem_wrapper", max_memory, runtime]], columns=["function", "max_memory","runtime"])])
df.to_csv("atlas_profile_hpc.csv", index=False)

#simple filter memory mode
start_time = time.time()
result = memory_usage(filter_annsql_inmem_memory_wrapper)
runtime = time.time() - start_time
max_memory = max(result) - min(result)
print("filter_annsql_memory_inmem_wrapper", max_memory, runtime)
df = pd.read_csv("atlas_profile_hpc.csv")
df = pd.concat([df, pd.DataFrame([["filter_annsql_memory_inmem_wrapper", max_memory, runtime]], columns=["function", "max_memory","runtime"])])
df.to_csv("atlas_profile_hpc.csv", index=False)

#simple filter memory mode
start_time = time.time()
result = memory_usage(filter_annsql_inmem_memory_wrapper)
runtime = time.time() - start_time
max_memory = max(result) - min(result)
print("filter_annsql_memory_inmem_wrapper", max_memory, runtime)
df = pd.read_csv("atlas_profile_hpc.csv")
df = pd.concat([df, pd.DataFrame([["filter_annsql_memory_inmem_wrapper", max_memory, runtime]], columns=["function", "max_memory","runtime"])])
df.to_csv("atlas_profile_hpc.csv", index=False)
