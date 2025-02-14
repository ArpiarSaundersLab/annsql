import scanpy as sc
from AnnSQL import AnnSQL
import time
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from memory_profiler import memory_usage
import os
import gc
import sys

#set seed for reproducibility
np.random.seed(0)

#generate a list of 10 random numbers between 0 and 9999
gene_ids = np.random.randint(0, 9999, 20)
gene_ids_2 = np.random.randint(0, 9999, 20)

#take a look at the gene_ids and confirm they are same as in gene_ids.txt (they are)
print(f"gene_ids: {gene_ids}")
print(f"gene_ids_2: {gene_ids_2}")

def add_to_comparisons(comparisons, size, method, memory, filter_type, gene_id):
	comparisons = pd.concat([
		comparisons, 
		pd.DataFrame([{"size": size, "method": method, "memory": memory, "filter_type": filter_type, "gene_id": gene_id}])
	], ignore_index=False)
	comparisons.to_csv(file_name, index=False, mode='a', header=False)
	comparisons = pd.DataFrame(columns=["size", "method", "memory", "filter_type","gene_id"])

def does_record_exist(size, method, filter_type,gene_id):
	record_check = pd.read_csv(file_name)
	if len(record_check[(record_check['size'] == size) & (record_check['method'] == method) & (record_check['filter_type'] == filter_type) & (record_check['gene_id'] == gene_id)]) > 0:
		return True
	return False

#functions for each of the filters necessary for memory analysis
def load_anndata_in_memory():
	global adata_memory
	adata_memory = sc.read(data_path+"data_" + str(i) + ".h5ad")

def load_anndata_on_disk():
	global adata_disk
	adata_disk = sc.read(data_path+"data_" + str(i) + ".h5ad", backed="r")

def load_annsql_in_memory():
	global adata_sql_memory
	adata_sql_memory = AnnSQL(adata=data_path+"data_" + str(i) + ".h5ad", print_output=False)

def load_annsql_on_disk():
	global adata_sql_disk
	adata_sql_disk = AnnSQL(db=db_path+"data_" + str(i) + ".asql", print_output=False)

def annsql_filter_1():
	adata_sql_memory.query(f"SELECT gene_{gene_id} FROM X WHERE gene_{gene_id} > 0.5 LIMIT 5")

def anndata_filter_1():
	pd.DataFrame(adata_memory[adata_memory[:, "gene_"+str(gene_id)].X > 0.5, "gene_"+str(gene_id)].X[:5])

def annsql_filter_2():
	adata_sql_memory.query(f"SELECT gene_{gene_id} FROM X WHERE gene_{gene_id} > 0.5 AND gene_{gene_id_2} < 0.5 LIMIT 5")

def anndata_filter_2():
	pd.DataFrame(adata_memory[(adata_memory[:, "gene_"+str(gene_id)].X > 0.5) & (adata_memory[:, "gene_"+str(gene_id_2)].X < 0.5), "gene_"+str(gene_id)].X[:5])

def annsql_filter_3():
	adata_sql_memory.query(f"SELECT gene_{gene_id} FROM X WHERE gene_{gene_id} > 0.5 AND gene_{gene_id_2} < 0.5 AND gene_{gene_id_3} > 0.5 LIMIT 5")

def anndata_filter_3():
	pd.DataFrame(adata_memory[(adata_memory[:, "gene_"+str(gene_id)].X > 0.5) & (adata_memory[:, "gene_"+str(gene_id_2)].X < 0.5) & (adata_memory[:, "gene_"+str(gene_id_3)].X > 0.5), "gene_"+str(gene_id)].X[:5])

def annsql_filter_4():
	adata_sql_memory.query(f"SELECT gene_{gene_id} FROM X WHERE gene_{gene_id} > 0.5 AND gene_{gene_id_2} < 0.5 AND gene_{gene_id_3} > 0.5 AND gene_{gene_id_4} < 0.5 LIMIT 5")

def anndata_filter_4():
	pd.DataFrame(adata_memory[(adata_memory[:, "gene_"+str(gene_id)].X > 0.5) & (adata_memory[:, "gene_"+str(gene_id_2)].X < 0.5) & (adata_memory[:, "gene_"+str(gene_id_3)].X > 0.5) & (adata_memory[:, "gene_"+str(gene_id_4)].X < 0.5), "gene_"+str(gene_id)].X[:5])

def annsql_filter_5():
	adata_sql_memory.query(f"SELECT gene_{gene_id}, gene_{gene_id_2}, gene_{gene_id_3}, gene_{gene_id_4}, (gene_{gene_id_3}/gene_{gene_id_4}) as gene_3_4	FROM X WHERE (gene_{gene_id} > 0.5 AND gene_{gene_id_2} < 0.3) AND (gene_{gene_id_3}/gene_{gene_id_4} > 0 AND gene_{gene_id_3}/gene_{gene_id_4} < 0.5) LIMIT 5")

def anndata_filter_5():
	gene_1_vals = adata_memory[:, "gene_"+str(gene_id)].X.flatten()
	gene_2_vals = adata_memory[:, "gene_"+str(gene_id_2)].X.flatten()
	gene_3_vals = adata_memory[:, "gene_"+str(gene_id_3)].X.flatten()
	gene_4_vals = adata_memory[:, "gene_"+str(gene_id_4)].X.flatten()
	d_3_4 = gene_3_vals / gene_4_vals
	condition = (gene_1_vals > 0.5) & (gene_2_vals < 0.3) & (d_3_4 > 0) & (d_3_4 < 0.5)
	filtered_indices = np.where(condition)[0]
	limited_indices = filtered_indices[:5]
	pd.DataFrame({
		"gene_"+str(gene_id): gene_1_vals[limited_indices],
		"gene_"+str(gene_id_2): gene_2_vals[limited_indices],
		"gene_"+str(gene_id_3): gene_3_vals[limited_indices],
		"gene_"+str(gene_id_4): gene_4_vals[limited_indices],
		"gene_3_4": d_3_4[limited_indices]
	})

def annsql_filter_6():
	adata_sql_memory.query(f"SELECT obs.cell_type, AVG(gene_{gene_id}) as avg_gene_{gene_id} FROM X INNER JOIN obs ON X.cell_id = obs.cell_id GROUP BY obs.cell_type")

def anndata_filter_6():
	pd.DataFrame({
		'cell_type': adata_memory.obs['cell_type'],
		'gene_'+str(gene_id): adata_memory[:, 'gene_'+str(gene_id)].X.flatten()
	}).groupby('cell_type')['gene_'+str(gene_id)].mean().reset_index()
	
def annsql_filter_1_disk():
	adata_sql_disk.query(f"SELECT gene_{gene_id} FROM X WHERE gene_{gene_id} > 0.5 LIMIT 5")

def anndata_filter_1_disk():
	pd.DataFrame(adata_disk[adata_disk[:, "gene_"+str(gene_id)].X > 0.5, "gene_"+str(gene_id)].X[:5])

def annsql_filter_2_disk():
	adata_sql_disk.query(f"SELECT gene_{gene_id} FROM X WHERE gene_{gene_id} > 0.5 AND gene_{gene_id_2} < 0.5 LIMIT 5")

def anndata_filter_2_disk():
	pd.DataFrame(adata_disk[(adata_disk[:, "gene_"+str(gene_id)].X > 0.5) & (adata_disk[:, "gene_"+str(gene_id_2)].X < 0.5), "gene_"+str(gene_id)].X[:5])

def annsql_filter_3_disk():
	adata_sql_disk.query(f"SELECT gene_{gene_id} FROM X WHERE gene_{gene_id} > 0.5 AND gene_{gene_id_2} < 0.5 AND gene_{gene_id_3} > 0.5 LIMIT 5")

def anndata_filter_3_disk():
	pd.DataFrame(adata_disk[(adata_disk[:, "gene_"+str(gene_id)].X > 0.5) & (adata_disk[:, "gene_"+str(gene_id_2)].X < 0.5) & (adata_disk[:, "gene_"+str(gene_id_3)].X > 0.5), "gene_"+str(gene_id)].X[:5])

def annsql_filter_4_disk():
	adata_sql_disk.query(f"SELECT gene_{gene_id} FROM X WHERE gene_{gene_id} > 0.5 AND gene_{gene_id_2} < 0.5 AND gene_{gene_id_3} > 0.5 AND gene_{gene_id_4} < 0.5 LIMIT 5")

def anndata_filter_4_disk():
	pd.DataFrame(adata_disk[(adata_disk[:, "gene_"+str(gene_id)].X > 0.5) & (adata_disk[:, "gene_"+str(gene_id_2)].X < 0.5) & (adata_disk[:, "gene_"+str(gene_id_3)].X > 0.5) & (adata_disk[:, "gene_"+str(gene_id_4)].X < 0.5), "gene_"+str(gene_id)].X[:5])

def annsql_filter_5_disk():
	adata_sql_disk.query(f"SELECT gene_{gene_id}, gene_{gene_id_2}, gene_{gene_id_3}, gene_{gene_id_4}, (gene_{gene_id_3}/gene_{gene_id_4}) as gene_3_4 FROM X WHERE (gene_{gene_id} > 0.5 AND gene_{gene_id_2} < 0.3) AND (gene_{gene_id_3}/gene_{gene_id_4} > 0 AND gene_{gene_id_3}/gene_{gene_id_4} < 0.5) LIMIT 5")

def anndata_filter_5_disk():
	gene_1_vals = adata_disk[:, "gene_"+str(gene_id)].X.flatten()
	gene_2_vals = adata_disk[:, "gene_"+str(gene_id_2)].X.flatten()
	gene_3_vals = adata_disk[:, "gene_"+str(gene_id_3)].X.flatten()
	gene_4_vals = adata_disk[:, "gene_"+str(gene_id_4)].X.flatten()
	d_3_4 = gene_3_vals / gene_4_vals
	condition = (gene_1_vals > 0.5) & (gene_2_vals < 0.3) & (d_3_4 > 0) & (d_3_4 < 0.5)
	filtered_indices = np.where(condition)[0]
	limited_indices = filtered_indices[:5]
	pd.DataFrame({
		"gene_"+str(gene_id): gene_1_vals[limited_indices],
		"gene_"+str(gene_id_2): gene_2_vals[limited_indices],
		"gene_"+str(gene_id_3): gene_3_vals[limited_indices],
		"gene_"+str(gene_id_4): gene_4_vals[limited_indices],
		"gene_3_4": d_3_4[limited_indices]
	})



def annsql_filter_6_disk():
	adata_sql_disk.query(f"SELECT obs.cell_type, AVG(gene_{gene_id}) as avg_gene_{gene_id} FROM X INNER JOIN obs ON X.cell_id = obs.cell_id GROUP BY obs.cell_type")

def anndata_filter_6_disk():
	pd.DataFrame({
		'cell_type': adata_disk.obs['cell_type'],
		'gene_'+str(gene_id): adata_disk[:, 'gene_'+str(gene_id)].X.flatten()
	}).groupby('cell_type')['gene_'+str(gene_id)].mean().reset_index()



data_path = "../data/splatter/"
db_path = "../db/splatter/"
file_name = f"../results/comparisons_memory_splatter_Pythonb.csv"
in_memory_high_filter = 100001
dataset_sizes = [1000, 5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000, 75000, 100000, 250000]
dataset_sizes = [9999]


for i in dataset_sizes:

	if not os.path.exists(data_path+"data_" + str(i) + ".h5ad"):
		continue

	print(f"Running for {i}")

	#in-memory vs non-backed
	if i <= in_memory_high_filter:
		mem_usage = memory_usage(load_anndata_in_memory)
		load_anndata_mem_max_memory = max(mem_usage) - min(mem_usage)
		mem_usage = memory_usage(load_annsql_in_memory)
		load_annsql_mem_max_memory = max(mem_usage) - min(mem_usage)
		

	#on-disk vs backed
	mem_usage = memory_usage(load_anndata_on_disk)
	load_anndata_disk_max_memory = max(mem_usage) - min(mem_usage)
	mem_usage = memory_usage(load_annsql_on_disk)
	load_annsql_disk_max_memory = max(mem_usage) - min(mem_usage)

	for gene in gene_ids:

		global gene_id, gene_id_2, gene_id_3, gene_id_4
		gene_id = gene
		gene_id_2 = gene_ids_2[np.random.randint(0, 20)]
		gene_id_3 = gene_ids[np.random.randint(0, 20)]
		gene_id_4 = gene_ids_2[np.random.randint(0, 20)]


		#open the comparison file does not exist
		if not os.path.exists(file_name):
			comparisons = pd.DataFrame(columns=["size", "method", "memory", "filter_type","gene_id"])
			comparisons.to_csv(file_name, index=False)
		else:
			comparisons = pd.DataFrame(columns=["size", "method", "memory", "filter_type","gene_id"])


		####################################################################################
		# COMPARISON 1 Simple Filter
		####################################################################################
		
		#in-memory vs non-backed
		if i <= in_memory_high_filter:

			#memory AnnSQL approach
			if does_record_exist(i, "AnnSql In-Memory", "Filter1", gene_id) == False:
				mem_usage = memory_usage(annsql_filter_1)
				max_memory = (max(mem_usage) - min(mem_usage)) + load_annsql_mem_max_memory
				comparisons = add_to_comparisons(comparisons, i, "AnnSql In-Memory", max_memory, "Filter1", gene_id)
				print(f"memory annsql filter 1: {max_memory} MB")

			#memory AnnData approach
			if does_record_exist(i, "AnnData Not-Backed", "Filter1", gene_id) == False:
				mem_usage = memory_usage(anndata_filter_1)
				max_memory = (max(mem_usage) - min(mem_usage)) + load_anndata_mem_max_memory
				comparisons = add_to_comparisons(comparisons, i, "AnnData Not-Backed", max_memory, "Filter1", gene_id)
				print(f"memory annData filter 1: {max_memory} MB")
				


		#disk AnnSQL approach
		if does_record_exist(i, "AnnSql On-Disk", "Filter1", gene_id) == False:
			mem_usage = memory_usage(annsql_filter_1_disk)
			max_memory = (max(mem_usage) - min(mem_usage)) + load_annsql_disk_max_memory
			comparisons = add_to_comparisons(comparisons, i, "AnnSql On-Disk", max_memory, "Filter1", gene_id)
			print(f"disk annsql filter 1: {max_memory} MB")
			

		#disk AnnData approach
		if does_record_exist(i, "AnnData Backed", "Filter1", gene_id) == False:
			mem_usage = memory_usage(anndata_filter_1_disk)
			max_memory = (max(mem_usage) - min(mem_usage)) + load_annsql_disk_max_memory
			comparisons = add_to_comparisons(comparisons, i, "AnnData Backed", max_memory, "Filter1", gene_id)
			print(f"disk annData filter 1: {max_memory} MB")
			

		
		####################################################################################
		# COMPARISON 2 Complex Filter
		####################################################################################

		if i <= in_memory_high_filter:
			
			#memory AnnSQL approach
			if does_record_exist(i, "AnnSql In-Memory", "Filter2", gene_id) == False:
				mem_usage = memory_usage(annsql_filter_2)
				max_memory = (max(mem_usage) - min(mem_usage)) + load_annsql_mem_max_memory
				comparisons = add_to_comparisons(comparisons, i, "AnnSql In-Memory", max_memory, "Filter2", gene_id)
				print(f"memory annsql filter 2: {max_memory} MB")

			#memory AnnData approach
			if does_record_exist(i, "AnnData Not-Backed", "Filter2", gene_id) == False:
				mem_usage = memory_usage(anndata_filter_2)
				max_memory = (max(mem_usage) - min(mem_usage)) + load_anndata_mem_max_memory
				comparisons = add_to_comparisons(comparisons, i, "AnnData Not-Backed", max_memory, "Filter2", gene_id)
				print(f"memory annData filter 2: {max_memory} MB")


		#disk AnnSQL approach
		if does_record_exist(i, "AnnSql On-Disk", "Filter2", gene_id) == False:
			mem_usage = memory_usage(annsql_filter_2_disk)
			max_memory = (max(mem_usage) - min(mem_usage)) + load_annsql_disk_max_memory
			comparisons = add_to_comparisons(comparisons, i, "AnnSql On-Disk", max_memory, "Filter2", gene_id)
			print(f"disk annsql filter 2: {max_memory} MB")

		#disk AnnData approach
		if does_record_exist(i, "AnnData Backed", "Filter2", gene_id) == False:
			mem_usage = memory_usage(anndata_filter_2_disk)
			max_memory = (max(mem_usage) - min(mem_usage)) + load_annsql_disk_max_memory
			comparisons = add_to_comparisons(comparisons, i, "AnnData Backed", max_memory, "Filter2", gene_id)
			print(f"disk annData filter 2: {max_memory} MB")

		####################################################################################
		# COMPARISON 3 Complex Filter
		####################################################################################
		
		if i <= in_memory_high_filter:
			#memory AnnSQL approach
			if does_record_exist(i, "AnnSql In-Memory", "Filter3", gene_id) == False:
				mem_usage = memory_usage(annsql_filter_3)
				max_memory = (max(mem_usage) - min(mem_usage)) + load_annsql_mem_max_memory
				comparisons = add_to_comparisons(comparisons, i, "AnnSql In-Memory", max_memory, "Filter3", gene_id)
				print(f"memory annsql filter 3: {max_memory} MB")


			#memory AnnData approach
			if does_record_exist(i, "AnnData Not-Backed", "Filter3", gene_id) == False:
				mem_usage = memory_usage(anndata_filter_3)
				max_memory = (max(mem_usage) - min(mem_usage)) + load_anndata_mem_max_memory
				comparisons = add_to_comparisons(comparisons, i, "AnnData Not-Backed", max_memory, "Filter3", gene_id)
				print(f"memory annData filter 3: {max_memory} MB")
			
		#disk AnnSQL approach
		if does_record_exist(i, "AnnSql On-Disk", "Filter3", gene_id) == False:
			mem_usage = memory_usage(annsql_filter_3_disk)
			max_memory = (max(mem_usage) - min(mem_usage)) + load_annsql_disk_max_memory
			comparisons = add_to_comparisons(comparisons, i, "AnnSql On-Disk", max_memory, "Filter3", gene_id)
			print(f"disk annsql filter 3: {max_memory} MB")

		#disk AnnData approach
		if does_record_exist(i, "AnnData Backed", "Filter3", gene_id) == False:
			mem_usage = memory_usage(anndata_filter_3_disk)
			max_memory = (max(mem_usage) - min(mem_usage)) + load_annsql_disk_max_memory
			comparisons = add_to_comparisons(comparisons, i, "AnnData Backed", max_memory, "Filter3", gene_id)
			print(f"disk annData filter 3: {max_memory} MB")


		####################################################################################
		# COMPARISON 4 Complex Filter
		####################################################################################
		
		if i <= in_memory_high_filter:
			#memory AnnSQL approach
			if does_record_exist(i, "AnnSql In-Memory", "Filter4", gene_id) == False:
				mem_usage = memory_usage(annsql_filter_4)
				max_memory = (max(mem_usage) - min(mem_usage)) + load_annsql_mem_max_memory
				comparisons = add_to_comparisons(comparisons, i, "AnnSql In-Memory", max_memory, "Filter4", gene_id)
				print(f"memory annsql filter 4: {max_memory} MB")

			#memory AnnData approach
			if does_record_exist(i, "AnnData Not-Backed", "Filter4", gene_id) == False:
				mem_usage = memory_usage(anndata_filter_4)
				max_memory = (max(mem_usage) - min(mem_usage)) + load_anndata_mem_max_memory
				comparisons = add_to_comparisons(comparisons, i, "AnnData Not-Backed", max_memory, "Filter4", gene_id)
				print(f"memory annData filter 4: {max_memory} MB")


		#disk AnnSQL approach
		if does_record_exist(i, "AnnSql On-Disk", "Filter4", gene_id) == False:
			mem_usage = memory_usage(annsql_filter_4_disk)
			max_memory = (max(mem_usage) - min(mem_usage)) + load_annsql_disk_max_memory
			comparisons = add_to_comparisons(comparisons, i, "AnnSql On-Disk", max_memory, "Filter4", gene_id)
			print(f"disk annsql filter 4: {max_memory} MB")


		#disk AnnData approach
		if does_record_exist(i, "AnnData Backed", "Filter4", gene_id) == False:
			mem_usage = memory_usage(anndata_filter_4_disk)
			max_memory = (max(mem_usage) - min(mem_usage)) + load_annsql_disk_max_memory
			comparisons = add_to_comparisons(comparisons, i, "AnnData Backed", max_memory, "Filter4", gene_id)
			print(f"disk annData filter 4: {max_memory} MB")

		####################################################################################
		# COMPARISON 5 Complex Filter
		####################################################################################
		
		if i <= in_memory_high_filter:

			#memory AnnSQL approach
			if does_record_exist(i, "AnnSql In-Memory", "Filter5", gene_id) == False:
				mem_usage = memory_usage(annsql_filter_5)
				max_memory = (max(mem_usage) - min(mem_usage)) + load_annsql_mem_max_memory
				comparisons = add_to_comparisons(comparisons, i, "AnnSql In-Memory", max_memory, "Filter5", gene_id)
				print(f"memory annsql filter 5: {max_memory} MB")



			#memory AnnData approach
			if does_record_exist(i, "AnnData Not-Backed", "Filter5", gene_id) == False:
				mem_usage = memory_usage(anndata_filter_5)
				max_memory = (max(mem_usage) - min(mem_usage)) + load_anndata_mem_max_memory
				comparisons = add_to_comparisons(comparisons, i, "AnnData Not-Backed", max_memory, "Filter5", gene_id)
				print(f"memory annData filter 5: {max_memory} MB")


		#disk AnnSQL approach
		if does_record_exist(i, "AnnSql On-Disk", "Filter5", gene_id) == False:
			mem_usage = memory_usage(annsql_filter_5_disk)
			max_memory = (max(mem_usage) - min(mem_usage)) + load_annsql_disk_max_memory
			comparisons = add_to_comparisons(comparisons, i, "AnnSql On-Disk", max_memory, "Filter5", gene_id)
			print(f"disk annsql filter 5: {max_memory} MB")

		#disk AnnData approach
		if does_record_exist(i, "AnnData Backed", "Filter5", gene_id) == False:
			mem_usage = memory_usage(anndata_filter_5_disk)
			max_memory = (max(mem_usage) - min(mem_usage)) + load_annsql_disk_max_memory
			comparisons = add_to_comparisons(comparisons, i, "AnnData Backed", max_memory, "Filter5", gene_id)
			print(f"disk annData filter 5: {max_memory} MB")

		####################################################################################
		# COMPARISON 6 Complex Filter
		####################################################################################
		
		if i <= in_memory_high_filter:
			#memory AnnSQL approach
			if does_record_exist(i, "AnnSql In-Memory", "Filter6", gene_id) == False:
				mem_usage = memory_usage(annsql_filter_6)
				max_memory = (max(mem_usage) - min(mem_usage)) + load_annsql_mem_max_memory
				comparisons = add_to_comparisons(comparisons, i, "AnnSql In-Memory", max_memory, "Filter6", gene_id)
				print(f"memory annsql filter 6: {max_memory} MB")

			#memory AnnData approach
			if does_record_exist(i, "AnnData Not-Backed", "Filter6", gene_id) == False:
				mem_usage = memory_usage(anndata_filter_6)
				max_memory = (max(mem_usage) - min(mem_usage)) + load_anndata_mem_max_memory
				comparisons = add_to_comparisons(comparisons, i, "AnnData Not-Backed", max_memory, "Filter6", gene_id)
				print(f"memory annData filter 6: {max_memory} MB")

		#disk AnnSQL approach
		if does_record_exist(i, "AnnSql On-Disk", "Filter6", gene_id) == False:
			mem_usage = memory_usage(annsql_filter_6_disk)
			max_memory = (max(mem_usage) - min(mem_usage)) + load_annsql_disk_max_memory
			comparisons = add_to_comparisons(comparisons, i, "AnnSql On-Disk", max_memory, "Filter6", gene_id)
			print(f"disk annsql filter 6: {max_memory} MB")

		#disk AnnData approach
		if does_record_exist(i, "AnnData Backed", "Filter6", gene_id) == False:
			mem_usage = memory_usage(anndata_filter_6_disk)
			max_memory = (max(mem_usage) - min(mem_usage)) + load_annsql_disk_max_memory
			comparisons = add_to_comparisons(comparisons, i, "AnnData Backed", max_memory, "Filter6", gene_id)
			print(f"disk annData filter 6: {max_memory} MB")

		#open the comparison file 
		if os.path.exists(file_name):
			comparisons = pd.read_csv(file_name)

		#ntural log scale for memory
		comparisons['memory_log'] = np.log(comparisons['memory'])

		#log10 scale for memory
		comparisons['memory_log10'] = np.log10(comparisons['memory'])

		#add log scale for size
		comparisons['size_log'] = np.log(comparisons['size'])

		#add log10 scale for size
		comparisons['size_log10'] = np.log10(comparisons['size'])

		#write all of the log data to the file data
		comparisons.to_csv(file_name, index=False)


	#clean up memory
	adata_memory = None
	adata_sql_memory = None
	adata_disk = None
	adata_sql_disk = None
	gc.collect()
