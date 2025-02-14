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

data_path = "../data/random/"
db_path = "../db/random/"
file_name = "../results/comparisons_random_7.csv"
in_memory_high_filter = 100001
dataset_sizes = [1000, 5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000, 75000, 100000, 250000]

#open the comparison file does not exist
if not os.path.exists(file_name):
	comparisons = pd.DataFrame(columns=["size", "method", "runtime", "filter_type"])
	comparisons.to_csv(file_name, index=False)
else:
	comparisons = pd.DataFrame(columns=["size", "method", "runtime", "filter_type"])

#function to add to the comparisons dataframe
def add_to_comparisons(comparisons, size, method, runtime, filter_type):
	comparisons = pd.concat([
		comparisons, 
		pd.DataFrame([{"size": size, "method": method, "runtime": runtime, "filter_type": filter_type}])
	], ignore_index=False)
	comparisons.to_csv(file_name, index=False, mode='a', header=False)
	comparisons = pd.DataFrame(columns=["size", "method", "runtime", "filter_type"])

def does_record_exist(size, method, filter_type):
	record_check = pd.read_csv(file_name)
	if len(record_check[(record_check['size'] == size) & (record_check['method'] == method) & (record_check['filter_type'] == filter_type)]) > 0:
		return True
	return False


for i in dataset_sizes:
	if not os.path.exists(data_path+"data_" + str(i) + ".h5ad"):
		continue

	print(f"Running for {i}")

	#in-memory vs non-backed
	if i <= in_memory_high_filter:
		adata_memory = sc.read(data_path+"data_" + str(i) + ".h5ad")
		adata_sql_memory = AnnSQL(adata=data_path+"data_" + str(i) + ".h5ad", print_output=False)

	#on-disk vs backed
	adata_disk = sc.read(data_path+"data_" + str(i) + ".h5ad", backed="r")
	adata_sql_disk = AnnSQL(db=db_path+"data_" + str(i) + ".asql", print_output=False)

	####################################################################################
	# COMPARISON 1 Simple Filter
	####################################################################################
	
	#in-memory vs non-backed
	if i <= in_memory_high_filter:

		#memory AnnSQL approach
		if does_record_exist(i, "AnnSql In-Memory", "Filter1") == False:
			start_time = time.time()
			adata_sql_memory.query("SELECT gene_1 FROM X WHERE gene_1 > 0.5 LIMIT 5")
			time1 = time.time() - start_time
			comparisons = add_to_comparisons(comparisons, i, "AnnSql In-Memory", time1, "Filter1")
			print(f"memory annsql filter 1: {time1} seconds")

		#memory AnnData approach
		if does_record_exist(i, "AnnData Not-Backed", "Filter1") == False:
			start_time = time.time()
			pd.DataFrame(adata_memory[adata_memory[:, "gene_1"].X > 0.5, "gene_1"].X[:5])
			time2 = time.time() - start_time
			comparisons = add_to_comparisons(comparisons, i, "AnnData Not-Backed", time2, "Filter1")
			print(f"memory annData filter 1: {time2} seconds")


	#disk AnnSQL approach
	if does_record_exist(i, "AnnSql On-Disk", "Filter1") == False:
		start_time = time.time()
		adata_sql_disk.query("SELECT gene_1 FROM X WHERE gene_1 > 0.5 LIMIT 5")
		time3 = time.time() - start_time
		comparisons = add_to_comparisons(comparisons, i, "AnnSql On-Disk", time3, "Filter1")
		print(f"disk annsql filter 1: {time3} seconds")

	#disk AnnData approach
	if does_record_exist(i, "AnnData Backed", "Filter1") == False:
		start_time = time.time()
		pd.DataFrame(adata_disk[adata_disk[:, "gene_1"].X > 0.5, "gene_1"].X[:5])
		time4 = time.time() - start_time
		comparisons = add_to_comparisons(comparisons, i, "AnnData Backed", time4, "Filter1")
		print(f"disk annData filter 1: {time4} seconds")

	
	####################################################################################
	# COMPARISON 2 Complex Filter
	####################################################################################

	if i <= in_memory_high_filter:
		
		#memory AnnSQL approach
		if does_record_exist(i, "AnnSql In-Memory", "Filter2") == False:
			start_time = time.time()
			adata_sql_memory.query("SELECT gene_1 FROM X WHERE gene_1 > 0.5 AND gene_2 < 0.5 LIMIT 5")
			time1 = time.time() - start_time
			comparisons = add_to_comparisons(comparisons, i, "AnnSql In-Memory", time1, "Filter2")
			print(f"memory annsql filter 2: {time1} seconds")

		#memory AnnData approach
		if does_record_exist(i, "AnnData Not-Backed", "Filter2") == False:
			start_time = time.time()
			pd.DataFrame(adata_memory[(adata_memory[:, "gene_1"].X > 0.5) & (adata_memory[:, "gene_2"].X < 0.5), "gene_1"].X[:5])
			time2 = time.time() - start_time
			comparisons = add_to_comparisons(comparisons, i, "AnnData Not-Backed", time2, "Filter2")
			print(f"memory annData filter 2: {time2} seconds")


	#disk AnnSQL approach
	if does_record_exist(i, "AnnSql On-Disk", "Filter2") == False:
		start_time = time.time()
		adata_sql_disk.query("SELECT gene_1 FROM X WHERE gene_1 > 0.5 AND gene_2 < 0.5 LIMIT 5")
		time3 = time.time() - start_time
		comparisons = add_to_comparisons(comparisons, i, "AnnSql On-Disk", time3, "Filter2")
		print(f"disk annsql filter 2: {time3} seconds")

	#disk AnnData approach
	if does_record_exist(i, "AnnData Backed", "Filter2") == False:
		start_time = time.time()
		pd.DataFrame(adata_disk[(adata_disk[:, "gene_1"].X > 0.5) & (adata_disk[:, "gene_2"].X < 0.5), "gene_1"].X[:5])
		time4 = time.time() - start_time
		comparisons = add_to_comparisons(comparisons, i, "AnnData Backed", time4, "Filter2")
		print(f"disk annData filter 2: {time4} seconds")

	####################################################################################
	# COMPARISON 3 Complex Filter
	####################################################################################
	
	if i <= in_memory_high_filter:
		#memory AnnSQL approach
		if does_record_exist(i, "AnnSql In-Memory", "Filter3") == False:
			start_time = time.time()
			adata_sql_memory.query("SELECT gene_1 FROM X WHERE gene_1 > 0.5 AND gene_2 < 0.5 AND gene_3 > 0.5 LIMIT 5")
			time1 = time.time() - start_time
			comparisons = add_to_comparisons(comparisons, i, "AnnSql In-Memory", time1, "Filter3")
			print(f"memory annsql filter 3: {time1} seconds")

		#memory AnnData approach
		if does_record_exist(i, "AnnData Not-Backed", "Filter3") == False:
			start_time = time.time()
			pd.DataFrame(adata_memory[(adata_memory[:, "gene_1"].X > 0.5) & (adata_memory[:, "gene_2"].X < 0.5) & (adata_memory[:, "gene_3"].X > 0.5), "gene_1"].X[:5])
			time2 = time.time() - start_time
			comparisons = add_to_comparisons(comparisons, i, "AnnData Not-Backed", time2, "Filter3")
			print(f"memory annData filter 3: {time2} seconds")
		
	#disk AnnSQL approach
	if does_record_exist(i, "AnnSql On-Disk", "Filter3") == False:
		start_time = time.time()
		adata_sql_disk.query("SELECT gene_1 FROM X WHERE gene_1 > 0.5 AND gene_2 < 0.5 AND gene_3 > 0.5 LIMIT 5")
		time3 = time.time() - start_time
		comparisons = add_to_comparisons(comparisons, i, "AnnSql On-Disk", time3, "Filter3")
		print(f"disk annsql filter 3: {time3} seconds")

	#disk AnnData approach
	if does_record_exist(i, "AnnData Backed", "Filter3") == False:
		start_time = time.time()
		pd.DataFrame(adata_disk[(adata_disk[:, "gene_1"].X > 0.5) & (adata_disk[:, "gene_2"].X < 0.5) & (adata_disk[:, "gene_3"].X > 0.5), "gene_1"].X[:5])
		time4 = time.time() - start_time
		comparisons = add_to_comparisons(comparisons, i, "AnnData Backed", time4, "Filter3")
		print(f"disk annData filter 3: {time4} seconds")



	####################################################################################
	# COMPARISON 4 Complex Filter
	####################################################################################
	
	if i <= in_memory_high_filter:
		#memory AnnSQL approach
		if does_record_exist(i, "AnnSql In-Memory", "Filter4") == False:
			start_time = time.time()
			adata_sql_memory.query("SELECT gene_1 FROM X WHERE gene_1 > 0.5 AND gene_2 < 0.5 AND gene_3 > 0.5 AND gene_4 < 0.5 LIMIT 5")
			time1 = time.time() - start_time
			comparisons = add_to_comparisons(comparisons, i, "AnnSql In-Memory", time1, "Filter4")
			print(f"memory annsql filter 4: {time1} seconds")

		#memory AnnData approach
		if does_record_exist(i, "AnnData Not-Backed", "Filter4") == False:
			start_time = time.time()
			pd.DataFrame(adata_memory[(adata_memory[:, "gene_1"].X > 0.5) & (adata_memory[:, "gene_2"].X < 0.5) & (adata_memory[:, "gene_3"].X > 0.5) & (adata_memory[:, "gene_4"].X < 0.5), "gene_1"].X[:5])
			time2 = time.time() - start_time
			comparisons = add_to_comparisons(comparisons, i, "AnnData Not-Backed", time2, "Filter4")
			print(f"memory annData filter 4: {time2} seconds")


	#disk AnnSQL approach
	if does_record_exist(i, "AnnSql On-Disk", "Filter4") == False:
		start_time = time.time()
		adata_sql_disk.query("SELECT gene_1 FROM X WHERE gene_1 > 0.5 AND gene_2 < 0.5 AND gene_3 > 0.5 AND gene_4 < 0.5 LIMIT 5")
		time3 = time.time() - start_time
		comparisons = add_to_comparisons(comparisons, i, "AnnSql On-Disk", time3, "Filter4")
		print(f"disk annsql filter 4: {time3} seconds")


	#disk AnnData approach
	if does_record_exist(i, "AnnData Backed", "Filter4") == False:
		start_time = time.time()
		pd.DataFrame(adata_disk[(adata_disk[:, "gene_1"].X > 0.5) & (adata_disk[:, "gene_2"].X < 0.5) & (adata_disk[:, "gene_3"].X > 0.5) & (adata_disk[:, "gene_4"].X < 0.5), "gene_1"].X[:5])
		time4 = time.time() - start_time
		comparisons = add_to_comparisons(comparisons, i, "AnnData Backed", time4, "Filter4")
		print(f"disk annData filter 4: {time4} seconds")


	####################################################################################
	# COMPARISON 5 Complex Filter
	####################################################################################
	
	if i <= in_memory_high_filter:

		#memory AnnSQL approach
		if does_record_exist(i, "AnnSql In-Memory", "Filter5") == False:
			start_time = time.time()
			adata_sql_memory.query("""SELECT 
										gene_1, gene_2, gene_3, gene_4, (gene_3/gene_4) as gene_3_4
									FROM 
										X 
									WHERE 
										(gene_1 > 0.5 AND gene_2 < 0.3) 
									AND 
										(gene_3/gene_4 > 0 AND gene_3/gene_4 < 0.5) 
									LIMIT 5""")
			time1 = time.time() - start_time
			comparisons = add_to_comparisons(comparisons, i, "AnnSql In-Memory", time1, "Filter5")
			print(f"memory annsql filter 5: {time1} seconds")



		#memory AnnData approach
		if does_record_exist(i, "AnnData Not-Backed", "Filter5") == False:
			start_time = time.time()
			gene_1_vals = adata_memory[:, "gene_1"].X.flatten()
			gene_2_vals = adata_memory[:, "gene_2"].X.flatten()
			gene_3_vals = adata_memory[:, "gene_3"].X.flatten()
			gene_4_vals = adata_memory[:, "gene_4"].X.flatten()
			d_3_4 = gene_3_vals / gene_4_vals
			condition = (gene_1_vals > 0.5) & (gene_2_vals < 0.3) & (d_3_4 > 0) & (d_3_4 < 0.5)
			filtered_indices = np.where(condition)[0]
			limited_indices = filtered_indices[:5]
			pd.DataFrame({
				"gene_1": gene_1_vals[limited_indices],
				"gene_2": gene_2_vals[limited_indices],
				"gene_3": gene_3_vals[limited_indices],
				"gene_4": gene_4_vals[limited_indices],
				"gene_3_4": d_3_4[limited_indices]
			})
			time2 = time.time() - start_time
			comparisons = add_to_comparisons(comparisons, i, "AnnData Not-Backed", time2, "Filter5")
			print(f"memory annData filter 5: {time2} seconds")


	#disk AnnSQL approach
	if does_record_exist(i, "AnnSql On-Disk", "Filter5") == False:
		start_time = time.time()
		adata_sql_disk.query("""SELECT 
								gene_1, gene_2, gene_3, gene_4, (gene_3/gene_4) as gene_3_4
							FROM 
								X 
							WHERE 
								(gene_1 > 0.5 AND gene_2 < 0.3) 
							AND 
								(gene_3/gene_4 > 0 AND gene_3/gene_4 < 0.5) 
							LIMIT 5""")
		time3 = time.time() - start_time
		comparisons = add_to_comparisons(comparisons, i, "AnnSql On-Disk", time3, "Filter5")
		print(f"disk annsql filter 5: {time3} seconds")

	#disk AnnData approach
	if does_record_exist(i, "AnnData Backed", "Filter5") == False:
		start_time = time.time()
		gene_1_vals = adata_disk[:, "gene_1"].X.flatten()
		gene_2_vals = adata_disk[:, "gene_2"].X.flatten()
		gene_3_vals = adata_disk[:, "gene_3"].X.flatten()
		gene_4_vals = adata_disk[:, "gene_4"].X.flatten()
		d_3_4 = gene_3_vals / gene_4_vals
		condition = (gene_1_vals > 0.5) & (gene_2_vals < 0.3) & (d_3_4 > 0) & (d_3_4 < 0.5)
		filtered_indices = np.where(condition)[0]
		limited_indices = filtered_indices[:5]
		pd.DataFrame({
			"gene_1": gene_1_vals[limited_indices],
			"gene_2": gene_2_vals[limited_indices],
			"gene_3": gene_3_vals[limited_indices],
			"gene_4": gene_4_vals[limited_indices],
			"gene_3_4": d_3_4[limited_indices]
		})
		end_time = time.time()
		time4 = end_time-start_time
		comparisons = add_to_comparisons(comparisons, i, "AnnData Backed", time4, "Filter5")
		print(f"disk annData filter 5: {time4} seconds")

	####################################################################################
	# COMPARISON 6 Complex Filter
	####################################################################################
	
	if i <= in_memory_high_filter:
		#memory AnnSQL approach
		if does_record_exist(i, "AnnSql In-Memory", "Filter6") == False:
			start_time = time.time()
			adata_sql_memory.query("SELECT obs.cell_type, AVG(gene_1) as avg_gene_1 FROM X INNER JOIN obs ON X.cell_id = obs.cell_id GROUP BY obs.cell_type")
			time1 = time.time() - start_time
			comparisons = add_to_comparisons(comparisons, i, "AnnSql In-Memory", time1, "Filter6")
			print(f"memory annsql filter 6: {time1} seconds")

		#memory AnnData approach
		if does_record_exist(i, "AnnData Not-Backed", "Filter6") == False:
			start_time = time.time()
			pd.DataFrame({
				'cell_type': adata_memory.obs['cell_type'],
				'gene_1': adata_memory[:, 'gene_1'].X.flatten()
			}).groupby('cell_type')['gene_1'].mean().reset_index()
			end_time = time.time()
			time2 = end_time-start_time
			comparisons = add_to_comparisons(comparisons, i, "AnnData Not-Backed", time2, "Filter6")
			print(f"memory annData filter 6: {time2} seconds")

	#disk AnnSQL approach
	if does_record_exist(i, "AnnSql On-Disk", "Filter6") == False:
		start_time = time.time()
		adata_sql_disk.query("SELECT obs.cell_type, AVG(gene_1) as avg_gene_1 FROM X INNER JOIN obs ON X.cell_id = obs.cell_id GROUP BY obs.cell_type")
		end_time = time.time()
		time3 = end_time-start_time
		comparisons = add_to_comparisons(comparisons, i, "AnnSql On-Disk", time3, "Filter6")
		print(f"disk annsql filter 6: {time3} seconds")

	#disk AnnData approach
	if does_record_exist(i, "AnnData Backed", "Filter6") == False:
		start_time = time.time()
		pd.DataFrame({
			'cell_type': adata_disk.obs['cell_type'],
			'gene_1': adata_disk[:, 'gene_1'].X.flatten()
		}).groupby('cell_type')['gene_1'].mean().reset_index()
		end_time = time.time()
		time4 = end_time-start_time
		comparisons = add_to_comparisons(comparisons, i, "AnnData Backed", time4, "Filter6")
		print(f"disk annData filter 6: {time4} seconds")

	#clean up memory
	adata_memory = None
	adata_sql_memory = None
	adata_disk = None
	adata_sql_disk = None
	gc.collect()


#open the comparison file 
if os.path.exists(file_name):
	comparisons = pd.read_csv(file_name)

#ntural log scale for runtime
comparisons['runtime_log'] = np.log(comparisons['runtime'])

#log10 scale for runtime
comparisons['runtime_log10'] = np.log10(comparisons['runtime'])

#add log scale for size
comparisons['size_log'] = np.log(comparisons['size'])

#add log10 scale for size
comparisons['size_log10'] = np.log10(comparisons['size'])

#write all of the log data to the file data
comparisons.to_csv(file_name, index=False)