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

data_path = "../data/splatter/"
db_path = "../db/splatter/"
in_memory_high_filter = 100001
comparison_records = []
dataset_sizes = [1000, 5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000, 75000, 100000, 250000]
file_name = "../results/comparisons_splatter.csv"

#set seed for reproducibility
np.random.seed(0)

#generate a list of 10 random numbers between 0 and 9999
gene_ids = np.random.randint(0, 9999, 20)
gene_ids_2 = np.random.randint(0, 9999, 20)

#write the gene ids to a file
with open("../results/gene_ids.txt", "w") as f:
	f.write("Set 1\n")
	for gene_id in gene_ids:
		f.write(str(gene_id) + "\n")
	f.write("\nSet 2\n")
	for gene_id in gene_ids_2:
		f.write(str(gene_id) + "\n")

#print the random nums
print(gene_ids)
print(gene_ids_2)

#open the comparison file does not exist
if not os.path.exists(file_name):
	comparisons = pd.DataFrame(columns=["size", "method", "runtime", "filter", "gene_id"])
	comparisons.to_csv(file_name, index=False)
else:
	comparisons = pd.DataFrame(columns=["size", "method", "runtime", "filter", "gene_id"])

#function to add to the comparisons dataframe
def add_to_comparisons(comparisons, size, method, runtime, filter, gene_id):
	comparisons = pd.concat([
		comparisons, 
		pd.DataFrame([{"size": size, "method": method, "runtime": runtime, "filter": filter, "gene_id": gene_id}])
	], ignore_index=False)
	comparisons.to_csv(file_name, index=False, mode='a', header=False)
	comparisons = pd.DataFrame(columns=["size", "method", "runtime", "filter", "gene_id"])

def does_record_exist(size, method, filter, gene_id):
	record_check = pd.read_csv(file_name)
	if len(record_check[(record_check['size'] == size) & (record_check['method'] == method) & (record_check['filter'] == filter) & (record_check['gene_id'] == gene_id)]) > 0:
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

	for gene_id in gene_ids:

		print(f"Size {i}:gene_id {gene_id}")

		####################################################################################
		# COMPARISON 1 Simple Filter
		####################################################################################
		
		#in-memory vs non-backed
		if i <= in_memory_high_filter:

			#memory AnnSQL approach
			if does_record_exist(i, "AnnSql In-Memory", "Filter1", gene_id) == False:
				start_time = time.time()
				adata_sql_memory.query(f"SELECT gene_{gene_id} FROM X WHERE gene_{gene_id} > 0.5 LIMIT 5")
				time1 = time.time() - start_time
				comparisons = add_to_comparisons(comparisons, i, "AnnSql In-Memory", time1, "Filter1", gene_id)
				print(f"memory annsql filter 1: {time1} seconds")

			#memory AnnData approach
			if does_record_exist(i, "AnnData Not-Backed", "Filter1", gene_id) == False:
				start_time = time.time()
				pd.DataFrame(adata_memory[adata_memory[:, "gene_"+str(gene_id)].X > 0.5, "gene_1"].X[:5])
				time2 = time.time() - start_time
				comparisons = add_to_comparisons(comparisons, i, "AnnData Not-Backed", time2, "Filter1", gene_id)
				print(f"memory annData filter 1: {time2} seconds")


		#disk AnnSQL approach
		if does_record_exist(i, "AnnSql On-Disk", "Filter1", gene_id) == False:
			start_time = time.time()
			adata_sql_disk.query(f"SELECT gene_{gene_id} FROM X WHERE gene_{gene_id} > 0.5 LIMIT 5")
			time3 = time.time() - start_time
			comparisons = add_to_comparisons(comparisons, i, "AnnSql On-Disk", time3, "Filter1", gene_id)
			print(f"disk annsql filter 1: {time3} seconds")

		#disk AnnData approach
		if does_record_exist(i, "AnnData Backed", "Filter1", gene_id) == False:
			start_time = time.time()
			pd.DataFrame(adata_disk[adata_disk[:, "gene_"+str(gene_id)].X > 0.5, "gene_1"].X[:5])
			time4 = time.time() - start_time
			comparisons = add_to_comparisons(comparisons, i, "AnnData Backed", time4, "Filter1", gene_id)
			print(f"disk annData filter 1: {time4} seconds")

		
		####################################################################################
		# COMPARISON 2 Complex Filter
		####################################################################################

		#select random pick from gene_ids_2
		gene_id_2 = gene_ids_2[np.random.randint(0, 20)]

		if i <= in_memory_high_filter:

			#memory AnnSQL approach
			if does_record_exist(i, "AnnSql In-Memory", "Filter2", gene_id) == False:
				start_time = time.time()
				adata_sql_memory.query(f"SELECT gene_{gene_id} FROM X WHERE gene_{gene_id} > 0.5 AND gene_{gene_id_2} < 0.5 LIMIT 5")
				time1 = time.time() - start_time
				comparisons = add_to_comparisons(comparisons, i, "AnnSql In-Memory", time1, "Filter2", gene_id)
				print(f"memory annsql filter 2: {time1} seconds")

			#memory AnnData approach
			if does_record_exist(i, "AnnData Not-Backed", "Filter2", gene_id) == False:
				start_time = time.time()
				pd.DataFrame(adata_memory[(adata_memory[:, "gene_"+str(gene_id)].X > 0.5) & (adata_memory[:, "gene_"+str(gene_id)].X < 0.5), "gene_1"].X[:5])
				time2 = time.time() - start_time
				comparisons = add_to_comparisons(comparisons, i, "AnnData Not-Backed", time2, "Filter2", gene_id)
				print(f"memory annData filter 2: {time2} seconds")


		#disk AnnSQL approach
		if does_record_exist(i, "AnnSql On-Disk", "Filter2", gene_id) == False:
			start_time = time.time()
			adata_sql_disk.query(f"SELECT gene_{gene_id} FROM X WHERE gene_{gene_id} > 0.5 AND gene_{gene_id_2} < 0.5 LIMIT 5")
			time3 = time.time() - start_time
			comparisons = add_to_comparisons(comparisons, i, "AnnSql On-Disk", time3, "Filter2", gene_id)
			print(f"disk annsql filter 2: {time3} seconds")

		#disk AnnData approach
		if does_record_exist(i, "AnnData Backed", "Filter2", gene_id) == False:
			start_time = time.time()
			pd.DataFrame(adata_disk[(adata_disk[:, "gene_"+str(gene_id)].X > 0.5) & (adata_disk[:, "gene_"+str(gene_id_2)].X < 0.5), "gene_1"].X[:5])
			time4 = time.time() - start_time
			comparisons = add_to_comparisons(comparisons, i, "AnnData Backed", time4, "Filter2", gene_id)
			print(f"disk annData filter 2: {time4} seconds")

		####################################################################################
		# COMPARISON 3 Complex Filter
		####################################################################################
		
		#select random picks from gene_ids_2
		gene_id_2 = gene_ids_2[np.random.randint(0, 20)]
		gene_id_3 = gene_ids_2[np.random.randint(0, 20)]

		if i <= in_memory_high_filter:
			#memory AnnSQL approach
			if does_record_exist(i, "AnnSql In-Memory", "Filter3", gene_id) == False:
				start_time = time.time()
				adata_sql_memory.query(f"SELECT gene_{gene_id} FROM X WHERE gene_{gene_id} > 0.5 AND gene_{gene_id_2} < 0.5 AND gene_{gene_id_3} > 0.5 LIMIT 5")
				time1 = time.time() - start_time
				comparisons = add_to_comparisons(comparisons, i, "AnnSql In-Memory", time1, "Filter3", gene_id)
				print(f"memory annsql filter 3: {time1} seconds")

			#memory AnnData approach
			if does_record_exist(i, "AnnData Not-Backed", "Filter3", gene_id) == False:
				start_time = time.time()
				pd.DataFrame(adata_memory[(adata_memory[:, "gene_"+str(gene_id)].X > 0.5) & (adata_memory[:, "gene_"+str(gene_id_2)].X < 0.5) & (adata_memory[:, "gene_"+str(gene_id_3)].X > 0.5), "gene_"+str(gene_id)].X[:5])
				time2 = time.time() - start_time
				comparisons = add_to_comparisons(comparisons, i, "AnnData Not-Backed", time2, "Filter3", gene_id)
				print(f"memory annData filter 3: {time2} seconds")
			
		#disk AnnSQL approach
		if does_record_exist(i, "AnnSql On-Disk", "Filter3", gene_id) == False:
			start_time = time.time()
			adata_sql_disk.query(f"SELECT gene_{gene_id} FROM X WHERE gene_{gene_id} > 0.5 AND gene_{gene_id_2} < 0.5 AND gene_{gene_id_3} > 0.5 LIMIT 5")
			time3 = time.time() - start_time
			comparisons = add_to_comparisons(comparisons, i, "AnnSql On-Disk", time3, "Filter3", gene_id)
			print(f"disk annsql filter 3: {time3} seconds")

		#disk AnnData approach
		if does_record_exist(i, "AnnData Backed", "Filter3", gene_id) == False:
			start_time = time.time()
			pd.DataFrame(adata_disk[(adata_disk[:, "gene_"+str(gene_id)].X > 0.5) & (adata_disk[:, "gene_"+str(gene_id_2)].X < 0.5) & (adata_disk[:, "gene_"+str(gene_id_3)].X > 0.5), "gene_"+str(gene_id)].X[:5])
			time4 = time.time() - start_time
			comparisons = add_to_comparisons(comparisons, i, "AnnData Backed", time4, "Filter3", gene_id)
			print(f"disk annData filter 3: {time4} seconds")



		####################################################################################
		# COMPARISON 4 Complex Filter
		####################################################################################
		
		#select random picks from gene_ids_2
		gene_id_2 = gene_ids_2[np.random.randint(0, 20)]
		gene_id_3 = gene_ids_2[np.random.randint(0, 20)]
		gene_id_4 = gene_ids_2[np.random.randint(0, 20)]

		if i <= in_memory_high_filter:
			#memory AnnSQL approach
			if does_record_exist(i, "AnnSql In-Memory", "Filter4", gene_id) == False:
				start_time = time.time()
				adata_sql_memory.query(f"SELECT gene_{gene_id} FROM X WHERE gene_{gene_id} > 0.5 AND gene_{gene_id_2} < 0.5 AND gene_{gene_id_3} > 0.5 AND gene_{gene_id_4} < 0.5 LIMIT 5")
				time1 = time.time() - start_time
				comparisons = add_to_comparisons(comparisons, i, "AnnSql In-Memory", time1, "Filter4", gene_id)
				print(f"memory annsql filter 4: {time1} seconds")

			#memory AnnData approach
			if does_record_exist(i, "AnnData Not-Backed", "Filter4", gene_id) == False:
				start_time = time.time()
				pd.DataFrame(adata_memory[(adata_memory[:, f"gene_{gene_id}"].X > 0.5) & (adata_memory[:, f"gene_{gene_id_2}"].X < 0.5) & (adata_memory[:, f"gene_{gene_id_3}"].X > 0.5) & (adata_memory[:, f"gene_{gene_id_4}"].X < 0.5), f"gene_{gene_id}"].X[:5])
				time2 = time.time() - start_time
				comparisons = add_to_comparisons(comparisons, i, "AnnData Not-Backed", time2, "Filter4", gene_id)
				print(f"memory annData filter 4: {time2} seconds")


		#disk AnnSQL approach
		if does_record_exist(i, "AnnSql On-Disk", "Filter4", gene_id) == False:
			start_time = time.time()
			adata_sql_disk.query(f"SELECT gene_{gene_id} FROM X WHERE gene_{gene_id} > 0.5 AND gene_{gene_id_2} < 0.5 AND gene_{gene_id_3} > 0.5 AND gene_{gene_id_4} < 0.5 LIMIT 5")
			time3 = time.time() - start_time
			comparisons = add_to_comparisons(comparisons, i, "AnnSql On-Disk", time3, "Filter4", gene_id)
			print(f"disk annsql filter 4: {time3} seconds")


		#disk AnnData approach
		if does_record_exist(i, "AnnData Backed", "Filter4", gene_id) == False:
			start_time = time.time()
			pd.DataFrame(adata_disk[(adata_disk[:, f"gene_{gene_id}"].X > 0.5) & (adata_disk[:, f"gene_{gene_id_2}"].X < 0.5) & (adata_disk[:, f"gene_{gene_id_3}"].X > 0.5) & (adata_disk[:, f"gene_{gene_id_4}"].X < 0.5), f"gene_{gene_id}"].X[:5])
			time4 = time.time() - start_time
			comparisons = add_to_comparisons(comparisons, i, "AnnData Backed", time4, "Filter4", gene_id)
			print(f"disk annData filter 4: {time4} seconds")


		####################################################################################
		# COMPARISON 5 Complex Filter
		####################################################################################
		
		#select random picks from gene_ids_2
		gene_id_2 = gene_ids_2[np.random.randint(0, 20)]
		gene_id_3 = gene_ids_2[np.random.randint(0, 20)]
		gene_id_4 = gene_ids_2[np.random.randint(0, 20)]

		if i <= in_memory_high_filter:

			#memory AnnSQL approach
			if does_record_exist(i, "AnnSql In-Memory", "Filter5", gene_id) == False:
				start_time = time.time()
				adata_sql_memory.query(f"""SELECT 
											gene_{gene_id}, gene_{gene_id_2}, gene_{gene_id_3}, gene_{gene_id_4}, (gene_{gene_id_3}/gene_{gene_id_4}) as gene_{gene_id_3}_{gene_id_4}
										FROM 
											X 
										WHERE 
											(gene_{gene_id} > 0.5 AND gene_{gene_id_2} < 0.3) 
										AND 
											(gene_{gene_id_3}/gene_{gene_id_4} > 0 AND gene_{gene_id_3}/gene_{gene_id_4} < 0.5) 
										LIMIT 5""")
				time1 = time.time() - start_time
				comparisons = add_to_comparisons(comparisons, i, "AnnSql In-Memory", time1, "Filter5", gene_id)
				print(f"memory annsql filter 5: {time1} seconds")

			#memory AnnData approach
			if does_record_exist(i, "AnnData Not-Backed", "Filter5", gene_id) == False:
				start_time = time.time()
				gene_1_vals = adata_memory[:, f"gene_{gene_id}"].X.flatten()
				gene_2_vals = adata_memory[:, f"gene_{gene_id_2}"].X.flatten()
				gene_3_vals = adata_memory[:, f"gene_{gene_id_3}"].X.flatten()
				gene_4_vals = adata_memory[:, f"gene_{gene_id_4}"].X.flatten()
				d_3_4 = gene_3_vals / gene_4_vals
				condition = (gene_1_vals > 0.5) & (gene_2_vals < 0.3) & (d_3_4 > 0) & (d_3_4 < 0.5)
				filtered_indices = np.where(condition)[0]
				limited_indices = filtered_indices[:5]
				pd.DataFrame({
					f"gene_{gene_id}": gene_1_vals[limited_indices],
					f"gene_{gene_id_2}": gene_2_vals[limited_indices],
					f"gene_{gene_id_3}": gene_3_vals[limited_indices],
					f"gene_{gene_id_4}": gene_4_vals[limited_indices],
					f"gene_{gene_id_3}_{gene_id_4}": d_3_4[limited_indices]
				})
				time2 = time.time() - start_time
				comparisons = add_to_comparisons(comparisons, i, "AnnData Not-Backed", time2, "Filter5", gene_id)
				print(f"memory annData filter 5: {time2} seconds")

		#disk AnnSQL approach
		if does_record_exist(i, "AnnSql On-Disk", "Filter5", gene_id) == False:
			start_time = time.time()
			adata_sql_disk.query(f"""SELECT 
										gene_{gene_id}, gene_{gene_id_2}, gene_{gene_id_3}, gene_{gene_id_4}, (gene_{gene_id_3}/gene_{gene_id_4}) as gene_{gene_id_3}_{gene_id_4}
									FROM 
										X 
									WHERE 
										(gene_{gene_id} > 0.5 AND gene_{gene_id_2} < 0.3) 
									AND 
										(gene_{gene_id_3}/gene_{gene_id_4} > 0 AND gene_{gene_id_3}/gene_{gene_id_4} < 0.5) 
									LIMIT 5""")
			time3 = time.time() - start_time
			comparisons = add_to_comparisons(comparisons, i, "AnnSql On-Disk", time3, "Filter5", gene_id)
			print(f"disk annsql filter 5: {time3} seconds")

		#disk AnnData approach
		if does_record_exist(i, "AnnData Backed", "Filter5", gene_id) == False:
			start_time = time.time()
			gene_1_vals = adata_disk[:, f"gene_{gene_id}"].X.flatten()
			gene_2_vals = adata_disk[:, f"gene_{gene_id_2}"].X.flatten()
			gene_3_vals = adata_disk[:, f"gene_{gene_id_3}"].X.flatten()
			gene_4_vals = adata_disk[:, f"gene_{gene_id_4}"].X.flatten()
			d_3_4 = gene_3_vals / gene_4_vals
			condition = (gene_1_vals > 0.5) & (gene_2_vals < 0.3) & (d_3_4 > 0) & (d_3_4 < 0.5)
			filtered_indices = np.where(condition)[0]
			limited_indices = filtered_indices[:5]
			pd.DataFrame({
				f"gene_{gene_id}": gene_1_vals[limited_indices],
				f"gene_{gene_id_2}": gene_2_vals[limited_indices],
				f"gene_{gene_id_3}": gene_3_vals[limited_indices],
				f"gene_{gene_id_4}": gene_4_vals[limited_indices],
				f"gene_{gene_id_3}_{gene_id_4}": d_3_4[limited_indices]
			})
			end_time = time.time()
			time4 = end_time - start_time
			comparisons = add_to_comparisons(comparisons, i, "AnnData Backed", time4, "Filter5", gene_id)
			print(f"disk annData filter 5: {time4} seconds")


		####################################################################################
		# COMPARISON 6 Complex Filter
		####################################################################################
		
		if i <= in_memory_high_filter:

			#memory AnnSQL approach
			if does_record_exist(i, "AnnSql In-Memory", "Filter6", gene_id) == False:
				start_time = time.time()
				adata_sql_memory.query(f"SELECT obs.cell_type, AVG(gene_{gene_id}) as avg_gene_{gene_id} FROM X INNER JOIN obs ON X.cell_id = obs.cell_id GROUP BY obs.cell_type")
				time1 = time.time() - start_time
				comparisons = add_to_comparisons(comparisons, i, "AnnSql In-Memory", time1, "Filter6", gene_id)
				print(f"memory annsql filter 6: {time1} seconds")

			#memory AnnData approach
			if does_record_exist(i, "AnnData Not-Backed", "Filter6", gene_id) == False:
				start_time = time.time()
				pd.DataFrame({
					'cell_type': adata_memory.obs['cell_type'],
					f'gene_{gene_id}': adata_memory[:, f'gene_{gene_id}'].X.flatten()
				}).groupby('cell_type')[f'gene_{gene_id}'].mean().reset_index()
				end_time = time.time()
				time2 = end_time - start_time
				comparisons = add_to_comparisons(comparisons, i, "AnnData Not-Backed", time2, "Filter6", gene_id)
				print(f"memory annData filter 6: {time2} seconds")

		#disk AnnSQL approach
		if does_record_exist(i, "AnnSql On-Disk", "Filter6", gene_id) == False:
			start_time = time.time()
			adata_sql_disk.query(f"SELECT obs.cell_type, AVG(gene_{gene_id}) as avg_gene_{gene_id} FROM X INNER JOIN obs ON X.cell_id = obs.cell_id GROUP BY obs.cell_type")
			end_time = time.time()
			time3 = end_time - start_time
			comparisons = add_to_comparisons(comparisons, i, "AnnSql On-Disk", time3, "Filter6", gene_id)
			print(f"disk annsql filter 6: {time3} seconds")

		#disk AnnData approach
		if does_record_exist(i, "AnnData Backed", "Filter6", gene_id) == False:
			start_time = time.time()
			pd.DataFrame({
				'cell_type': adata_disk.obs['cell_type'],
				f'gene_{gene_id}': adata_disk[:, f'gene_{gene_id}'].X.flatten()
			}).groupby('cell_type')[f'gene_{gene_id}'].mean().reset_index()
			end_time = time.time()
			time4 = end_time - start_time
			comparisons = add_to_comparisons(comparisons, i, "AnnData Backed", time4, "Filter6", gene_id)
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