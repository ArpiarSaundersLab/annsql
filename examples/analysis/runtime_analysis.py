import scanpy as sc
from AnnSQL import AnnSQL
import time
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import gc

in_memory_high_filter = 100001
comparison_records = []
dataset_sizes = [1000, 5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000, 75000, 100000, 250000]

for i in dataset_sizes:
	if not os.path.exists("../data/random_data_" + str(i) + ".h5ad"):
		continue

	print(f"Running for {i}")

	#in-memory vs non-backed
	if i <= in_memory_high_filter:
		adata_memory = sc.read("../data/random_data_" + str(i) + ".h5ad")
		adata_sql_memory = AnnSQL(adata="../data/random_data_" + str(i) + ".h5ad")

	#on-disk vs backed
	adata_disk = sc.read("../data/random_data_" + str(i) + ".h5ad", backed="r")
	adata_sql_disk = AnnSQL(db="../db/random_data_" + str(i) + ".asql")

	####################################################################################
	# COMPARISON 1 Simple Filter
	####################################################################################
	
	#in-memory vs non-backed
	if i <= in_memory_high_filter:

		#memory AnnSQL approach
		start_time = time.time()
		adata_sql_memory.query("SELECT gene_1 FROM X WHERE gene_1 > 0.5 LIMIT 5")
		time1 = time.time() - start_time
		print(f"memory annsql filter 1: {time1} seconds")

		#memory AnnData approach
		# Why wrap in pd.DataFrame? AnnSQL returns a DataFrame, so we attempt to compare apples to apples
		start_time = time.time()
		pd.DataFrame(adata_memory[adata_memory[:, "gene_1"].X > 0.5, "gene_1"].X[:5], columns=["gene_1"])
		time2 = time.time() - start_time
		print(f"memory annData filter 1: {time2} seconds")

		#Store comparison results
		comparison_records.append([i, "AnnSql In-Memory", time1, "Filter1"])
		comparison_records.append([i, "AnnData Not-Backed", time2, "Filter1"])

	#disk AnnSQL approach
	start_time = time.time()
	adata_sql_disk.query("SELECT gene_1 FROM X WHERE gene_1 > 0.5 LIMIT 5")
	time3 = time.time() - start_time
	print(f"disk annsql filter 1: {time3} seconds")

	#disk AnnData approach
	start_time = time.time()
	pd.DataFrame(adata_disk[adata_disk[:, "gene_1"].X > 0.5, "gene_1"].X[:5], columns=["gene_1"])
	time4 = time.time() - start_time
	print(f"disk annData filter 1: {time4} seconds")

	#Store comparison results
	comparison_records.append([i, "AnnSql On-Disk", time3, "Filter1"])
	comparison_records.append([i, "AnnData Backed", time4, "Filter1"])

	####################################################################################
	# COMPARISON 2 Complex Filter
	####################################################################################

	if i <= in_memory_high_filter:
		
		#memory AnnSQL approach
		start_time = time.time()
		adata_sql_memory.query("SELECT gene_1 FROM X WHERE gene_1 > 0.5 AND gene_2 < 0.5 LIMIT 5")
		time1 = time.time() - start_time
		print(f"memory annsql filter 2: {time1} seconds")

		#memory AnnData approach
		start_time = time.time()
		pd.DataFrame(adata_memory[(adata_memory[:, "gene_1"].X > 0.5) & (adata_memory[:, "gene_2"].X < 0.5), "gene_1"].X[:5], columns=["gene_1"])
		time2 = time.time() - start_time
		print(f"memory annData filter 2: {time2} seconds")

		#Store comparison results
		comparison_records.append([i, "AnnSql In-Memory", time1, "Filter2"])
		comparison_records.append([i, "AnnData Not-Backed", time2, "Filter2"])


	#disk AnnSQL approach
	start_time = time.time()
	adata_sql_disk.query("SELECT gene_1 FROM X WHERE gene_1 > 0.5 AND gene_2 < 0.5 LIMIT 5")
	time3 = time.time() - start_time
	print(f"disk annsql filter 2: {time3} seconds")

	#disk AnnData approach
	start_time = time.time()
	pd.DataFrame(adata_disk[(adata_disk[:, "gene_1"].X > 0.5) & (adata_disk[:, "gene_2"].X < 0.5), "gene_1"].X[:5], columns=["gene_1"])
	time4 = time.time() - start_time
	print(f"disk annData filter 2: {time4} seconds")

	#Store comparison results
	comparison_records.append([i, "AnnSql On-Disk", time3, "Filter2"])
	comparison_records.append([i, "AnnData Backed", time4, "Filter2"])

	####################################################################################
	# COMPARISON 3 Complex Filter
	####################################################################################
	
	if i <= in_memory_high_filter:
		#memory AnnSQL approach
		start_time = time.time()
		adata_sql_memory.query("SELECT gene_1 FROM X WHERE gene_1 > 0.5 AND gene_2 < 0.5 AND gene_3 > 0.5 LIMIT 5")
		time1 = time.time() - start_time
		print(f"memory annsql filter 3: {time1} seconds")

		#memory AnnData approach
		start_time = time.time()
		pd.DataFrame(adata_memory[(adata_memory[:, "gene_1"].X > 0.5) & (adata_memory[:, "gene_2"].X < 0.5) & (adata_memory[:, "gene_3"].X > 0.5), "gene_1"].X[:5], columns=["gene_1"])
		time2 = time.time() - start_time
		print(f"memory annData filter 3: {time2} seconds")

		#Store comparison results
		comparison_records.append([i, "AnnSql In-Memory", time1, "Filter3"])
		comparison_records.append([i, "AnnData Not-Backed", time2, "Filter3"])

	#disk AnnSQL approach
	start_time = time.time()
	adata_sql_disk.query("SELECT gene_1 FROM X WHERE gene_1 > 0.5 AND gene_2 < 0.5 AND gene_3 > 0.5 LIMIT 5")
	time3 = time.time() - start_time
	print(f"disk annsql filter 3: {time3} seconds")

	#disk AnnData approach
	start_time = time.time()
	pd.DataFrame(adata_disk[(adata_disk[:, "gene_1"].X > 0.5) & (adata_disk[:, "gene_2"].X < 0.5) & (adata_disk[:, "gene_3"].X > 0.5), "gene_1"].X[:5], columns=["gene_1"])
	time4 = time.time() - start_time
	print(f"disk annData filter 3: {time4} seconds")

	#Store comparison results
	comparison_records.append([i, "AnnSql On-Disk", time3, "Filter3"])
	comparison_records.append([i, "AnnData Backed", time4, "Filter3"])


	####################################################################################
	# COMPARISON 4 Complex Filter
	####################################################################################
	
	if i <= in_memory_high_filter:
		#memory AnnSQL approach
		start_time = time.time()
		adata_sql_memory.query("SELECT gene_1 FROM X WHERE gene_1 > 0.5 AND gene_2 < 0.5 AND gene_3 > 0.5 AND gene_4 < 0.5 LIMIT 5")
		time1 = time.time() - start_time
		print(f"memory annsql filter 4: {time1} seconds")

		#memory AnnData approach
		start_time = time.time()
		pd.DataFrame(adata_memory[(adata_memory[:, "gene_1"].X > 0.5) & (adata_memory[:, "gene_2"].X < 0.5) & (adata_memory[:, "gene_3"].X > 0.5) & (adata_memory[:, "gene_4"].X < 0.5), "gene_1"].X[:5], columns=["gene_1"])
		time2 = time.time() - start_time
		print(f"memory annData filter 4: {time2} seconds")

		#Store comparison results
		comparison_records.append([i, "AnnSql In-Memory", time1, "Filter4"])
		comparison_records.append([i, "AnnData Not-Backed", time2, "Filter4"])

	#disk AnnSQL approach
	start_time = time.time()
	adata_sql_disk.query("SELECT gene_1 FROM X WHERE gene_1 > 0.5 AND gene_2 < 0.5 AND gene_3 > 0.5 AND gene_4 < 0.5 LIMIT 5")
	time3 = time.time() - start_time
	print(f"disk annsql filter 4: {time3} seconds")

	#disk AnnData approach
	start_time = time.time()
	pd.DataFrame(adata_disk[(adata_disk[:, "gene_1"].X > 0.5) & (adata_disk[:, "gene_2"].X < 0.5) & (adata_disk[:, "gene_3"].X > 0.5) & (adata_disk[:, "gene_4"].X < 0.5), "gene_1"].X[:5], columns=["gene_1"])
	time4 = time.time() - start_time
	print(f"disk annData filter 4: {time4} seconds")

	#Store comparison results
	comparison_records.append([i, "AnnSql On-Disk", time3, "Filter4"])
	comparison_records.append([i, "AnnData Backed", time4, "Filter4"])

	####################################################################################
	# COMPARISON 5 Complex Filter
	####################################################################################
	
	if i <= in_memory_high_filter:
		#memory AnnSQL approach
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
		end_time = time.time()
		time1 = end_time-start_time
		print(f"memory annsql filter 5: {time1} seconds")

		#memory AnnData approach
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
		end_time = time.time()
		time2 = end_time-start_time
		print(f"memory annData filter 5: {time2} seconds")

		#Store comparison results
		comparison_records.append([i, "AnnSql In-Memory", time1, "Filter5"])
		comparison_records.append([i, "AnnData Not-Backed", time2, "Filter5"])


	#disk AnnSQL approach
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
	end_time = time.time()
	time3 = end_time-start_time
	print(f"disk annsql filter 5: {time3} seconds")

	#disk AnnData approach
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
	print(f"disk annData filter 5: {time4} seconds")

	#Store comparison results
	comparison_records.append([i, "AnnSql On-Disk", time3, "Filter5"])
	comparison_records.append([i, "AnnData Backed", time4, "Filter5"])

	####################################################################################
	# COMPARISON 6 Complex Filter
	####################################################################################
	
	if i <= in_memory_high_filter:
		#memory AnnSQL approach
		start_time = time.time()
		adata_sql_memory.query("SELECT obs.cell_type, AVG(gene_1) as avg_gene_1 FROM X INNER JOIN obs ON X.cell_id = obs.cell_id GROUP BY obs.cell_type")
		end_time = time.time()
		time1 = end_time-start_time
		print(f"memory annsql filter 6: {time1} seconds")

		#memory AnnData approach
		start_time = time.time()
		pd.DataFrame({
			'cell_type': adata_memory.obs['cell_type'],
			'gene_1': adata_memory[:, 'gene_1'].X.flatten()
		}).groupby('cell_type')['gene_1'].mean().reset_index()
		end_time = time.time()
		time2 = end_time-start_time
		print(f"memory annData filter 6: {time2} seconds")

		#Store comparison results
		comparison_records.append([i, "AnnSql In-Memory", time1, "Filter6"])
		comparison_records.append([i, "AnnData Not-Backed", time2, "Filter6"])

	#disk AnnSQL approach
	start_time = time.time()
	adata_sql_disk.query("SELECT obs.cell_type, AVG(gene_1) as avg_gene_1 FROM X INNER JOIN obs ON X.cell_id = obs.cell_id GROUP BY obs.cell_type")
	end_time = time.time()
	time3 = end_time-start_time
	print(f"disk annsql filter 6: {time3} seconds")

	#disk AnnData approach
	start_time = time.time()
	pd.DataFrame({
		'cell_type': adata_disk.obs['cell_type'],
		'gene_1': adata_disk[:, 'gene_1'].X.flatten()
	}).groupby('cell_type')['gene_1'].mean().reset_index()
	end_time = time.time()
	time4 = end_time-start_time
	print(f"disk annData filter 6: {time4} seconds")

	#Store comparison results
	comparison_records.append([i, "AnnSql On-Disk", time3, "Filter6"])
	comparison_records.append([i, "AnnData Backed", time4, "Filter6"])

	#clean up memory
	adata_memory = None
	adata_sql_memory = None
	adata_disk = None
	adata_sql_disk = None
	gc.collect()


# Convert the list to DataFrame after the loop
comparisons = pd.DataFrame(comparison_records, columns=["size", "type", "runtime", "filter"])

#ntural log scale for runtime
comparisons['runtime_log'] = np.log(comparisons['runtime'])

#log10 scale for runtime
comparisons['runtime_log10'] = np.log10(comparisons['runtime'])

#add log scale for size
comparisons['size_log'] = np.log(comparisons['size'])

#add log10 scale for size
comparisons['size_log10'] = np.log10(comparisons['size'])

#store the data
comparisons.to_csv("../results/comparisons.csv", index=False)

#load the data
comparisons = pd.read_csv("../results/comparisons.csv")

#set the colors of the plots (ansql1, anndata in-mem, ansql2, anndata backed)
colors = ["#07b88e", "#a4a6a4", "#07b88e", "#a4a6a4"]

#plot aggregation of all filters runtime
sns.set(style="whitegrid")
plt.figure(dpi=1200)
sns.lineplot(data=comparisons, x='size', y='runtime_log', hue='type', palette=colors)
plt.ylabel("")
plt.xlabel("")
plt.xticks(rotation=45)
plt.xticks(np.arange(0, 275000, 25000))
plt.xlim(0, 250000)
plt.yticks(np.log([0.01,0.1, 1, 10, 100, 1000]), [0.01,0.1, 1, 10, 100, 1000])
plt.xticks(fontsize=20, fontweight='bold')
plt.yticks(fontsize=20, fontweight='bold')	
plt.grid(False)
plt.legend().remove()
plt.show()

#mean runtime for each type with a 250,000 cell library 
comparisons[comparisons['size'] == 250000].groupby('type')['runtime'].mean()

#comparision 1
plt.figure(dpi=950)
sns.lineplot(data=comparisons[comparisons['filter'] == 'Filter1'], x='size', y='runtime_log', hue='type',linewidth=4, palette=colors)
plt.ylabel("")
plt.xlabel("")
plt.yticks(np.log([0.01,0.1, 1, 10,100]), [0.01,0.1, 1, 10,100])
plt.xticks(rotation=45)
plt.legend().remove()
plt.grid(False)
plt.xticks([])

#comparision 2
plt.figure(dpi=950)
sns.lineplot(data=comparisons[comparisons['filter'] == 'Filter2'], x='size', y='runtime_log', hue='type',linewidth=4, palette=colors)
plt.ylabel("")
plt.xlabel("")
plt.yticks(np.log([0.01,0.1, 1, 10,100]), [0.01,0.1, 1, 10,100])
plt.xticks(rotation=45)
plt.legend().remove()
plt.grid(False)
plt.xticks([])

#comparision 3
plt.figure(dpi=950)
sns.lineplot(data=comparisons[comparisons['filter'] == 'Filter3'], x='size', y='runtime_log', hue='type',linewidth=4, palette=colors)
plt.ylabel("")
plt.xlabel("")
plt.yticks(np.log([0.01,0.1, 1, 10,100]), [0.01,0.1, 1, 10,100])
plt.xticks(rotation=45)
plt.legend().remove()
plt.grid(False)
plt.xticks([])

#comparision 4
plt.figure(dpi=950)
sns.lineplot(data=comparisons[comparisons['filter'] == 'Filter4'], x='size', y='runtime_log', hue='type',linewidth=4, palette=colors)
plt.ylabel("")
plt.xlabel("")
plt.yticks(np.log([0.01,0.1, 1, 10,100]), [0.01,0.1, 1, 10,100])
plt.xticks(rotation=45)
plt.legend().remove()
plt.grid(False)
plt.xticks([])

#comparision 5
plt.figure(dpi=950)
sns.lineplot(data=comparisons[comparisons['filter'] == 'Filter5'], x='size', y='runtime_log', hue='type',linewidth=4, palette=colors)
plt.ylabel("")
plt.xlabel("")
plt.yticks(np.log([0.01,0.1, 1, 10,100]), [0.01,0.1, 1, 10,100])
plt.xticks(rotation=45)
plt.legend().remove()
plt.grid(False)
plt.xticks([])

#comparision 6
plt.figure(dpi=950)
sns.lineplot(data=comparisons[comparisons['filter'] == 'Filter6'], x='size', y='runtime_log', hue='type',linewidth=4, palette=colors)
plt.ylabel("")
plt.xlabel("")
plt.yticks(np.log([0.01,0.1, 1, 10,100]), [0.01,0.1, 1, 10,100])
plt.xticks(rotation=45)
plt.legend().remove()
plt.grid(False)