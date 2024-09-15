#import libraries
import scanpy as sc
from AnnSQL import AnnSQL
import time
import pandas as pd
from MakeDb import MakeDb
import os
import numpy as np
import gc
import seaborn as sns
import matplotlib.pyplot as plt

#set the iterations
low_range = 1000
high_range = 102000

#iterate from 1k to x
for i in range(low_range, high_range):
	if i % 2000 == 0:
	#if i == 2000:
		print(f"Running for {i} rows")

		#the amount of columns (genes)
		columns = 10000

		#generate a cellxgene object using random data
		adata = sc.AnnData(X=np.random.rand(i,columns), 
							obs=pd.DataFrame(index=[f"cell_{i}" for i in range(i)]), 
							var=pd.DataFrame(data=[f"gene_{i}" for i in range(columns)],columns=["gene_name"]))
		adata.var.index = [f"gene_{i}" for i in range(adata.shape[1])]
		adata.obs["cell_type"] = np.random.choice(["A","B","C","D"], adata.shape[0])

		#save the object
		adata.write("data/random_data_"+str(i)+".h5ad")

		#open the object
		adata = sc.read("data/random_data_"+str(i)+".h5ad", backed="r")

		#make the database using backed mode
		MakeDb(adata, db_path="db/", db_name="random_data_"+str(i), create_all_indexes=False, add_uns=False, convenience_view=False)


#comparisons
comparisons = {}

#iterate from 1k to x
for i in range(low_range, high_range):
	if i % 2000 == 0:
		if not os.path.exists("data/random_data_"+str(i)+".h5ad"):
			continue

		print(f"Running for {i}")
		#open the object
		adata = sc.read("data/random_data_"+str(i)+".h5ad", backed="r")
		#adata = sc.read("data/random_data_"+str(i)+".h5ad")

		#open the database
		adata_sql = AnnSQL(db="db/random_data_"+str(i)+".asql")

		#add to the comparisons
		comparisons[i] = []

		####################################################################################
		#COMPARISON 1 Simple Filter
		####################################################################################
		#Query approach
		start_time = time.time()
		adata_sql.query("SELECT gene_1 FROM X WHERE gene_1 > 0.5 LIMIT 5")
		end_time = time.time()
		time1 = end_time-start_time	
		print(f"query1 time: {end_time-start_time}")

		#AnnData approach
		start_time = time.time()
		pd.DataFrame(adata[adata[:,"gene_1"].X > 0.5,"gene_1"].X[:5], columns=["gene_1"])
		end_time = time.time()
		time2 = end_time-start_time
		print(f"adata1 time: {end_time-start_time}")
		comparisons[i].append([time1,time2])

		####################################################################################
		#COMPARISON 2 Multiple Filter
		####################################################################################
		#Query approach
		start_time = time.time()
		adata_sql.query("SELECT gene_1, gene_2 FROM X WHERE gene_1 > 0.5 AND gene_2 < 0.3 LIMIT 5")
		end_time = time.time()
		print(f"query2 time: {end_time-start_time}")
		time1 = end_time-start_time	

		#must be copied to memory to compare. backed mode doesn't work
		start_time = time.time()
		filtered_rows = (adata[:, "gene_1"].X > 0.5) & (adata[:, "gene_2"].X < 0.3)
		if adata.isbacked:
			filtered_data = adata[filtered_rows, :].to_memory() #for backed mode	
		else:
			filtered_data = adata[filtered_rows, :]
		filtered_result = filtered_data[:, ["gene_1", "gene_2"]].X[:5]
		pd.DataFrame(filtered_result, columns=["gene_1", "gene_2"]) #to be fair, must be df
		end_time = time.time()
		time2 = end_time-start_time
		print(f"adata2 time: {end_time-start_time}")
		comparisons[i].append([time1,time2])

		####################################################################################
		#COMPARISON 3 Complex Filter
		####################################################################################
		#Query approach
		start_time = time.time()
		adata_sql.query("""SELECT 
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
		print(f"query3 time: {end_time-start_time}")

		#AnnData approach
		start_time = time.time()
		gene_1_vals = adata[:, "gene_1"].X.flatten()
		gene_2_vals = adata[:, "gene_2"].X.flatten()
		gene_3_vals = adata[:, "gene_3"].X.flatten()
		gene_4_vals = adata[:, "gene_4"].X.flatten()
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
		print(f"adata3 time: {end_time-start_time}")

		comparisons[i].append([time1,time2])

		####################################################################################
		#COMPARISON 4 Complex Filter
		####################################################################################
		start_time = time.time()
		adata_sql.query("SELECT obs.cell_type, AVG(gene_1) as avg_gene_1 FROM X INNER JOIN obs ON X.cell_id = obs.cell_id GROUP BY obs.cell_type")
		end_time = time.time()
		time1 = end_time-start_time
		print(f"query4: {end_time-start_time}")

		start_time = time.time()
		pd.DataFrame({
			'cell_type': adata.obs['cell_type'],
			'gene_1': adata[:, 'gene_1'].X.flatten()
		}).groupby('cell_type')['gene_1'].mean().reset_index()
		end_time = time.time()
		time2 = end_time-start_time
		print(f"adata4: {end_time-start_time}")
		comparisons[i].append([time1,time2])


		#clear the adata object
		adata = None
		adata_sql=None
		gc.collect()



#save the comparisons in a dataframe
data = []
for group, values in comparisons.items():
	inc=1
	for value in values:
		data.append([group, value[0], 'AnnSQL', inc])
		data.append([group, value[1], 'AnnData', inc])
		inc+=1
df = pd.DataFrame(data, columns=['size', 'time', 'type', 'filter'])
df['size'] = df['size'].astype(str)

#save the dataframe
#df.to_csv("results/runtime_comparison_backed.csv", index=False)

#open the data
df = pd.read_csv("results/runtime_comparison_not_backed.csv")
df = pd.read_csv("results/runtime_comparison_backed.csv")

#log scale time for better visualization
df['time_log'] = np.log(df['time'])

#plot the data as a lineplot
sns.set(style="whitegrid")
sns.lineplot(data=df, x='size', y='time', hue='type')
plt.ylabel("runtime seconds")
plt.xlabel("Cell Library Count")
plt.xticks(rotation=90)
plt.title("Query vs AnnData (Backed) \nALL FILTERS Runtime Combined")
plt.show()

#create a barplot for each filter
for i in range(1,5):
	df_filter = df[df['filter'] == i]
	sns.set(style="whitegrid")
	sns.barplot(data=df_filter, x='size', y='time', hue='type')
	plt.ylabel("runtime (sec)")
	plt.xlabel("Cell Library Count")
	plt.xticks(rotation=90)
	plt.title(f"Filter {i}: Query vs AnnData")
	plt.show()