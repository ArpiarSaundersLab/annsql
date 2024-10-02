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

#dataset sizes to generate
dataset_sizes = [1000, 5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000, 75000, 100000, 250000]

#iterates and creates data_sizes defined above
for i in dataset_sizes:

	print(f"Running for {i}")

	#the amount of columns (genes)
	columns = 10000

	#generate a cellxgene object using random data
	adata = sc.AnnData(X=np.random.rand(i,columns), 
						obs=pd.DataFrame(index=[f"cell_{i}" for i in range(i)]), 
						var=pd.DataFrame(data=[f"gene_{i}" for i in range(columns)],columns=["gene_name"]))
	adata.var.index = [f"gene_{i}" for i in range(adata.shape[1])]
	adata.obs["cell_type"] = np.random.choice(["A","B","C"], adata.shape[0])

	#save the object
	adata.write("../data/random_data_"+str(i)+".h5ad")

	#open the object
	adata = sc.read("../data/random_data_"+str(i)+".h5ad", backed="r")

	#make the database using backed mode
	start_time = time.time()
	MakeDb(adata, db_path="../db/", db_name="random_data_"+str(i), create_all_indexes=False, convenience_view=False)
	print("--- %s seconds ---" % (time.time() - start_time))

	#clear mem
	gc.collect()
	adata = None
