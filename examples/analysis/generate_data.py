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
high_range = 250001

#iterate from 1k to x
for i in range(low_range, high_range):
	if i % 2000 == 0:		
		#skip from 100k to 250k. Too much data to store locally for all the steps.
		if i > 100000 and i < 250000:
			continue

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
		adata.write("../data/random_data_"+str(i)+".h5ad")

		#open the object
		adata = sc.read("../data/random_data_"+str(i)+".h5ad", backed="r")

		#make the database using backed mode
		MakeDb(adata, db_path="../db/", db_name="random_data_"+str(i), create_all_indexes=False, convenience_view=False)