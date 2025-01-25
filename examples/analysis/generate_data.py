#import libraries
import scanpy as sc
from AnnSQL import AnnSQL
import time
import pandas as pd
from AnnSQL.MakeDb import MakeDb
import os
import numpy as np
import gc
import seaborn as sns
import matplotlib.pyplot as plt
from scsim import scsim


#dataset sizes to generate
dataset_sizes = [1000, 5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000, 75000, 100000, 250000]

#iterates and creates data_sizes defined above
for i in dataset_sizes:

	print(f"Running for {i}")

	columns = 10000

	simulator = scsim(ngenes=columns, ncells=i, ngroups=5, libloc=7.64, libscale=0.78,
					mean_rate=7.68,mean_shape=0.34, expoutprob=0.00286,
					expoutloc=6.15, expoutscale=0.49,
					diffexpprob=0.025, diffexpdownprob=0., diffexploc=1.0, diffexpscale=1.0,
					bcv_dispersion=0.448, bcv_dof=22.087, ndoublets=0,
					nproggenes=400, progdownprob=0., progdeloc=1.0,
					progdescale=1.0, progcellfrac=0.35, proggoups= list(range(1, int(10/3)+1)),
					minprogusage=.1, maxprogusage=.7, seed=21)


	start_time = time.time()
	simulator.simulate()
	end_time = time.time()
	print('minutes elapsed for seed:',((end_time-start_time)/60))

	#generate a cellxgene object using random data
	adata = sc.AnnData(X=simulator.counts, 
						obs=simulator.cellparams, 
						)
	adata.var.index = [f"gene_{i}" for i in range(adata.shape[1])]
	adata.obs["cell_type"] = adata.obs["group"]
	adata.obs = adata.obs.drop(columns=["group","has_program","program_usage","libsize"])

	#save the object
	adata.write("../data/random_data_"+str(i)+".h5ad")

	#open the object
	adata = sc.read("../data/random_data_"+str(i)+".h5ad", backed="r")

	#make the database using backed mode
	os.system(f"rm ../db/random_data_{i}.asql")
	start_time = time.time()
	MakeDb(adata, db_path="../db/", db_name="random_data_"+str(i), create_all_indexes=False, convenience_view=False)
	print("--- %s seconds ---" % (time.time() - start_time))

	#clear mem
	gc.collect()
	adata = None