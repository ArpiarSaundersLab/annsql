import argparse
import os, sys
from scsim import scsim
import numpy as np
import time
import anndata as ad
import scanpy as sc
import pandas as pd


simulator = scsim(ngenes=10000, ncells=1000, ngroups=10, libloc=7.64, libscale=0.78,
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
adata.write("test.h5ad")

#open the object
adata = sc.read("test.h5ad")
adata.to_df()
adata.obs
