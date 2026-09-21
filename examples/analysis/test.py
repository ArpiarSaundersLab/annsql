from AnnSQL import AnnSQL
from AnnSQL.MakeDb import MakeDb
import scanpy as sc
import os
import scipy.sparse

#open a h5ad file
adata = sc.read_h5ad("../data/random/data_1000.h5ad")

#visual check if dense or sparse
print(adata.X)

#Check if the data is sparse using scipy
is_sparse = scipy.sparse.issparse(adata.X)
print("Sparse:", is_sparse)


