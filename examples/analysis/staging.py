import scanpy as sc
import pandas as pd
import numpy as np
import time
import duckdb
import seaborn as sns
import matplotlib.pyplot as plt
from numpy.linalg import eig
import gc
import polars as pl
from memory_profiler import memory_usage
import os
from AnnSQL import AnnSQL
from AnnSQL.MakeDb import MakeDb
from scipy import stats
import warnings
warnings.filterwarnings('ignore')
from scipy.sparse import issparse

#TODO debug duckdb memory leak on daabases with >30k genes. 
def run_test(chunk):

	# Df with 1k rows and 30k columns. 
	# This is a very similar matrix to a cell x gene matrix 
	# which typically will have >50k rows and >30k columns.
	n_rows = 10
	n_columns = 30000
	df = np.random.rand(n_rows, n_columns)
	df = pd.DataFrame(df)
	df.columns = [f"gene_{i}" for i in range(n_columns)]
	df.index = [f"{i}" for i in range(n_rows)]
	df.index.name = 'cell_id'
	df = df.reset_index()
	
	#open the door
	#conn = duckdb.connect(database='test.db', config={'memory_limit':'3GB'})
	conn = duckdb.connect(database='test.db')

	#does the table exist? If not, create it.
	if "X" not in conn.sql("SHOW TABLES;").df()["name"].to_list():
		conn.execute("CREATE OR REPLACE TABLE X ({})".format(', '.join([f"{col} FLOAT" for col in df.columns])))

	#display the memory occupied by df	(should be roughly 8 bytes x 30k genes x 1k cells = 240Mb)
	print("\tTotal df memory: ",df.memory_usage(index=True).sum() / 1024**2,"Mb")

	#insert the data (will be iteratived to simulate chunks)
	start_time = time.time()
	conn.execute("INSERT INTO X SELECT * FROM df")
	print(f"\tInsert time: {time.time() - start_time} seconds")

	#total database size
	print("\tTotal database size: ",os.path.getsize('test.db') / 1024**2,"Mb")

	#close the door
	conn.close()


if __name__ == "__main__":
	#run the test with 5 chunks
	for i in range(1, 6):
		print(f"Chunk: {i}")
		mem_usage = memory_usage((run_test, (i,), {}))
		print('\tMax memory usage: {:.2f} MB'.format(max(mem_usage)))
