#import libraries
import scanpy as sc
import pandas as pd
from AnnSQL import AnnSQL
from AnnSQL.MakeDb import MakeDb
import time


#open the database
asql = AnnSQL(db="../db/pbmc3k.asql")

asql.show_tables()

asql.query("Describe X")


# #filter the ENSMUSG00000070880 the same way as the query but in scanpy
# start_time = time.time()

# print("Filtering runtime no df:", time.time() - start_time)
