#import libraries
import scanpy as sc
import pandas as pd
from AnnSQL import AnnSQL
from AnnSQL.MakeDb import MakeDb
from statsmodels.stats.multitest import multipletests
import time

#open the database
asql = AnnSQL(db="../db/random_data_1000.asql")
asql.filter_by_cell_counts(min_cell_count=4000, max_cell_count=5000)