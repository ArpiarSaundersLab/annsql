#import libraries
import scanpy as sc
import pandas as pd
from AnnSQL import AnnSQL

#open the database
asql = AnnSQL(db="../db/random_data_10000.asql")

#show tables
asql.show_tables()

#query the obs to see what categories are available
asql.query("SELECT * FROM obs")

#build the meta cells table
asql.build_meta_cells(primary_cluster="cell_type", print_progress=True)

#query the database to see the meta_cells table
asql.query("SELECT * FROM meta_cells")

#generate a csv file
asql.query("SELECT * FROM meta_cells").to_csv("meta_cells.csv", index=False)