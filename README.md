<center><img src="examples/images/logo.png" width=500></center>
<br />

# Query AnnData Objects with SQL
The Python based AnnSQL package enables SQL-based queries on [AnnData](https://anndata.readthedocs.io/en/latest/) objects, returning results as either a [Pandas](https://pandas.pydata.org/) DataFrame, an AnnData object, or a [Parquet](https://parquet.apache.org/) file that can easily be imported into a variety of data analysis tools. Behind the scenes, AnnSQL converts the layers of an AnnData object into a relational [DuckDB](https://duckdb.org/) database. Each layer is stored as an individual table, allowing for simple or complex SQL queries, including table joins.

## Features
- Query AnnData with **SQL**.
- **Fast** for complex queries and aggregative functions.
- Return query results as **Pandas** Dataframes, **Parquet** files, or **AnnData** objects.
- Create in-memory or on-disk databases directly from AnnData objects.
- Open AnnSQL databases in **R**. *No conversions necessary*. <a href="https://docs.annsql.com/R_usage/" target="_blank">Learn more</a>

<br>

## Full Documentation

<h3> <a href="https://docs.annsql.com">docs.annsql.com</a></h3>

<br>

## Quick Setup
```
pip install annsql
```

## Basic Usage (In-Memory)
Ideal for smaller datasets.
```python
from AnnSQL import AnnSQL
import scanpy as sc

#read sample data
adata = sc.datasets.pbmc68k_reduced()

#instantiate the AnnData object (you may also pass a h5ad file to the adata parameter)
asql = AnnSQL(adata=adata)

#query the expression table. Returns Pandas Dataframe by Default
asql.query("SELECT * FROM adata LIMIT 10")
```


## Basic Usage (On-Disk)
For larger datasets, AnnSQL can create a local database (asql) from the AnnData object. This database is stored on-disk, can be queried, and is persistent.
```python
import scanpy as sc
from AnnSQL import AnnSQL
from AnnSQL.MakeDb import MakeDb

#read sample data
adata = sc.datasets.pbmc68k_reduced()

#build the AnnSQL database
MakeDb(adata=adata, db_name="pbmc3k_reduced", db_path="db/")

#open the AnnSQL database
asql = AnnSQL(db="db/pbmc3k_reduced.asql")

#query the expression table
asql.query("SELECT * FROM adata LIMIT 5")
```

## Advanced Queries and Usage
```python
from AnnSQL import AnnSQL
import scanpy as sc

#read sample data
adata = sc.datasets.pbmc68k_reduced()

#pass the AnnData object to the AnnSQL class
asql = AnnSQL(adata=adata)

#group and count all labels
asql.query("SELECT obs.bulk_labels, COUNT(*) FROM obs GROUP BY obs.bulk_labels")

#take the log10 of a value
asql.query("SELECT LOG10(HES4) FROM X WHERE HES4 > 0")

#sum all gene counts | Memory intensive | See method calculate_gene_counts for chunked approach.
asql.query("SELECT SUM(COLUMNS(*)) FROM (SELECT * EXCLUDE (cell_id) FROM X)")

#taking the correlation of genes ITGB2 and SSU72 in dendritic cells that express either gene > 0
asql.query("SELECT corr(ITGB2,SSU72) as correlation FROM adata WHERE bulk_labels = 'Dendritic' AND (ITGB2 > 0 OR SSU72 >0)")
```

<br>
<br>


## Reference
AnnSQL: A Python SQL-based package for large-scale single-cell genomics analysis on a laptop<br />
Kenny Pavan, Arpiar Saunders<br />
bioRxiv 2024.11.02.621676; [doi: https://doi.org/10.1101/2024.11.02.621676](https://www.biorxiv.org/content/10.1101/2024.11.02.621676)

<br>
<br>