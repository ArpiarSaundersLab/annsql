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
For larger datasets, AnnSQL can create a local database (asql) from the AnnData object. This database is stored on-disk and queried.
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

## Advanced Queries and usage
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

## Accessing and processing 4.4 million cells on a laptop
To illustrate how AnnSQL can be used to access atlas sized datasets on a local computer, we examine the single nuclei dataset presented in "The molecular cytoarchitecture of the adult mouse brain" by <a href='https://www.nature.com/articles/s41586-023-06818-7' target="_blank">Langlieb et al 2023</a>. First, we opened the <a href='https://docs.braincelldata.org/downloads/index.html' target="_blank">atlas</a> AnnData object in backed mode and created a `asql` database using the `MakeDb` class provided with AnnSQL. Next, we performed some basic querying of the data to return subsets. We then calculated total counts per gene which we accomplished entirely in SQL; even with the non-optimized schema. Next, we blended Python and SQL to normalize and log the counts per library. Lastly, we calculated highly variable genes in the entire dataset using two SQL queries which: (1) provide a list of all gene names in the X table, then (2) use those gene names to calculate the variance for each gene and return a list of the top 2000. Our results demonstrate AnnSQL is a capable tool for basic (and possibly more advanced) analyses of atlas scale datasets. 

```python
#import libraries
from MakeDb import MakeDb
from AnnSQL import AnnSQL

#load the atlas dataset in backed mode
adata = sc.read_h5ad("Macosko_Mouse_Atlas_Single_Nuclei.Use_Backed.h5ad", backed="r+")

#build the asql database | Runtime 7hr 10min
MakeDb(adata=adata, db_name="Macosko_Mouse_Atlas", db_path="../db/", layers=["X", "obs"])

#query example | Runtime: 0.24sec
asql.query("SELECT ENSMUSG00000070880 FROM X WHERE ENSMUSG00000070880 > 0")

#count the number of cells in each cluster | Runtime: 0.35sec
asql.query("SELECT ClusterNm, COUNT(cell_id) AS num_cells FROM obs GROUP BY ClusterNm ORDER BY num_cells DESC")

#determine the total counts per cell library | Runtime: 4min 30sec
asql.calculate_total_counts(chunk_size=950)

#normalize umi counts to 10k per cell | Runtime: 1hr 48mins
asql.expression_normalize(total_counts_per_cell=1e4, chunk_size=300) 

#log scale the normalized counts | Runtime: 59mins 13sec
asql.expression_log(log_type="LN", chunk_size=250)

```

### Laptop system details for runtimes above
- **Memory:**                                      40.0 GiB
- **Processor:**                                   12th Gen Intel® Core™ i7-1255U × 12
- **Disk Capacity:**                               1.0 TB
- **OS:**                                     	   Ubuntu 24.04.1 LTS
- **Python Version:**                              3.12


<br>
<br>


## Reference
AnnSQL: A Python SQL-based package for large-scale single-cell genomics analysis on a laptop<br />
Kenny Pavan, Arpiar Saunders<br />
bioRxiv 2024.11.02.621676; [doi: https://doi.org/10.1101/2024.11.02.621676](https://www.biorxiv.org/content/10.1101/2024.11.02.621676)

<br>
<br>