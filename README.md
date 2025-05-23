<p align="center">
	<img src="examples/images/logo.png" width=500>
</p>

<br>

<p align="center">
	<a href="https://github.com/ArpiarSaundersLab/annsql/tree/main/tests"><img src="https://img.shields.io/badge/build-passing-brightgreen"></a>
	<a href="https://img.shields.io/github/v/release/ArpiarSaundersLab/annsql"><img src="https://img.shields.io/github/v/release/ArpiarSaundersLab/annsql"></a>
	<a href="https://static.pepy.tech/badge/annsql/month"><img src="https://static.pepy.tech/badge/annsql/month"></a>
	<a href="https://static.pepy.tech/badge/annsql"><img src="https://static.pepy.tech/badge/annsql"></a>
</p>

<br>

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

############################################################################
# Extended AnnSQL methods (See: https://docs.annsql.com/preprocessing)
# These methods are either SQL based or Python/SQL hybrid implementations.
############################################################################

#basic QC on the dataset
asql.calculate_total_counts()
asql.filter_by_cell_counts(min_cell_count=1000, max_cell_count=50000)
asql.filter_by_gene_counts(min_gene_counts=100, max_gene_counts=10000)

#normalize & log umi counts
asql.expression_normalize(total_counts_per_cell=10000)
asql.expression_log(log_type="LN")

#select highly variable genes
asql.calculate_variable_genes(save_var_names=True, top_variable_genes=1000)

#run pca
asql.calculate_pca(n_pcs=50, top_variable_genes=1000, zero_center=False)

#umap, cluster, then and plot.
asql.calculate_umap()
asql.calculate_leiden_clusters(resolution=0.25, n_neighbors=5)
asql.plot_umap(color_by="leiden_clusters", annotate=True)
```

<br>
<br>


## Reference
Kenny Pavan, Arpiar Saunders, AnnSQL: A Python SQL-based package for fast large-scale single-cell genomics analysis using minimal computational resources<br>
Bioinformatics Advances, 2025; vbaf105, [https://doi.org/10.1093/bioadv/vbaf105](https://doi.org/10.1093/bioadv/vbaf105)
<br>
<br>