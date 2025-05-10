<center><img src="https://github.com/ArpiarSaundersLab/annsql/raw/main/examples/images/logo.png" width=500></center>
<br />
The Python based AnnSQL package enables SQL-based queries on [AnnData](https://anndata.readthedocs.io/en/latest/) objects, returning results as either a [Pandas](https://pandas.pydata.org/) DataFrame, an AnnData object, or a [Parquet](https://parquet.apache.org/) file that can easily be imported into a variety of data analysis tools. Behind the scenes, AnnSQL converts the layers of an AnnData object into a relational [DuckDB](https://duckdb.org/) database. Each layer is stored as an individual table, allowing for simple or complex SQL queries, including table joins.

## Features
- Query AnnData with **SQL**.
- **Fast** for complex queries and aggregative functions.
- Return query results as **Pandas** Dataframes, **Parquet** files, or **AnnData** objects.
- Create in-memory or on-disk databases directly from AnnData objects.
- Open AnnSQL databases in **R**. *No conversions necessary*. <a href="inteRoperability/">Learn more</a>

## Installation
**Note:** Higher memory consumption using Apple M-Series is expected when building AnnSQL databases. Please use the *make_buffer_file* parameter when using the *MakeDb* class if you experience memory related issues. 
```
pip install annsql
```
<br>

## Basic Usage 
### In-Memory
Ideal for smaller datasets.
```python
from AnnSQL import AnnSQL
import scanpy as sc

#read sample data
adata = sc.datasets.pbmc68k_reduced()

#instantiate the AnnData object (you may also pass a h5ad file to the adata parameter)
asql = AnnSQL(adata=adata)
```
<br>

### On-Disk
For larger datasets, AnnSQL can create a local database (asql) from the AnnData object. This database is stored on-disk and queried. Storage requirements are similar to the original AnnData h5ad filesize; however, complex aggregative functions can now benefit from the DuckDb engine with full SQL support. 
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
<br>

## Entity Relationship Diagram
Using the [Scanpy](https://scanpy.readthedocs.io/) sample dataset, "Processed 3k PBMCs from 10x Genomics," the following ERD was generated from the DuckDB database constructed via **AnnSQL**, based on the corresponding AnnData object. **Note:** The database structure is not optimized for performance. Instead, the tables are designed to closely mirror the familiar structure of the AnnData object for ease of use.
```bash
AnnData object with n_obs × n_vars = 700 × 765
	obs: 'bulk_labels', 'n_genes', 'percent_mito', 'n_counts', 'S_score', 'phase', 'louvain'
	var: 'n_counts', 'means', 'dispersions', 'dispersions_norm', 'highly_variable'
	uns: 'bulk_labels_colors', 'louvain', 'louvain_colors', 'neighbors', 'pca', 'rank_genes_groups'
	obsm: 'X_pca', 'X_umap'
	varm: 'PCs'
	obsp: 'distances', 'connectivities'
```

<img src="https://github.com/ArpiarSaundersLab/annsql/raw/main/examples/images/erd.png">

<br>

## Runtime Comparision
There are two key reasons to use **AnnSQL**: (1) if you prefer SQL's expressive syntax for filtering and querying your data, or (2) if you're working with datasets that exceed memory limits and require loading AnnData in backed mode. Using backed mode in AnnData can limit available functions, especially aggregate operations, and slow down data access. **AnnSQL** offers a solution by enabling SQL-style queries that may perform more efficiently in these scenarios. Below are rough runtime comparisons between AnnData and AnnSQL after a database has been built. Running AnnSQL locally for datasets that are larger than memory, that would typically require AnnData in backed mode see substantial runtime improvements for a variety of filtering operations. 
<img src="https://github.com/ArpiarSaundersLab/annsql/raw/main/examples/images/comparision.png">

<br>
<br>
<br>


## Citation
Kenny Pavan, Arpiar Saunders, AnnSQL: A Python SQL-based package for fast large-scale single-cell genomics analysis using minimal computational resources<br>
Bioinformatics Advances, 2025; vbaf105, [https://doi.org/10.1093/bioadv/vbaf105](https://doi.org/10.1093/bioadv/vbaf105)

<br>
<br>