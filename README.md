<center><img src="examples/images/logo.png" width=500></center>
<br />

# Query AnnData Objects with SQL
The AnnSQL package enables SQL-based queries on [AnnData](https://anndata.readthedocs.io/en/latest/) objects, returning results as either a [Pandas](https://pandas.pydata.org/) DataFrame, an AnnData object, or a [Parquet](https://parquet.apache.org/) file that can easily be imported into a variety of data analysis tools. Behind the scenes, AnnSQL converts the layers of an AnnData object into a relational [DuckDB](https://duckdb.org/) database. Each layer is stored as an individual table, allowing for simple or complex SQL queries, including table joins.

To get started, see the usage notes below or refer to the example scripts in the `examples` directory.

## Features
- Query AnnData with **SQL**.
- Return query results as **Pandas** Dataframes, **Parquet** files, or **AnnData** objects.
- Create in-memory or on-Disk databases directly from AnnData objects.
- **Fast** for complex queries and aggregative functions.

## Installation
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

#instantiate the AnnData object
adata_sql = AnnSQL(adata=adata)

#query the expression table
adata_sql.query("SELECT * FROM X")

#query the observation table
adata_sql.query("SELECT * FROM obs")

#query the join of 'X' and 'obs' table
adata_sql.query("SELECT * FROM adata")
```


## Basic Usage (On-Disk)
For larger datasets, AnnSQL can create a local database (asql) from the AnnData object. This database is stored on-disk and queried. Storage requirements are similar to the original AnnData h5ad filesize; however, complex aggregative functions can now benefit from the DuckDb engine with full SQL support. Please see <a href="#">manuscript</a> for benchmarks.
```python
import scanpy as sc
from MakeDb import MakeDb

#read sample data
adata = sc.datasets.pbmc68k_reduced()

#build the AnnSQL database
MakeDb(adata=adata, db_name="pbmc3k_reduced", db_path="db/")

#open the AnnSQL database
adata_sql = AnnSQL(db="db/pbmc3k_reduced.asql")

#query the expression table
adata_sql.query("SELECT * FROM adata")
```

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

<img src="examples/images/erd.png">

## Advanced Queries and usage
```python
from AnnSQL import AnnSQL
import scanpy as sc

#read sample data
adata = sc.datasets.pbmc68k_reduced()

#pass the AnnData object to the AnnSQL class
adata_sql = AnnSQL(adata=adata)

#group and count all labels
adata_sql.query("SELECT obs.bulk_labels, COUNT(*) FROM obs GROUP BY obs.bulk_labels")

#take the log10 of a value
adata_sql.query("SELECT LOG10(HES4) FROM X WHERE HES4 > 0")

#sum all gene counts
adata_sql.query("SELECT SUM(COLUMNS(*)) FROM (SELECT * EXCLUDE (cell_id) FROM X)")

#taking the correlation of genes ITGB2 and SSU72 in dendritic cells that express either gene > 0
adata_sql.query("SELECT corr(ITGB2,SSU72) as correlation FROM adata WHERE bulk_labels = 'Dendritic' AND (ITGB2 > 0 OR SSU72 >0)")

```


## AnnSQL Class

<table>
  <thead>
    <tr>
      <th>Method</th>
      <th>Parameters</th>
      <th>Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><code>__init__(adata, db)</code></td>
      <td>
          <li><code>adata</code>: AnnData object (optional)</li>
          <li><code>db</code>: Path to DuckDB database (optional)</li>
          <li><code>create_all_indexes</code>: Boolean (optional. default: False). <i>Warning: Runtime can be significant when building.</i></li>
      </td>
      <td>Initializes the AnnSQL object. Requires either a AnnData object (<code>adata</code>) or a DuckDB database path (<code>db</code>).</td>
    </tr>
    <tr>
      <td><code>query(query, return_type)</code></td>
      <td>
          <li><code>query</code>: SQL query string</li>
          <li><code>return_type</code>: 'pandas', 'adata', or 'parquet' (default: 'pandas')</li>
      </td>
      <td>Executes a SELECT SQL query. Returns results as a pandas DataFrame, AnnData object, or parquet file.</td>
    </tr>
    <tr>
      <td><code>query_raw(query)</code></td>
      <td>
          <li><code>query</code>: SQL query string</li>
      </td>
      <td>Executes a raw SQL query without restrictions on the type of query. Returns the raw result of the query.</td>
    </tr>
    <tr>
      <td><code>update_query(query)</code></td>
      <td>
          <li><code>query</code>: SQL update, delete, or insert query</li>
      </td>
      <td>Executes an UPDATE, DELETE, or INSERT SQL query. Raises an error if a SELECT query is detected.</td>
    </tr>
    <tr>
      <td><code>show_tables()</code></td>
      <td><li>None</li></td>
      <td>Displays the list of all tables in the DuckDB instance.</td>
    </tr>
    <tr>
      <td><code>show_settings()</code></td>
      <td><li>None</li></td>
      <td>Returns the current DuckDB settings in a pandas DataFrame format.</td>
    </tr>
    <tr>
      <td><code>export_parquet()</code></td>
      <td><li>None</li></td>
      <td>Exports all tables in the DuckDB database to individual Parquet files, saved in the <code>parquet_files</code> folder.</td>
    </tr>
  </tbody>
</table>

## MakeDb Class
<table>
  <thead>
    <tr>
      <th>Method</th>
      <th>Parameters</th>
      <th>Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><code>__init__(adata, db_name, db_path, create_all_indexes)</code></td>
      <td>
          <li><code>adata</code>: AnnData object (required)</li>
          <li><code>db_name</code>: Name for the database (required)</li>
          <li><code>db_path</code>: Path to store the database (default: 'db/')</li>
          <li><code>create_all_indexes</code>: Boolean (optional. default: False). <i>Warning: Runtime can be significant when building.</i></li>
      </td>
      <td>Initializes the MakeDb object and validates parameters, then proceeds to build the DuckDB database.</td>
    </tr>
  </tbody>
</table>


<br>
<br>

## Notes
There are two key reasons to use **AnnSQL**: (1) if you prefer SQL's expressive syntax for filtering and querying your data, or (2) if you're working with datasets that exceed memory limits and require loading AnnData in backed mode. Using backed mode in AnnData can limit available functions, especially aggregate operations, and slow down data access. **AnnSQL** offers a solution by enabling SQL-style queries that may perform more efficiently in these scenarios. Below are rough runtime comparisons between AnnData (backed mode) and AnnSQL once the database has been built. 


<img src="examples/images/comparision_backed.png">

<br>
<br>

## Citation
Coming soon...