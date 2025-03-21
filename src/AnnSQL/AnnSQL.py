from .BuildDb import BuildDb
import scanpy as sc
import pandas as pd
import numpy as np
from numpy.linalg import eig
import duckdb
import warnings
import logging
import time
import os
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors
from scipy import stats
from scipy.stats import t as tdist
import umap

class AnnSQL:
	def __init__(self, adata=None, db=None, create_all_indexes=False, create_basic_indexes=False, print_output=True, layers=["X", "obs", "var", "var_names", "obsm", "varm", "obsp", "uns"], memory_limit=None, db_config={}):
		"""
		Initializes an instance of the AnnSQL class. This class is used to query and update a database created from an AnnData object. 
		it also provides methods for preprocessing and basic analysis. The in-process database engine is DuckDB and the database is 
		stored in memory by default, However, the database can be loaded from a file path by providing the db parameter. Databases can be
		built from an AnnData object by using the MakeDb class.

		Parameters:
			adata (AnnData or None): An AnnData object containing the data to be stored in the database. If None, an empty AnnData object will be created.
			db (str or None): The path to an existing database file. 
			create_all_indexes (bool): Whether to create indexes for all columns in the database. Memory intensive. Default is False.
			create_basic_indexes (bool): Whether to create indexes for basic columns. Memory intensive. Default is False.
			print_output (bool): Whether to print output messages for database creation. Default is True.
			layers (list): A list of layer names to be stored in the database. Default is ["X", "obs", "var", "var_names", "obsm", "varm", "obsp", "uns"].
			memory_limit (str): The memory limit for the DuckDB database. Default is None.

		Returns:
			AnnSQL (Object): An instance of the AnnSQL class.
		"""		
		self.adata = self.open_anndata(adata)
		self.db = db
		self.create_basic_indexes = create_basic_indexes
		self.create_all_indexes = create_all_indexes
		self.layers = layers
		self.validate_params()
		self.is_open = False
		self.print_output = print_output
		self.memory_limit = memory_limit
		self.db_config = db_config
		if self.db is None:
			self.build_db()
		else:
			self.open_db()

	def validate_params(self):
		"""
		Validates the parameters of the AnnSQL object.

		Raises:
			ValueError: If both `adata` and `db` parameters are not defined.
			ValueError: If `adata` is defined but not an instance of `scanpy.AnnData`.
			ValueError: If `db` is defined but the file does not exist.
		"""

		if self.adata is None and self.db is None:
			raise ValueError('Both adata and db parameters not defined. Select an option')
		if self.adata is not None:
			if not isinstance(self.adata, sc.AnnData):
				raise ValueError('adata must be a scanpy AnnData object')
		if self.db is not None:
			if not os.path.exists(self.db):
				raise ValueError('The db provided doesn\'t exist. Please check the path')

	def open_anndata(self,adata):
		"""
		Opens an AnnData object.

		Parameters:
			adata (AnnData or str): The AnnData object to be opened or the file path to the AnnData object.

		Returns:
			adata (AnnData): The opened AnnData object instance.
		"""

		if not isinstance(adata, sc.AnnData) and isinstance(adata, str):	
			return sc.read_h5ad(adata)
		else:
			return adata
			
	def open_db(self):
		"""
		Opens a connection to the database specified by `self.db`.
		If a connection is already open, it will be closed before opening a new one.

		Parameters:
			None
		Returns:
			None
		"""

		if self.db is not None:
			self.conn = duckdb.connect(self.db, config=self.db_config)
			self.is_open = True

	def close_db(self):
		"""
		Closes the database connection.
		If the database connection is open, this method closes the connection and sets the `is_open` attribute to False.

		Parameters:
			None
		Returns:
			None
		"""

		if self.db is not None:
			self.conn.close()
			self.is_open = False

	def asql_register(self, table_name, df):
		"""
		Registers a table in the database from a  pandas df.

		Parameters:
			table_name (str): The name of the table to be registered.
			df (pandas.DataFrame): The DataFrame containing the data to be registered.
		Returns:
			None
		"""

		self.open_db()
		self.conn.register(table_name, df)
		self.close_db()

	def build_db(self):
		"""
		Builds the database connection and initializes the necessary tables and indexes.
		
		Parameters:
			None
		Returns:
			None
		"""
		self.conn = duckdb.connect(':memory:', config=self.db_config)
		db = BuildDb(adata=self.adata, conn=self.conn, create_all_indexes=self.create_all_indexes, create_basic_indexes=self.create_basic_indexes, layers=self.layers, print_output=self.print_output)
		self.conn = db.conn

	def query(self, query, return_type='pandas'):
		"""
		Executes the given SQL query and returns the result based on the specified return type.

		Parameters:
			query (str): The SQL query to be executed.
			return_type (str, optional): The desired return type of the query result. Options are 'pandas', 'adata', and 'parquet'. Defaults to 'pandas'.

		Returns:
			results (pandas df, AnnData, or parquet): The result of the query based on the specified return type.

		Raises:
			ValueError: If the return_type is not one of 'pandas', 'adata', or 'parquet'.
			ValueError: If the query contains 'UPDATE', 'DELETE', or 'INSERT' statements. Use update_query() instead for such statements.
		
		Examples:
			>>> # Query the X layer and return the result as a pandas DataFrame
			>>> asql.query("SELECT * FROM X LIMIT 5")
			>>> # Query the X layer and return the result as an AnnData object
			>>> asql.query("SELECT * FROM X LIMIT 5", return_type='parquet')
			
		"""

		if return_type not in ['pandas', 'adata', 'parquet']:
			raise ValueError('return_type must be either pandas, parquet or adata')
		if 'UPDATE' in query.upper() or 'DELETE' in query.upper() or 'INSERT' in query.upper():
			raise ValueError('UPDATE, DELETE, and INSERT detected. Please use update_query() instead')

		self.open_db()

		if self.memory_limit is not None:			
			self.conn.execute(f"SET memory_limit = '{self.memory_limit}';")

		if return_type == 'parquet' and 'SELECT' in query.upper():
			query = "COPY ("+query+") TO 'output.parquet' (FORMAT PARQUET);"
			self.conn.execute(query)
			logging.info("Query results saved as 'query.parquet' file in the current directory")
		else:
			result_df = self.conn.execute(query).df()
		self.close_db()

		if return_type == 'pandas':
			return result_df
		elif return_type == 'adata':
			if self.db is not None and self.adata is None:
				print('Warning: No adata object provided. return_type="adata" is disabled.')
				return result_df
			return self.adata[result_df["cell_id"]]

	def query_raw(self, query):
		"""
		Executes a raw SQL query on the database.

		Parameters:
			query (str): The SQL query to be executed.

		Returns:
			result (DuckDb Object): The result of the query execution.

		Examples:
			>>> asql.query_raw("SELECT * FROM X LIMIT 5")
		"""

		self.open_db()
		if self.memory_limit is not None:			
			self.conn.execute(f"SET memory_limit = '{self.memory_limit}';")
		result = self.conn.execute(query)
		self.close_db()
		return result

	def update_query(self, query, suppress_message=False):
		"""
		Executes an update query on the database.

		Parameters:
			query (str): The SQL query to be executed.
			suppress_message (bool, optional): Whether to suppress the success message. Defaults to False.

		Raises:
			ValueError: If the query contains 'SELECT' or 'DELETE' statements.

		Returns:
			result (DuckDb Object): The result of the query execution.
		
		Examples:
			>>> asql.update_query("UPDATE obs SET cell_type = 'Dendritic' WHERE leiden_cluster = 0")
		"""


		if 'SELECT' in query.upper() or 'DELETE' in query.upper():
			raise ValueError('SELECT detected. Please use query() instead')
		try:
			self.open_db()
			self.conn.execute(query)
			self.close_db()
			if suppress_message == False:
				print("Query Successful")
		except Exception as e:
			print("Update Query Error:", e)

	def delete_query(self, query, suppress_message=False):
		"""
		Executes a delete query on the database.

		Parameters:
			query (str): The delete query to be executed.
			suppress_message (bool, optional): Whether to suppress the success message. Defaults to False.

		Raises:
			ValueError: If the query contains 'SELECT' keyword.

		Returns:
			result (DuckDb Object): The result of the delete.

		Examples:
			>>> asql.delete_query("DELETE FROM X WHERE cell_id IN (SELECT cell_id FROM obs WHERE cell_type = 'Dendritic')")	
		"""

		if 'DELETE' not in query.upper():
			raise ValueError('SELECT detected. Please use query() instead')
		try:
			self.open_db()
			self.conn.execute(query)
			self.close_db()
			if suppress_message == False:
				print("Delete Query Successful")
		except Exception as e:
			print("Delete Query Error:", e)

	def show_tables(self):
		"""
		A simple helper method to retrieve a list of table names from the 'main' schema in the database.

		Returns:
			result (pandas.DataFrame): A DataFrame containing the table names.

		Examples:
			>>> asql.show_tables()
		"""

		self.open_db()
		result = self.conn.execute("SELECT table_name FROM information_schema.tables  WHERE table_schema='main'").df()
		self.close_db()
		return result

	def show_settings(self):
		"""
		A simple helper method to retrieve a list of the duckdb database settings

		Returns:
			result (pandas.DataFrame): A DataFrame containing the configuration options.

		Examples:
			>>> asql.show_settings()
		"""

		self.open_db()
		result = self.conn.execute("SELECT * FROM duckdb_settings()").df()
		self.close_db()
		return result

	def export_parquet(self):
		"""
		This method exports all tables as parquet files.
		This method retrieves the list of tables using the `show_tables` method and exports each table as a parquet file.
		The parquet files are saved in the 'parquet_files' directory.
		"""

		tables = self.show_tables()
		if not os.path.exists("parquet_files"):
			os.mkdir("parquet_files")
		for table in tables["table_name"]:
			query = "SELECT * FROM "+table
			query = "COPY ("+query+") TO 'parquet_files/"+table+".parquet' (FORMAT PARQUET);"
			self.open_db()
			self.conn.execute(query)
			self.close_db()
		logging.info("All tables exported as parquet files in the 'parquet_files' directory")
	
	def replace_special_chars(self, string):
		"""
		Replaces special characters in a string with underscores. Is useful when creating tables from AnnData objects as certain characters in gene names can cause issues with column names.

		Parameters:
			String (str): The input string.

		Returns:
			String (str): The modified string with special characters replaced by underscores.
		"""

		return string.replace("-", "_").replace(".", "_")

	def expression_normalize(self, total_counts_per_cell=10000, chunk_size=200, print_progress=False):
		"""
		Normalize the expression values in the dataset to a desired total counts per cell.

		Parameters:
			total_counts_per_cell (int, optional): The desired total counts per cell after normalization. Defaults to 10000.
			chunk_size (int, optional): The number of genes to process in each chunk. Defaults to 200.
			print_progress (bool, optional): Whether to print progress information. Defaults to False.

		Notes:
			- This method normalizes the expression values in the dataset by dividing each gene's expression value by the total counts and then multiplying it by the desired total counts per cell.
			- The normalization is performed in chunks to decrease memory usage.
			- If 'total_counts' column is not found in the dataset, it will be calculated using the calculate_total_counts method.

		Example:
			asql = AnnSQL(db='db/pbmc.asql')			
			asql.expression_normalize(total_counts_per_cell=10000, chunk_size=100, print_progress=True)
		"""

		self.check_chunk_size(chunk_size)

		#check if total_counts column exists
		if 'total_counts' not in self.query(f"Describe obs")['column_name'].values:
			self.calculate_total_counts(chunk_size=chunk_size,print_progress=print_progress)

		print("Expression Normalization Started")
		gene_names = self.query(f"Describe X")['column_name'][1:].values
		if 'total_counts' in gene_names:
			gene_names = gene_names[:-1]
		for i in range(0, len(gene_names), chunk_size):
			updates = []
			chunk = gene_names[i:i + chunk_size]
			for gene in chunk:
				if gene == 'total_counts':
					continue
				updates.append(f"{gene} = (({gene} / total_counts) * {total_counts_per_cell})")
			update_query = f"UPDATE X SET {', '.join(updates)}"
			self.update_query(update_query, suppress_message=True)
			if print_progress == True:
				print(f"Processed chunk {i // chunk_size + 1}")
		print("Expression Normalization Complete")


	def expression_log(self, log_type="LN", chunk_size=200, print_progress=False):
		"""
		Log transform the expression values of genes in the dataset.

		Parameters:
			log_type (str, optional): The type of logarithm to use for the transformation. 
				Possible values are "LN" (natural logarithm), "LOG" (base 10 logarithm), 
				"LOG2" (base 2 logarithm), and "LOG10" (base 10 logarithm). 
				Defaults to "LN".
			chunk_size (int, optional): The number of genes to process in each chunk. 
				Defaults to 200.
			print_progress (bool, optional): Whether to print progress information during the transformation. 
				Defaults to False.

		Notes:
			- This method log-transforms the expression values of genes in the dataset.
			- The log transformation is performed in chunks to decrease memory usage.
			- The log transformation is applied to each gene individually.
		
		Example:
			asql = AnnSQL(db='db/pbmc.asql')			
			asql.expression_log(log_type="LN", chunk_size=100)
		"""



		#log_type can be LN, LOG (LOG2 alias), LOG2, LOG10
		self.check_chunk_size(chunk_size)
		gene_names = self.query(f"Describe X")['column_name'][1:].values
		if 'total_counts' in gene_names:
			gene_names = gene_names[:-1]
		
		print("Log Transform Started")
		for i in range(0, len(gene_names), chunk_size):
			updates = []
			chunk = gene_names[i:i + chunk_size]
			for gene in chunk:
				if gene == 'total_counts':
					continue
				updates.append(f"{gene} = {log_type}({gene}+1)") #handle zero values like scanpy
			update_query = f"UPDATE X SET {', '.join(updates)}"
			self.update_query(update_query, suppress_message=True)
			if print_progress == True:
				print(f"Processed chunk {i // chunk_size + 1}")
		print("Log Transform Complete")

	
	def calculate_total_counts(self, chunk_size=200, print_progress=False):
		"""
		Calculate the total counts for each gene in the X table and update the total_counts column.
		Also update the total_counts column in the obs table based on the corresponding values in the X table.

		Parameters:
			chunk_size (int): The number of gene names to process in each chunk. Default is 200.
			print_progress (bool): Whether to print progress information during the calculation. Default is False.

		Notes:
			- This method calculates the total counts for each gene in the X table by summing the expression values of all cells for each gene.
			- The method also updates the 'total_counts' column in the obs table based on the corresponding values in the X table.

		Example:
			asql = AnnSQL(db='db/pbmc.asql')			
			asql.calculate_total_counts(chunk_size=100, print_progress=True)
		"""

		self.check_chunk_size(chunk_size)
		gene_names = self.query(f"Describe X")['column_name'][1:].values
		
		if "total_counts" in gene_names:
			self.update_query(f"UPDATE X SET total_counts = 0;")
			gene_names = gene_names[:-1] 
		else:
			self.query(f"ALTER TABLE X ADD COLUMN total_counts FLOAT DEFAULT 0;")
		
		print("Total Counts Calculation Started")
		for i in range(0, len(gene_names), chunk_size):
			chunk = gene_names[i:i+chunk_size]
			chunk = " + ".join(chunk) + " + total_counts"
			self.update_query(f"UPDATE X SET total_counts = ({chunk});", suppress_message=True)
			if print_progress == True:
				print(f"Processed chunk {i // chunk_size + 1}")

		#set obs total_counts
		if 'total_counts' not in self.query("SELECT * FROM obs LIMIT 1").columns:
			self.query_raw("ALTER TABLE obs ADD COLUMN total_counts FLOAT DEFAULT 0;")
		self.query_raw("UPDATE obs SET total_counts = (SELECT total_counts FROM X WHERE obs.cell_id = X.cell_id)")
		print("Total Counts Calculation Complete")

	def calculate_gene_counts(self, chunk_size=200, print_progress=False, gene_field="gene_names"):
		"""
		Calculate gene counts and gene means for each gene in the var table.

		Parameters:
			chunk_size (int): The number of genes to process in each chunk. Default is 200.
			print_progress (bool): Whether to print progress information. Default is False.
			gene_field (str): The name of the column in the var table that contains gene names. Default is "gene_names".

		Notes:
			- This method calculates the gene counts and gene means for each gene in the var table by summing the expression values of all cells for each gene.
			- The method updates the 'gene_counts' and 'gene_mean' columns in the var table based on the calculated values.
		
		Example:
			asql = AnnSQL(db='db/pbmc.asql')			
			asql.calculate_gene_counts(chunk_size=100)
		"""
		
		self.check_chunk_size(chunk_size)
		gene_names_df = self.query(f"SELECT {gene_field} FROM var")
		gene_names_df["gene_counts"] = 0.0	
		gene_names_df["gene_mean"] = 0.0
		gene_names_df = gene_names_df.reset_index(drop=True)
		
		var_table = self.query("SELECT * FROM var LIMIT 1")
		if var_table.shape[0] == 0:
			print("Creating Var Table")
			self.open_db()
			self.conn.register("gene_names_df", gene_names_df)
			self.conn.execute(f"CREATE TABLE var AS SELECT * FROM gene_names_df")
		else:
			print("Updating Var Table")
			if "gene_counts" not in var_table.columns or "gene_mean" not in var_table.columns:
				#drop the columns if they exist
				self.update_query("ALTER TABLE var DROP COLUMN IF EXISTS gene_counts;", suppress_message=True)
				self.update_query("ALTER TABLE var DROP COLUMN IF EXISTS gene_mean;", suppress_message=True)
				self.update_query("ALTER TABLE var ADD COLUMN gene_counts FLOAT DEFAULT 0;", suppress_message=True)
				self.update_query("ALTER TABLE var ADD COLUMN gene_mean FLOAT DEFAULT 0;", suppress_message=True)
			else:
				self.update_query("UPDATE var SET gene_counts = 0.0;", suppress_message=True)
				self.update_query("UPDATE var SET gene_mean = 0.0;", suppress_message=True)

		print("Gene Counts Calculation Started")
		gene_counts = []
		gene_means = []
		for i in range(0, len(gene_names_df), chunk_size):
			chunk = gene_names_df[gene_field][i:i+chunk_size]
			query = f"SELECT {', '.join([f'SUM({gene}) as {gene}' for gene in chunk])} FROM X;"
			query2 = f"SELECT {', '.join([f'AVG({gene}) as {gene}' for gene in chunk])} FROM X;"
			counts_chunk = self.query(query)
			counts_chunk2 = self.query(query2)
			gene_counts.extend(counts_chunk.values.flatten())
			gene_means.extend(counts_chunk2.values.flatten())
			if print_progress == True:
				print(f"Processed chunk {i // chunk_size + 1}")

		#insert these values into the var table matching on the index.
		gene_counts_df = pd.DataFrame({"gene_counts": gene_counts,"gene_mean": gene_means})
		gene_counts_df[gene_field] = gene_names_df[gene_field]

		#update the var table with the gene_counts values
		self.open_db()
		self.conn.execute("DROP TABLE IF EXISTS gene_counts_df")
		self.conn.register("gene_counts_df", gene_counts_df)
		self.conn.execute(f"CREATE TABLE gene_counts_df AS SELECT * FROM gene_counts_df")
		
		#self.conn.execute(f"UPDATE var SET gene_counts = (SELECT gene_counts FROM gene_counts_df WHERE var.{gene_field} = gene_counts_df.{gene_field})")
		self.conn.execute(f"""
			UPDATE var 
			SET gene_counts = (SELECT gene_counts FROM gene_counts_df WHERE var.{gene_field} = gene_counts_df.{gene_field}),
				gene_mean = (SELECT gene_mean FROM gene_counts_df WHERE var.{gene_field} = gene_counts_df.{gene_field})
		""")

		self.conn.execute("DROP VIEW IF EXISTS gene_counts_df")
		self.query("DROP TABLE IF EXISTS gene_counts_df CASCADE")
		print("Gene Counts Calculation Complete")


	def calculate_variable_genes(self, chunk_size=100, print_progress=False, gene_field="gene_names", save_var_names=True,save_top_variable_genes=2000):
		"""
		Calculates variable genes based on the given parameters. This method uses the duckdb VARIANCE function to calculate the variance of each gene 
		in the X table and updates the 'variance' column in the var table.

		Parameters:
			chunk_size (int, optional): The size of each chunk for processing. Defaults to 100.
			print_progress (bool, optional): Whether to print progress while processing. Defaults to False.
			gene_field (str, optional): The field name for gene names in the database table. Defaults to "gene_names".
			save_top_variable_genes (int, optional): The number of top variable genes to save. Defaults to 2000.
			save_var_names (bool, optional): Whether to save the variable gene names to the X table and remove all others. Defaults to True.

		Notes:
			- This method calculates the variance of each gene in the X table using the VARIANCE function.
			- The method updates the 'variance' column in the var table based on the calculated values.
			- If 'variance' column is not found in the var table, it will be created.
			- If 'variance' column already exists in the var table, it will be updated.
			- If 'variance' column is not found in the var table, it will be created.

		Example:
			asql = AnnSQL(db='db/pbmc.asql')			
			asql.calculate_variable_genes(chunk_size=250)
		"""

		
		self.check_chunk_size(chunk_size)
		gene_names_df = self.query(f"SELECT {gene_field} FROM var")
		gene_names_df["variance"] = 0.0		
		gene_names_df = gene_names_df.reset_index(drop=True)
		
		var_table = self.query("SELECT * FROM var LIMIT 1")
		if var_table.shape[0] == 0:
			print("Creating Var Table")
			self.open_db()
			self.conn.register("gene_names_df", gene_names_df)
			self.conn.execute(f"CREATE TABLE var AS SELECT * FROM gene_names_df")
			self.update_query("ALTER TABLE var ADD COLUMN variance FLOAT DEFAULT 0;", suppress_message=True)
		else:
			print("Updating Var Table")
			if "variance" not in var_table.columns:
				self.update_query("ALTER TABLE var ADD COLUMN variance FLOAT DEFAULT 0;", suppress_message=True)
			else:
				self.update_query("UPDATE var SET variance = 0.0;", suppress_message=True)	

		variance_values = []
		for i in range(0, len(gene_names_df), chunk_size):
			chunk = gene_names_df[gene_field][i:i+chunk_size]
			query = f"SELECT {', '.join([f'VARIANCE({gene}) as {gene}' for gene in chunk])} FROM X;"
			variance_chunk = self.query(query)
			variance_values.extend(variance_chunk.values.flatten())
			if print_progress == True:
				print(f"Processed chunk {i // chunk_size + 1}")

		#insert these values into the var table matching on the index.
		variance_df = pd.DataFrame({"variance": variance_values})
		variance_df[gene_field] = gene_names_df[gene_field]

		#update the var table with the variance values
		self.open_db()
		self.conn.execute("DROP TABLE IF EXISTS variance_df")
		self.conn.register("variance_df", variance_df)
		self.conn.execute(f"CREATE TABLE variance_df AS SELECT * FROM variance_df")
		self.conn.execute(f"UPDATE var SET variance = (SELECT variance FROM variance_df WHERE var.{gene_field} = variance_df.{gene_field})")
		self.conn.execute("DROP VIEW IF EXISTS variance_df")

		if save_var_names == True:
			self.save_highly_variable_genes(top_variable_genes=save_top_variable_genes, gene_field=gene_field)

		print("Variance Calculation Complete")


	def build_meta_cells(self, primary_cluster=None, secondary_cluster=None, aggregate_type="AVG", table_name="meta_cells", chunk_size=100, print_progress=False):
		"""
		Builds a meta_cells table by aggregating data from the X and obs tables. This will group cells by the primary and secondary cluster columns and aggregate 
		the expression values using the specified aggregate type.

		Parameters:
			primary_cluster (str, optional): The name of the primary cluster column. Defaults to None.
			secondary_cluster (str, optional): The name of the secondary cluster column. Defaults to None.
			aggregate_type (str, optional): The type of aggregation to perform. Defaults to "AVG".
			table_name (str, optional): The name of the table to create. Defaults to "meta_cells".
			chunk_size (int, optional): The size of each processing chunk. Defaults to 100.
			print_progress (bool, optional): Whether to print progress information. Defaults to False.

		Notes:
			- This method creates a new table with the specified name and aggregates the data from the X and obs tables based on the primary and secondary cluster columns.
			- The method aggregates the expression values using the specified aggregate type (e.g., AVG, SUM, etc.).
			- The method processes the data in chunks to decrease memory usage.
		
		The example below builds a table of average expression values of each gene for each cell type.

		Example:
			asql = AnnSQL(db='db/pbmc.asql')
			asql.build_meta_cells(primary_cluster="cell_type", aggregate_type="AVG", table_name="meta_cells", chunk_size=250)		
		"""

		self.check_chunk_size(chunk_size)
		columns = self.query("DESCRIBE X")[1:]["column_name"].tolist()
		self.query_raw(f"DROP TABLE IF EXISTS {table_name}")
		self.query_raw(f"""
						CREATE TABLE {table_name} AS 
							SELECT CAST('' as Varchar) as {primary_cluster},
							CAST('' as Varchar) as {secondary_cluster},
							CAST(0 as Float) as cell_count,
							* 
						FROM X WHERE FALSE;""")

		if primary_cluster is not None and secondary_cluster is not None:
			self.query_raw(f"""
				INSERT INTO {table_name} ({primary_cluster}, {secondary_cluster}, cell_count)
				SELECT 
					obs.{primary_cluster} as {primary_cluster},
					obs.{secondary_cluster} as {secondary_cluster},
					0 as cell_count
				FROM obs
				GROUP BY obs.{primary_cluster}, obs.{secondary_cluster}
			""")

			#process in chunks
			for i in range(0, len(columns), chunk_size):
				if print_progress == True:
					print(f"Processing chunk {i + 1} of {i + chunk_size}")
				chunk_columns = columns[i:i+chunk_size]
				chunk_query = ", ".join([f"{aggregate_type}(X.{col}) as {col}" for col in chunk_columns])
				update_query = f"""
					UPDATE {table_name}
					SET 
						cell_count = sub.cell_count,
						{", ".join([f"{col} = sub.{col}" for col in chunk_columns])}
					FROM (
						SELECT 
							obs.{primary_cluster} as {primary_cluster},
							obs.{secondary_cluster} as {secondary_cluster},
							COUNT(X.cell_id) as cell_count,
							{chunk_query}
						FROM X
						INNER JOIN obs ON X.cell_id = obs.cell_id
						GROUP BY obs.{secondary_cluster}, obs.{primary_cluster}
					) as sub
					WHERE 
						{table_name}.{primary_cluster} = sub.{primary_cluster} AND
						{table_name}.{secondary_cluster} = sub.{secondary_cluster}
				"""
				self.query_raw(update_query)
		else:
			self.query_raw(f"""
				INSERT INTO {table_name} ({primary_cluster}, cell_count)
				SELECT 
					obs.{primary_cluster} as {primary_cluster},
					0 as cell_count
				FROM obs
				GROUP BY obs.{primary_cluster}
			""")

			#process in chunks
			for i in range(0, len(columns), chunk_size):
				if print_progress == True:
					print(f"Processing chunk {i + 1} of {i + chunk_size}")
				chunk_columns = columns[i:i+chunk_size]
				chunk_query = ", ".join([f"{aggregate_type}(X.{col}) as {col}" for col in chunk_columns])
				update_query = f"""
					UPDATE {table_name}
					SET 
						cell_count = sub.cell_count,
						{", ".join([f"{col} = sub.{col}" for col in chunk_columns])}
					FROM (
						SELECT 
							obs.{primary_cluster} as {primary_cluster},
							COUNT(X.cell_id) as cell_count,
							{chunk_query}
						FROM X
						INNER JOIN obs ON X.cell_id = obs.cell_id
						GROUP BY obs.{primary_cluster}
					) as sub
					WHERE 
						{table_name}.{primary_cluster} = sub.{primary_cluster}
				"""
				self.query_raw(update_query)
		
		print(f"{table_name} table created. You may now query the table for results.")


	def check_chunk_size(self, chunk_size):
		"""
		Check if the given chunk size is valid. DuckDb imposes a limit of 1000 (999) on the operations that can be performed in a single query. 
		it can be exceeded in some cases, but it's not recommended and uses more memory. We've implemented a check to ensure the chunk size is 
		within the limit.

		Parameters:
			chunk_size (int): The chunk size to be checked.

		Raises:
			ValueError: If the chunk size is greater than 999.
		"""

		if chunk_size > 999:
			raise ValueError('chunk_size must be less than 1000. DuckDb limitation')

	
	def filter_by_cell_counts(self, min_cell_count=None, max_cell_count=None):
		"""
		Filter cells based on their total counts exisiting in the obs table in the total_counts column.

		Parameters:
			min_cell_count (int, optional): Minimum total count threshold. Cells with total counts less than this value will be removed. Defaults to None.
			max_cell_count (int, optional): Maximum total count threshold. Cells with total counts greater than this value will be removed. Defaults to None.

		Notes:
			- This method removes cells from the X and obs tables based on the specified total count thresholds.
			- If 'total_counts' column is not found in the obs table, it will be calculated using the calculate_total_counts method.
		
		Example:
			asql = AnnSQL(db='db/pbmc.asql')			
			asql.filter_by_cell_counts(min_cell_count=1000, max_cell_count=45000)
		"""

		if 'total_counts' not in self.query("SELECT * FROM obs LIMIT 1").columns:
			print("Total counts not found. Running total counts...")
			self.calculate_total_counts()
		if min_cell_count >= 0 and max_cell_count is None:
			query_x = f"DELETE FROM X WHERE cell_id IN (SELECT cell_id FROM obs WHERE total_counts < {min_cell_count})"
			query_obs = f"DELETE FROM obs WHERE total_counts < {min_cell_count}"
			self.delete_query(query_x, suppress_message=True)
			self.delete_query(query_obs, suppress_message=True)
			print(f"Cells with total counts less than {min_cell_count} removed")
		elif min_cell_count is None and max_cell_count >= 0:
			query_x = f"DELETE FROM X WHERE cell_id IN (SELECT cell_id FROM obs WHERE total_counts > {max_cell_count})"
			query_obs = f"DELETE FROM obs WHERE total_counts > {max_cell_count}"
			self.delete_query(query_x, suppress_message=True)
			self.delete_query(query_obs, suppress_message=True)
			print(f"Cells with total counts greater than {max_cell_count} removed")
		elif min_cell_count >= 0 and max_cell_count >= 0:
			query_x = f"DELETE FROM X WHERE cell_id IN (SELECT cell_id FROM obs WHERE total_counts < {min_cell_count} OR total_counts > {max_cell_count})"
			query_obs = f"DELETE FROM obs WHERE total_counts < {min_cell_count} OR total_counts > {max_cell_count}"
			self.delete_query(query_x, suppress_message=True)
			self.delete_query(query_obs, suppress_message=True)
			print(f"Cells with total counts less than {min_cell_count} and greater than {max_cell_count} removed")


	def calculate_pca(self, n_pcs=50, 
						table_name="X", 
						chunk_size=100, 
						print_progress=False, 
						zero_center=False, 
						top_variable_genes=2000,
						max_cells_memory_threshold=25000,
						gene_field="gene_names"):

		"""
		Calculates the Principal Component Analysis (PCA) for the data stored in the specified table. This method uses the top variable genes 
		based on variance to perform PCA. The PCA calculation is performed as a hybrid of SQL and python. 

		Parameters:
			n_pcs (int, default=50): Number of principal components to calculate.
			table_name (str, default="X"): Name of the table containing the data.
			chunk_size (int, default=100): Size of the data chunks to process.
			print_progress (bool, default=False): Whether to print progress messages.
			zero_center (bool, default=False): Whether to zero-center the data before PCA.
			top_variable_genes (int, default=2000): Number of top variable genes to use for PCA.
			max_cells_memory_threshold (int, default=5000): Maximum number of cells to hold in memory before using SQL. Beyond this threshold, covariance matrix calculation is done in SQL to be more memory efficient.
			gene_field (str, default="gene_names"): The field name for gene names in the var database table.

		Functionality:
			1. Checks if the specified table exists. If not, raises a ValueError.
			2. Checks if the 'variance' column exists in the var table. If not, calls calculate_variable_genes to compute it.
			3. Retrieves the top variable genes based on variance.
			4. Constructs a query to standardize the data, with an option for zero-centering.
			5. Creates a wide standardized table.
			6. If the number of cells is greater than the specified threshold, calculates the covariance matrix using SQL.
			7. Calculates the eigenvalues and eigenvectors using numpy as this is a small matrix and does not require SQL.
			8. Creates a table for the eigenvalues.
			9. Creates a table for the eigenvectors.
			10. Creates a table for the principal components.
		
		Example:
			asql = AnnSQL(db='db/pbmc.asql')			
			asql.calculate_pca(n_pcs=50, table_name="X", chunk_size=250, zero_center=False, top_variable_genes=2000, max_cells_memory_threshold=5000)
		"""


		#does the table exist?
		if table_name not in self.show_tables()['table_name'].tolist():
			raise ValueError(f"{table_name} table not found.")
		
		#does the variance column exist?
		if 'variance' not in self.query("SELECT * FROM var LIMIT 1").columns:
			print("Variance not found. Running calculate_variable_genes...")
			self.calculate_variable_genes(chunk_size=chunk_size, print_progress=print_progress, 
										save_var_names=True, save_top_variable_genes=top_variable_genes)
		
		print("PCA Calculation Started (4 Steps)")

		#get the top genes to use
		genes_df = self.query(f"SELECT {gene_field} FROM var ORDER BY variance DESC")
		genes = genes_df[gene_field].tolist()[:top_variable_genes]
		
		#build query for wide standardized table with two options for zero-centering
		print("Step 1: Building Wide Standardized Table")
		col_exprs = []
		for gene in genes:
			if zero_center:
				expr = f"(({gene} - AVG({gene}) OVER ()) / STDDEV({gene}) OVER ()) AS {gene}"
			else:
				expr = f"({gene} - AVG({gene}) OVER ()) AS {gene}"
			col_exprs.append(expr)
		cols_sql = ", ".join(col_exprs)
		
		#insert the wide standardized table
		self.query_raw("DROP TABLE IF EXISTS X_standard_wide;")
		self.query_raw(f"CREATE TABLE X_standard_wide AS SELECT cell_id, {cols_sql} FROM {table_name};")

		print("Step 2: Calculating Covariance Matrix")

		#if there's less than n cells in X_standard_wide use np.cov method. 
		#the SQL method isn't as fast, but is more memory efficient.
		if self.query("SELECT COUNT(*) as total FROM X_standard_wide")["total"][0] < max_cells_memory_threshold:
		
			#get the gene data
			wide_df = self.query("SELECT * FROM X_standard_wide ORDER BY cell_id")
			
			#column order matters for covariance matrix
			wide_df = wide_df[['cell_id'] + genes]
			
			#gene data as a numpy array.
			gene_data = wide_df[genes].to_numpy()

			#covariance matrix for the genes.
			cov_matrix_np = np.cov(gene_data, rowvar=False)


		else:
			
			covariance_matrix_df = pd.DataFrame(index=genes, columns=genes, dtype=float)

			inc=1
			for i, gene_1 in enumerate(genes): 
				start_time = time.time()
				gene_2_list = genes[i:]  #upper triangle
				for j in range(0, len(gene_2_list), chunk_size):
					gene_2_chunk = gene_2_list[j:j + chunk_size] 
					genes_clause = ", ".join([f"covar_samp({gene_1}, {gene_2}) AS cov_{gene_1}_{gene_2}" for gene_2 in gene_2_chunk])
					result = self.query(f"SELECT {genes_clause} FROM X_standard_wide")
					for k, gene_2 in enumerate(gene_2_chunk):
						covariance_matrix_df.at[gene_1, gene_2] = result.iloc[0, k]
				if print_progress == True:
					end_time = time.time() - start_time
					print(f"Covariance: {inc} genes out of {len(genes)}: {str(end_time)} seconds")
				inc += 1

			covariance_matrix_df = covariance_matrix_df.fillna(0)	

			#convert to np (small matrix, okay to represent as numpy)
			cov_matrix_np = covariance_matrix_df.to_numpy()
			
			self.query_raw("DROP TABLE IF EXISTS X_covariance;")
			self.open_db()
			self.conn.register("covariance_matrix_df", covariance_matrix_df)
			self.conn.execute("CREATE TABLE X_covariance AS SELECT * FROM covariance_matrix_df")
			self.close_db()


		print("Step 3: Calculating Eigenvalues and Eigenvectors")

		#use linalg.eigh for the eigenvalues and eigenvectors (cov matrix is small)
		eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix_np)
		sorted_idx = np.argsort(eigenvalues)[::-1]
		eigenvalues_sorted = eigenvalues[sorted_idx]
		eigenvectors_sorted = eigenvectors[:, sorted_idx][:, :n_pcs]  # take only top n_pcs
		
		#create a table for the eigenvalues
		eigenvalues_df = pd.DataFrame({
			"pc": np.arange(n_pcs),
			"eigenvalue": eigenvalues_sorted[:n_pcs]
		})

		self.query_raw("DROP TABLE IF EXISTS PC_eigenvalues;")
		self.open_db() #The way is shut. It was made by those who are...
		self.conn.execute("CREATE TABLE PC_eigenvalues AS SELECT * FROM eigenvalues_df;")
		self.close_db()

		#df for PC loadings
		pc_loadings_df = pd.DataFrame({
			"gene": np.repeat(genes, n_pcs),
			"pc": np.tile(np.arange(n_pcs), len(genes)),
			"loading": eigenvectors_sorted.flatten()
		})
		
		#insert the loadings
		self.query_raw("DROP TABLE IF EXISTS PC_loadings;")
		self.open_db()
		self.conn.execute("CREATE TABLE PC_loadings AS SELECT * FROM pc_loadings_df;")
		self.close_db()
		

		print("Step 4: Calculating PC Scores")

		#if there's less than 10k cells in X_standard_wide use np.dot method. SQL method isn't as fast, but is more memory efficient.
		if self.query("SELECT COUNT(*) as total FROM X_standard_wide")["total"][0] < max_cells_memory_threshold:

			#dot product of each cell's standardized gene vector with the loadings.
			pc_scores = gene_data.dot(eigenvectors_sorted) 
			
			#long-form df for PC scores
			cell_ids = wide_df['cell_id'].tolist()
			pc_scores_list = []
			n_cells = len(cell_ids)
			for i in range(n_cells):
				for pc in range(n_pcs):
					pc_scores_list.append({
						"cell_id": cell_ids[i],
						"pc": pc,
						"pc_score": pc_scores[i, pc]
					})

			pc_scores_df = pd.DataFrame(pc_scores_list)

			#convert to matrix
			pc_scores_df = pc_scores_df.pivot(index="cell_id", columns="pc", values="pc_score").fillna(0)
			pc_scores_df.reset_index(inplace=True)

			
			#insert PC scores
			self.query_raw("DROP TABLE IF EXISTS PC_scores;")
			self.open_db()
			self.conn.register("pc_scores_df", pc_scores_df)
			self.conn.execute("CREATE TABLE PC_scores AS SELECT * FROM pc_scores_df;")
			self.close_db()

			
		else:

			#drop temp table
			self.query_raw("DROP TABLE IF EXISTS PC_scores_temp;")
		
			num_rows = self.query(f"SELECT {genes[0]} FROM X_standard_wide").shape[0]
			matrix = np.zeros((num_rows, n_pcs))

			#load the pcs into a dictionary so we can do O(1) lookups
			pc_loadings = {}
			for gene in genes:
				pc_loadings[gene] = pc_loadings_df[pc_loadings_df['gene'] == gene]["loading"].to_numpy()

			#iterate each gene and multiply the gene expression by the pc loading
			inc=1
			for gene in genes[:]:
				start_time = time.time()
				#gene_stand = self.query(f"SELECT {gene} FROM X_standard_wide").to_numpy()
				gene_stand = self.query(f"SELECT {gene} FROM X_standard_wide ORDER BY cell_id").to_numpy()

				gene_stand = gene_stand.reshape(-1) 
				matrix += gene_stand[:, None] * pc_loadings[gene] 
				if print_progress == True:
					print(f"PCs Processed: {inc} of {len(genes)} genes: {str(time.time() - start_time)} seconds")
				inc += 1

			#cell_ids = self.query(f"SELECT cell_id FROM obs")
			cell_ids = self.query(f"SELECT cell_id FROM X_standard_wide ORDER BY cell_id")["cell_id"].tolist()

			pc_loadings_temp_df = pd.DataFrame(matrix)
			pc_loadings_temp_df.insert(0, "cell_id", cell_ids)

			#save to PC_scores table
			self.query("DROP TABLE IF EXISTS PC_scores")
			self.open_db()
			self.conn.register("PC_scores_temp", pc_loadings_temp_df)
			self.conn.execute("CREATE TABLE PC_scores AS SELECT * FROM PC_scores_temp")
			self.close_db()


		print("\nPCA Calculation Complete\n")
			


	def save_highly_variable_genes(self, top_variable_genes=1000, gene_field="gene_names"):
		"""
		Save only the top highly variable genes from the 'var' table into table 'X'.
		
		Parameters:
			top_variable_genes (int): The number of top variable genes to save. Default is 1000.

		Notes:
			- The 'X' table is updated with only the top highly variable genes.
			- Consider running save_raw() before running this method.

		Example:
			asql = AnnSQL(db='db/pbmc.asql')			
			asql.save_highly_variable_genes(top_variable_genes=1000)
		"""

		
		genes = self.query(f"SELECT {gene_field} FROM var ORDER BY variance DESC LIMIT {top_variable_genes}")
		genes = genes[gene_field].tolist()

		query = f"""
		CREATE TABLE X_buffer AS
		SELECT cell_id, {', '.join(genes)}
		FROM X;
		"""
		self.query_raw(query)
		self.query_raw(f"DROP TABLE IF EXISTS X;")
		self.query_raw(f"ALTER TABLE X_buffer RENAME TO X")
		self.query_raw(f"DELETE FROM var WHERE {gene_field} NOT IN ({', '.join([f'\'{gene}\'' for gene in genes])});")
		print(f"X table updated with only HV genes.")


	def filter_by_gene_counts(self, min_gene_counts=None, max_gene_counts=None, gene_field="gene_names"):
		"""
		Filter the data by gene counts. This method removes cells with gene counts below the minimum threshold and above the maximum threshold.

		Parameters:
			min_gene_counts (int, optional): The minimum gene counts threshold. Genes with counts below this threshold will be removed. Defaults to None.
			max_gene_counts (int, optional): The maximum gene counts threshold. Genes with counts above this threshold will be removed. Defaults to None.
			gene_field (str, optional): The field name for gene names in the var database table. Defaults to "gene_names".

		Notes:
			- This method removes genes with counts below the minimum threshold and above the maximum threshold from the 'X' table.
			- The method also updates the 'var' table to reflect the changes in the 'X' table
			- Consider running save_raw() before running this method.

		Example:
			asql = AnnSQL(db='db/pbmc.asql')			
			asql.filter_by_gene_counts(min_gene_counts=200, max_gene_counts=45000)
		"""

		if min_gene_counts != None and max_gene_counts == None:
			genes = self.query(f"SELECT {gene_field} FROM var WHERE gene_counts > {min_gene_counts}")
			genes = genes[gene_field].tolist()
			
			query = f"""
			CREATE TABLE X_buffer AS
			SELECT cell_id, {', '.join(genes)}
			FROM X;
			"""
			self.query_raw(query)
			self.query_raw(f"DROP TABLE IF EXISTS X;")
			self.query_raw(f"ALTER TABLE X_buffer RENAME TO X")
			self.query_raw(f"DELETE FROM var WHERE gene_counts <= {min_gene_counts};")
			print(f"Removed genes with less than {min_gene_counts} from X table.")

		elif min_gene_counts != None and max_gene_counts != None:
			genes = self.query(f"SELECT {gene_field} FROM var WHERE gene_counts > {min_gene_counts} AND gene_counts < {max_gene_counts}")
			genes = genes[gene_field].tolist()

			#check X columns for total_counts. We need to keep this column for meow
			if 'total_counts' in self.query("SELECT * FROM X LIMIT 1").columns:
				genes.append('total_counts')

			query = f"""
			CREATE TABLE X_buffer AS
			SELECT cell_id, {', '.join(genes)}
			FROM X;
			"""
			self.query_raw(query)
			self.query_raw(f"DROP TABLE IF EXISTS X;")
			self.query_raw(f"ALTER TABLE X_buffer RENAME TO X")
			self.query_raw(f"DELETE FROM var WHERE gene_counts <= {min_gene_counts} OR gene_counts >= {max_gene_counts};")
			print(f"Removed genes with less than {min_gene_counts} and greater than {max_gene_counts} from X table.")

		
	def return_pca_scores_matrix(self, legacy=False):		
		"""
		Returns the PCA scores matrix in the form of a matrix with cell IDs as index, PCs as columns, and PC scores as values.

		Parameters:
			legacy (bool, optional): Whether to return the legacy PCA scores matrix. For version<1 compatibility. Defaults to False.

		Raises:
			ValueError: If the 'PC_scores' table is not found. Run calculate_pca() first.

		Returns:
			result (DataFrame): The PCA scores matrix with cell IDs as index, PCs as columns, and PC scores as values.
		"""

		if 'PC_scores' not in self.show_tables()['table_name'].tolist():
			raise ValueError('PC_scores table not found. Run calculate_pca() first')
		
		if legacy == True:
			return self.query("SELECT * FROM PC_scores").pivot(index="cell_id", columns="pc", values="pc_score").fillna(0)
		else:
			return self.query("SELECT * FROM PC_scores")

	def save_raw(self, table_name="X_raw"):
		"""
		Saves the raw data from the 'X' table into a new table named 'X_raw'. Helps in preserving the raw data before any modifications.

		Parameters:
			table_name (str): The name of the table to save the raw data into. Defaults to "X_raw".

		Example:
			asql = AnnSQL(db='db/pbmc.asql')			
			asql.save_raw()
		"""

		self.query_raw(f"DROP TABLE IF EXISTS X_raw;")
		self.query_raw(f"CREATE TABLE X_raw AS SELECT * FROM X;")
		print("X_raw table created from X.")

	def raw_to_X(self, table_name="X_raw"):
		"""
		Replaces the 'X' table with the contents of the specified raw data table and updates the 'var' table with gene names.

		Parameters:
			table_name (str): The name of the table containing the raw data to be renamed to 'X'. Defaults to "X_raw".

		Example:
			asql = AnnSQL(db='db/pbmc.asql')			
			asql.raw_to_X()
		"""

		self.query_raw(f"DROP TABLE IF EXISTS X;")
		self.query_raw(f"ALTER TABLE {table_name} RENAME TO X;")

		#empy the var table
		self.query_raw("DELETE FROM var;")

		#get all of the column names from the X table
		columns = self.query("DESCRIBE X")[1:]["column_name"].tolist()
		
		#insert the gene names into the var table
		values = [f"('{col}', '{col}')" for col in columns if col != "cell_id"]

		if values:
			query = f"INSERT INTO var (gene_names, gene_names_orig) VALUES {', '.join(values)};"
			self.query_raw(query)

		#we need to reset the obs table as well
		self.query_raw("DELETE FROM obs;")
		self.query_raw("INSERT INTO obs (cell_id) SELECT cell_id FROM X;")

		print("X table created from X_raw. Please note: X_raw table has been deleted.")

	def pca_variance_explained(self, show_plot=True, return_values=False):
		"""
		Calculate the variance explained by each principal component in a PCA analysis.

		Parameters:
			show_plot (bool, optional): Whether to display a bar plot showing the variance explained by each principal component. Default is True.
			return_values (bool, optional): Whether to return the variance explained values as a pandas Series. Default is False.

		Raises:
			ValueError: If PCA has not been calculated yet. Run calculate_pca() first.

		Notes:
			- The variance explained by each principal component is calculated as the ratio of the eigenvalue of the component to the total sum of eigenvalues.
			- Plotting the variance explained by each principal component can help in determining the number of components to retain in the analysis.

		Returns:
			pandas Series (optional): The variance explained by each principal component, if return_values is True.
		"""

		if 'PC_eigenvalues' not in self.show_tables()['table_name'].tolist():
			raise ValueError('PCA not found. Run calculate_pca() first')

		eigen_df = self.query("SELECT * FROM PC_eigenvalues ORDER BY pc;")
		total_variance = eigen_df['eigenvalue'].sum()
		eigen_df['variance_explained'] = eigen_df['eigenvalue'] / total_variance

		#make the pc 1-based
		eigen_df['pc'] = eigen_df['pc'] + 1

		if show_plot == True:
			plt.figure(figsize=(12, 6))
			g = sns.barplot(x='pc', y='variance_explained', data=eigen_df)
			g.set_title("PCA Variance Explained")
			g.set_xlabel("Principal Component")
			g.set_ylabel("Variance Explained")
			g.set_xticklabels(g.get_xticklabels(), rotation=45, ha="right", fontsize=10)
			plt.show()
		
		if return_values == True:
			return eigen_df['variance_explained']

	def plot_pca(self,PcX=1, PcY=2):
		"""
		Plots the PCA scores for two principal components.

		Parameters:
			PcX (int): The index of the first principal component (default: 1).
			PcY (int): The index of the second principal component (default: 2).

		Notes: 
			- Indexes are 1-based. (offset by +1 for convienence, or confusion. Whichever way you look at it.)

		Raises:
			ValueError: If 'PC_eigenvalues' table is not found. Run calculate_pca() first.

		Example:
			asql = AnnSQL(db='db/pbmc.asql')			
			asql.plot_pca(PcX=1, PcY=2)		
		"""

		if 'PC_eigenvalues' not in self.show_tables()['table_name'].tolist():
			raise ValueError('PCA not found. Run calculate_pca() first')

		#get the scores in a matrix
		pca_scores = self.return_pca_scores_matrix()

		#match the 1-based index
		pc1 = PcX-1
		pc2 = PcY-1

		g = sns.scatterplot(x=pca_scores[pc1], y=pca_scores[pc2], data=pca_scores)
		g.set_title(f"PC{PcY} vs PC{PcX}")
		g.set_xlabel(f"PC{PcX}")
		g.set_ylabel(f"PC{PcY}")


	def plot_highly_variable_genes(self, top_variable_genes=2000):
		"""
		Plots the highly variable genes based on their variance and mean expression.

		This method checks if the 'variance' column exists in the 'var' table. If not, it calculates the variable genes.
		It then ensures the 'hv' column exists in the 'var' table, creating or resetting it as necessary.
		Finally, it updates the 'hv' column to highlight the top variable genes and plots the variance vs mean expression.

		Parameters:
			top_variable_genes (int): The number of top variable genes to highlight. Default is 2000.

		Example:
			asql = AnnSQL(db='db/pbmc.asql')			
			asql.plot_highly_variable_genes(top_variable_genes=2000)
		"""

		#does the variance column exist in var?
		if 'variance' not in self.query("SELECT * FROM var LIMIT 1").columns:
			print("Variance not found. Running calculate_variable_genes...")
			self.calculate_variable_genes(chunk_size=200, print_progress=False, 
										save_var_names=False, save_top_variable_genes=top_variable_genes)

		#does the column hv exist in var?
		if 'hv' not in self.query("SELECT * FROM var LIMIT 1").columns:
			self.query_raw("ALTER TABLE var ADD COLUMN hv INT DEFAULT 0")
		else:
			self.query_raw("ALTER TABLE var DROP COLUMN hv")
			self.query_raw("ALTER TABLE var ADD COLUMN hv INT DEFAULT 0")

		self.query_raw(f"UPDATE var SET hv=1 WHERE gene_names IN (SELECT gene_names FROM var ORDER BY variance DESC LIMIT {top_variable_genes})")
		g=sns.scatterplot(data=self.query("SELECT gene_mean, log10(variance) as variance, hv FROM var"), x="gene_mean", y="variance", hue="hv", s=2)
		g.set_title(f"Showing all genes variance vs mean expression\nTop {top_variable_genes} genes highlighted")
		g.set_xlabel("Mean Expression")
		g.set_ylabel("Variance (log10)")

	def plot_total_counts(self):
		"""
		Plots the total UMI counts for each cell.

		Notes:
			- This method plots a violin plot showing the distribution of total UMI counts across cells.
			- The 'total_counts' column is used from the 'obs

		Example:
			asql = AnnSQL(db='db/pbmc.asql')			
			asql.plot_total_counts()
		"""

		if 'total_counts' not in self.query("SELECT * FROM obs LIMIT 1").columns:
			print("Total counts not found. Run method calculate_total_counts() first.")

		g=sns.violinplot(x="total_counts", data=self.query("SELECT total_counts FROM obs"))
		g.set_title("Total UMI Counts")
		g.set_xlabel("Total UMI Counts")
		
	def plot_gene_counts(self):
		"""
		Plots the gene counts for each cell.

		Notes:
			- This method plots a violin plot showing the distribution of gene counts across cells.
			- The 'gene_counts' column is used from the 'var' table.

		Example:
			asql = AnnSQL(db='db/pbmc.asql')			
			asql.plot_gene_counts()
		"""

		if 'gene_counts' not in self.query("SELECT * FROM var LIMIT 1").columns:
			print("Gene counts not found. Run method calculate_gene_counts() first.")

		g=sns.violinplot(x="gene_counts", data=self.query("SELECT gene_counts FROM var"))
		g.set_title("Gene Counts")
		g.set_xlabel("Gene Counts")
	

	def calculate_umap(self, n_neighbors=15, min_dist=0.5, n_components=2, metric='euclidean'):
		"""
		Calculates the Uniform Manifold Approximation and Projection (UMAP) for the data stored in the 'PC_scores' table.

		Parameters:
			n_neighbors (int, default=15): The size of local neighborhood (in terms of number of neighboring sample points) used for manifold approximation.
			min_dist (float, default=0.5): The effective minimum distance between embedded points.
			n_components (int, default=2): The number of dimensions of the UMAP embedding.
			metric (str, default='euclidean'): The metric to use for distance computation.
		
		Notes:
			- The 'PC_scores' table is used to calculate the UMAP embedding.
			- The UMAP embedding is saved in a new table named 'umap_embeddings'.
			- The 'umap_embeddings' table contains the UMAP1 and UMAP2 coordinates for each cell.
		
		Example:
			asql = AnnSQL(db='db/pbmc.asql')			
			asql.calculate_umap(n_neighbors=15, min_dist=0.5, n_components=2, metric='euclidean')
		"""

		pca_scores = self.return_pca_scores_matrix().drop(columns=['cell_id']).values

		#drop column cell_id
		#pca_scores = pca_scores.drop(columns=['cell_id'])

		reducer = umap.UMAP(n_neighbors=n_neighbors,
							min_dist=min_dist,
							n_components=n_components,
							metric=metric)
		self.umap_embedding = pd.DataFrame(reducer.fit_transform(pca_scores))
		self.umap_embedding.rename(columns={0: "UMAP1", 1: "UMAP2"}, inplace=True)		
		
		self.query_raw("DROP TABLE IF EXISTS umap_embeddings;")
		self.open_db()
		self.conn.register("umap_embeddings_df", self.umap_embedding)
		self.conn.execute("CREATE TABLE umap_embeddings AS SELECT * FROM umap_embeddings_df;")
		self.close_db()
		print("UMAP embedding calculated.")


	def plot_umap(self, color_by=None, palette='viridis', title=None, legend_location=None, annotate=False, counts_table="X"):
		"""
		Plots the UMAP projection of the data from the 'umap_embeddings' table.

		Parameters:
			color_by (str, optional): The column name in the 'obs' or 'var' table to use for coloring the cells. Defaults to None.
			palette (str, optional): The color palette to use for coloring the cells. Defaults to 'viridis'.
			title (str, optional): The title of the plot. Defaults to None.
			legend_location (str, optional): The location of the legend in the plot. Defaults to None.
			annotate (bool, optional): Whether to annotate the plot with cluster centers. Defaults to False.
		
		Raises:
			ValueError: If 'umap_embeddings' table is not found.

		Notes:
			- The 'umap_embeddings' table is used to plot the UMAP projection.
			- The 'color_by' parameter can be a column name in the 'obs' or 'var' table to color the cells based on that column.
			- The 'palette' parameter can be any valid seaborn color palette.

		Example:
			asql = AnnSQL(db='db/pbmc.asql')			
			asql.plot_umap(color_by='leiden_clusters' title='UMAP Projection', annotate=True)
		"""


		if 'umap_embeddings' not in self.show_tables()['table_name'].tolist():
			raise ValueError('UMAP embedding not found. Run calculate_umap() first')
		
		umap_embedding = self.query("SELECT * FROM umap_embeddings")
		
		obs_values = None
		if counts_table == "X":
			if color_by in self.query("SELECT * FROM obs LIMIT 1").columns:
				obs_values = self.query(f"SELECT {color_by} FROM obs ORDER BY cell_id")
			elif color_by in self.query("SELECT gene_names FROM var")["gene_names"].values:
				obs_values = self.query(f"SELECT {color_by} FROM X ORDER BY cell_id")
		else:
			obs_values = self.query(f"SELECT {color_by} FROM {counts_table} WHERE cell_id IN (SELECT cell_id FROM obs) ORDER BY cell_id")
		
		df = umap_embedding.copy()
		if obs_values is not None:
			df[color_by] = obs_values.iloc[:, 0].values
		
		plt.figure(figsize=(8, 6))
		
		g = sns.scatterplot(data=df, x="UMAP1", y="UMAP2", hue=color_by, palette=palette)
		
		if legend_location is not None:
			plt.legend(loc=legend_location)
		else:
			plt.legend([],[], frameon=False)

		#continuous values, add a colorbar .
		if obs_values is not None and pd.api.types.is_numeric_dtype(df[color_by]):
			sc = g.collections[0]
			plt.colorbar(sc)
		
		if annotate:
			import matplotlib.patheffects as path_effects
			cluster_centers = df.groupby(color_by)[['UMAP1', 'UMAP2']].mean()
			for label, row in cluster_centers.iterrows():
				text = plt.text(row['UMAP1'], row['UMAP2'], str(label), fontsize=12, ha='center', va='center', color='white')
				text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])

		plt.title(title if title is not None else (f"UMAP Projection | {color_by}" if color_by else "UMAP Projection"))
		plt.xlabel("UMAP1")
		plt.ylabel("UMAP2")
		plt.show()


	def calculate_leiden_clusters(self, resolution=1.0, n_neighbors=30):
		"""
		Performs Leiden clustering on the data using the UMAP embeddings.

		Parameters:
			resolution (float, default=1.0): The resolution parameter for the Leiden algorithm.
			n_neighbors (int, default=30): The size of the neighborhood to consider for the Leiden algorithm.
		
		Raises:
			ValueError: If 'umap_embeddings' table is not found. Run calculate_umap() first.

		Notes:
			- The Leiden clustering results are saved in the 'obs' table under the column 'leiden_clusters'.
			- The 'leiden_clusters' column is cast as a string for plotting purposes.
		
		Example:
			asql = AnnSQL(db='db/pbmc.asql')			
			asql.calculate_leiden_clusters(resolution=1.0, n_neighbors=30)
		"""


		if 'PC_scores' not in self.show_tables()['table_name'].tolist():
			raise ValueError('PCA scores table not found. Run calculate_pca() first')
		
		from sklearn.neighbors import kneighbors_graph
		from sknetwork.clustering import Leiden
		
		pca_scores = self.return_pca_scores_matrix().drop(columns=['cell_id']).values
		cell_ids = self.return_pca_scores_matrix()["cell_id"].tolist()

		knn_graph = kneighbors_graph(pca_scores, n_neighbors=30, mode="connectivity", include_self=False)
		leiden = Leiden(resolution=resolution, random_state=42)
		labels = leiden.fit_predict(knn_graph)

		leiden_clusters = pd.DataFrame({
			"cell_id": cell_ids,
			"leiden_clusters": labels
		})

		self.query_raw("DROP TABLE IF EXISTS leiden_clusters;")
		self.open_db()
		self.conn.register("leiden_clusters_df", leiden_clusters)
		self.conn.execute("CREATE TABLE leiden_clusters AS SELECT * FROM leiden_clusters_df;")
		self.conn.execute("""
			ALTER TABLE obs ADD COLUMN IF NOT EXISTS leiden_clusters INT;
			UPDATE obs 
			SET leiden_clusters = (
				SELECT leiden_clusters 
				FROM leiden_clusters 
				WHERE obs.cell_id = leiden_clusters.cell_id
			);
		""")
		self.conn.execute("DROP TABLE IF EXISTS leiden_clusters;")
		#cast as string for plotting :shrug:
		self.query_raw("ALTER TABLE obs ALTER COLUMN leiden_clusters TYPE TEXT;")
		self.close_db()

		print("Leiden clustering complete. Clusters saved in 'obs' as 'leiden_clusters'.")


	def add_observations(self, obs_key, obs_values={}, match_on="ledien_clusters"):
		"""
		Adds observation values to the 'obs' table based on a key-value pair.

		Parameters:
			obs_key (str): The key to add to the 'obs' table.
			obs_values (dict): A dictionary of cell IDs and corresponding values for the observation key.
			match_on (str, optional): The column in the 'obs' table to match the cell IDs on. Defaults to 'leiden_clusters'.
		
		Notes:
			- The 'obs' table is updated with the new observation key and values.
			- The 'obs_values' dictionary should have cell IDs as keys and the corresponding values for the observation key.
			- The 'match_on' parameter specifies the column in the 'obs' table to match the cell IDs on.
		
		Example:
			asql = AnnSQL(db='db/pbmc.asql')	
			cell_annotations = {'0': 'T cell', '1': 'B cell'}		
			asql.add_observations(obs_key='cell_type', match_on="leiden_clusters" obs_values=cell_annotations)
		"""

		#if the obs_key doesn't exist, add it.
		if obs_key not in self.query("SELECT * FROM obs LIMIT 1").columns:
			self.query_raw(f"ALTER TABLE obs ADD COLUMN {obs_key} TEXT;")
		else:
			#drop the column if it exists
			self.query_raw(f"ALTER TABLE obs DROP COLUMN {obs_key};")
			self.query_raw(f"ALTER TABLE obs ADD COLUMN {obs_key} TEXT;")


		#iterate over the obs_values in key, value pairs
		for key, value in obs_values.items():
			self.query_raw(f"UPDATE obs SET {obs_key} = '{value}' WHERE {match_on} = '{key}'")
		
		print(f"{obs_key} added to obs table!\nobs_values keys matched on {match_on}.")


	@staticmethod
	def t_cdf(t_val: float, df: float) -> float:
		return tdist.cdf(t_val, df)

	def adjusted_p_value(self):
		"""
		Adjusts the p-values for multiple testing using the Benjamini-Hochberg procedure. 

		Notes:
			- The adjusted p-values are saved in the 'diff_expression' table under the column 'adj_pval'.
			- The 'diff_expression' table must have the columns 'gene' and 'pval' for this method to work.
		
		Example:
			asql = AnnSQL(db='db/pbmc.asql')			
			asql.adjusted_p_value()
		"""


		results = self.query("SELECT * FROM diff_expression")
		adj_p_val_results = stats.false_discovery_control(results["pval"].to_numpy())
		results
		adj_p_val_results
		df_adj = pd.DataFrame({
			"gene": results["gene"],
			"adj_pval": adj_p_val_results
		})

		self.open_db()
		self.conn.register("temp_adj", df_adj)
		self.conn.execute("""
		UPDATE diff_expression
		SET adj_pval = temp_adj.adj_pval
		FROM temp_adj
		WHERE diff_expression.gene = temp_adj.gene;
		""")
		self.conn.unregister("temp_adj")
		

	def calculate_differential_expression(self, obs_key=None, group1_value=None, group2_value=None, name="None", drop_table=False, marker_genes=False, gene_field="gene_names"):
		"""
		Calculates differential expression between two groups of cells based on the 'X' table and the 'obs' table. This is performed using a t-test.

		Parameters:
			obs_key (str): The observation key to use for grouping the cells.
			group1_value (str): The value of the observation key for group 1.
			group2_value (str): The value of the observation key for group 2.
			name (str, optional): The name of the differential expression analysis. Defaults to "None".
			drop_table (bool, optional): Whether to drop the 'diff_expression' table if it already exists. Defaults to False.
			marker_genes (bool, optional): Whether to calculate marker genes for the groups. This is used for the calculate_marker_genes method. Defaults to False.
		
		Raises:
			ValueError: If 'obs_key', 'group1_value', and 'group2_value' are not provided.
			ValueError: If 'obs_key' is not provided for marker_genes=True.
		
		Notes:
			- The differential expression results are saved in the 'diff_expression' table.
			- The 'diff_expression' table contains the columns 'name', 'group1', 'group2', 'gene', 'tstat', 'logfc', 'df', 'pval', and 'adj_pval'.
			- The 'name' parameter is used to identify the differential expression analysis.
			- The 'marker_genes' parameter can be set to True to calculate marker genes for the groups.
		
		Example:
			asql = AnnSQL(db='db/pbmc.asql')			
			asql.calculate_differential_expression(obs_key='leiden_clusters', group1_value='0', group2_value='1')
		"""


		if marker_genes == False and (group1_value is None or group2_value is None or obs_key is None):
			raise ValueError("obs_key, group1 and group2 must be provided.")
		elif marker_genes == True and (obs_key is None or group1_value is None):
			raise ValueError("obs_key must be provided for marker_genes=True.")

		genes = self.query(f"SELECT {gene_field} FROM var;")[gene_field].values.tolist()

		if drop_table == True:
			self.query_raw("DROP TABLE IF EXISTS diff_expression;")
		
		self.query_raw("CREATE TABLE IF NOT EXISTS diff_expression (name TEXT, group1 TEXT, group2 TEXT, gene TEXT, tstat DOUBLE, logfc DOUBLE, df DOUBLE, pval DOUBLE, adj_pval DOUBLE);")
		self.open_db()
		self.conn.create_function("t_cdf", AnnSQL.t_cdf, [float, float], float)
		
		if (marker_genes == True):
			groups = self.conn.execute(f"SELECT DISTINCT {obs_key} FROM obs").df()[obs_key].values.tolist()
			groups = [str(group) for group in groups]
			groups_str = "', '".join(groups)
			in_condition = f"IN ('{groups_str}')"
			group_by_condition = f"GROUP BY CASE WHEN obs.{obs_key} = '{group1_value}' THEN '{group1_value}' ELSE 'ALL' END"
			select_condition = f"CASE WHEN obs.{obs_key} = '{group1_value}' THEN '{group1_value}' ELSE 'ALL' END AS group_label,"
			group2_value = "ALL"
		else:
			in_condition = f"IN ('{group1_value}', '{group2_value}')"
			group_by_condition = f"GROUP BY obs.{obs_key}"
			select_condition = f"obs.{obs_key} AS group_label,"
		i=1
		for gene in genes:
			print(f"Complete gene {i} of {len(genes)}")
			i=i+1
			query = f"""
			WITH stats AS (
				SELECT 
					{select_condition}
					AVG(X.{gene})  AS mean_1,
					COUNT(X.{gene}) AS count_1,
					var_samp(X.{gene}) AS variance_1
				FROM X
				INNER JOIN obs ON X.cell_id = obs.cell_id 
				WHERE {obs_key} {in_condition}
				{group_by_condition}
			),
			calc AS (
				SELECT
					'{name}' AS name,
					'{group1_value}' AS group1,
					'{group2_value}' AS group2,
					'{gene}' AS gene,
					(s1.mean_1 - s2.mean_1) / SQRT(s1.variance_1/s1.count_1 + s2.variance_1/s2.count_1) AS tstat,
					LOG((s1.mean_1 + 1e-10) / (s2.mean_1+ + 1e-10)) / LOG(2) AS logfc,
					POWER(s1.variance_1/s1.count_1 + s2.variance_1/s2.count_1, 2) /
					(
						POWER(s1.variance_1/s1.count_1, 2)/(s1.count_1 - 1) +
						POWER(s2.variance_1/s2.count_1, 2)/(s2.count_1 - 1)
					) AS df,
					2.0 * (1.0 - t_cdf(ABS(tstat), df)) AS pval,
				FROM stats s1
				CROSS JOIN stats s2
				WHERE s1.group_label = '{group1_value}'
				AND s2.group_label = '{group2_value}'
			)
			INSERT INTO diff_expression SELECT *, 1 FROM calc;
			"""
			self.conn.execute(query)
		self.close_db()
		self.adjusted_p_value()
		print(f"DE Calculation Complete.")

	def plot_differential_expression(self, pvalue_threshold=0.05, logfc_threshold=1.0, group1=None, group2=None, title=None, filter_name=None):
		"""
		Plots the differential expression as a volcano plot from the 'diff_expression' table.

		Parameters:
			pvalue_threshold (float, optional): The p-value threshold for significance. Default is 0.05.
			logfc_threshold (float, optional): The log fold change threshold for significance. Default is 1.0.
			group1 (str, optional): The value of group 1 for differential expression. Defaults to None.
			group2 (str, optional): The value of group 2 for differential expression. Defaults to None.
			title (str, optional): The title of the plot. Defaults to None.
			filter_name (str, optional): The name of the differential expression analysis to filter on. Defaults to None.
		
		Notes:
			- The 'diff_expression' table must contain the columns 'logfc', 'adj_pval', 'group1', and 'group2'.
			- The volcano plot shows the log fold change on the x-axis and the negative log10 of the adjusted p-value on the y-axis.
			- Significant genes are highlighted in red based on the p-value and log fold change thresholds.
		
		Example:
			asql = AnnSQL(db='db/pbmc.asql')			
			asql.plot_differential_expression(pvalue_threshold=0.05, logfc_threshold=1.0, group1='0', group2='1', title='DE Analysis')
		"""


		if group1 is None or group2 is None:
			df = self.query("SELECT * FROM diff_expression")
		else:
			df = self.query(f"SELECT * FROM diff_expression WHERE group1 = '{group1}' AND group2 = '{group2}'")

		df['neg_log10_adj_pval'] = -np.log10(df['adj_pval'])
		significant = (df['adj_pval'] < pvalue_threshold) & (np.abs(df['logfc']) >= logfc_threshold)
		plt.figure(figsize=(8, 6))
		plt.scatter(df['logfc'], df['neg_log10_adj_pval'], color='grey', alpha=0.7, label='Not significant')
		plt.scatter(df.loc[significant, 'logfc'], df.loc[significant, 'neg_log10_adj_pval'], 
					color='red', alpha=0.7, label='Significant')
		plt.axhline(-np.log10(pvalue_threshold), color='blue', linestyle='--', label=f'p={pvalue_threshold}')
		plt.axvline(logfc_threshold, color='blue', linestyle='--', label=f'logFC={logfc_threshold}')
		plt.axvline(-(logfc_threshold), color='blue', linestyle='--', label=f'logFC=-{logfc_threshold}')
		plt.xlabel('Log2 Fold Change')
		plt.ylabel('-Log10 Adjusted p-value')
		plt.title(f'Differential Expression | {group1} and {group2}')
		if title is not None:
			plt.title(title)
		plt.legend()
		plt.show()
		print(f"Query the results with:\n\"SELECT * FROM diff_expression WHERE group1='{group1}' and group2='{group2}'\".")


	def calculate_marker_genes(self, obs_key="leiden_clusters", table_name="X"):
		"""
		Calculates marker genes for each group based on differential expression analysis.

		Parameters:
			obs_key (str): The observation key to use for grouping the cells.
			table_name (str, optional): The name of the table to use for the analysis. Defaults to "X".
		
		Notes:
			- This method calculates marker genes for each group based on the differential expression analysis  using a t-test.
			- The marker genes are saved in the 'diff_expression' table under the name 'Markers'.
		
		Example:
			asql = AnnSQL(db='db/pbmc.asql')			
			asql.calculate_marker_genes(obs_key='leiden_clusters', table_name='X')
		"""

		self.open_db()
		groups = self.conn.execute(f"SELECT DISTINCT {obs_key} FROM obs").df()[obs_key].values.tolist()
		for group in groups:
			print(f"Calculating marker genes for {obs_key}: {group}")
			self.calculate_differential_expression(obs_key=obs_key, group1_value=group, group2_value="ALL", name="Markers", drop_table=False, marker_genes=True)
		print(f"Marker genes calculation complete.")
		print(f"Query the results with:\n\"SELECT * FROM diff_expression WHERE name='Markers'\".")


	def plot_marker_genes(self, obs_key="leiden_clusters", columns=2):
		"""
		Plots the top marker genes for each cluster based on the differential expression analysis.

		Parameters:
			obs_key (str): The observation key to use for grouping the cells.
			columns (int, optional): The number of columns in the plot. Defaults to 2.

		Raises:
			ValueError: If 'diff_expression' table is not found. Run calculate_marker_genes() first.

		Notes:
			- This method plots the top marker genes for each cluster based on the differential expression analysis.
			- The 'diff_expression' table must contain the columns 'gene', 'tstat', and 'adj_pval'.
		
		Example:
			asql = AnnSQL(db='db/pbmc.asql')			
			asql.plot_marker_genes(obs_key='leiden_clusters', columns=2)
		"""

		if 'diff_expression' not in self.show_tables()['table_name'].tolist():
			raise ValueError('Differential expression not found. Run calculate_differential_expression() first')
		
		groups = self.query(f"SELECT DISTINCT({obs_key}) FROM obs ORDER BY {obs_key};")[obs_key].tolist()

		num_rows = (len(groups) + 2) // columns
		fig, axes = plt.subplots(num_rows, columns, figsize=(15, 5 * num_rows))
		axes = axes.flatten()

		for idx, group in enumerate(groups):
			results = self.query(f"""
			SELECT row_number() OVER (ORDER BY tstat DESC) as id, * FROM diff_expression 
				WHERE name='Markers' AND group1 = '{group}' AND group2 = 'ALL' 
				AND adj_pval < 0.05 
				AND (logfc > 0.5 OR logfc < 0.5) 
				ORDER BY tstat DESC LIMIT 20;
			""")
			
			ax = axes[idx]
			sns.scatterplot(x=results["id"], y=results["tstat"], ax=ax)
			for i, txt in enumerate(results["gene"]):
				ax.annotate(txt, (results["id"][i], results["tstat"][i]), xytext=(-5, 12), textcoords='offset points', rotation=90)
			max_tstat = results["tstat"].max()
			ax.set_ylim(0, max_tstat + 20)
			ax.set_ylabel("t-statistic")
			ax.set_xlabel("Rank")
			ax.set_title(f"Top marker genes for cluster {group} vs all")

		#remove unused subplots
		for j in range(idx + 1, len(axes)):
			fig.delaxes(axes[j])

		plt.tight_layout()
		plt.show()

	def get_marker_genes(self, obs_key="leiden_clusters", group=None):
		"""
		Returns the top marker genes for a specific group based on the differential expression analysis. 

		Parameters:
			obs_key (str): The observation key to use for grouping the cells.
			group (str): The value of the observation key for the group.
		
		Raises:
			ValueError: If 'diff_expression' table is not found. Run calculate_marker_genes() first.
			ValueError: If 'group' is not provided.

		Returns:
			pandas DataFrame: The top marker genes for the specified group.

		Example:
			asql = AnnSQL(db='db/pbmc.asql')			
			asql.get_marker_genes(obs_key='leiden_clusters', group='0')
		"""

		if group is None:
			raise ValueError("Group must be provided.")
		if 'diff_expression' not in self.show_tables()['table_name'].tolist():
			raise ValueError('Differential expression not found. Run calculate_marker_genes() first')

		return self.query(f"""
		SELECT * FROM diff_expression 
		WHERE name='Markers' AND group1 = '{group}' AND group2 = 'ALL' 
		AND adj_pval < 0.05 
		AND (logfc > 0.5 OR logfc < 0.5) 
		ORDER BY tstat DESC;
		""")


	def write_adata(self, filename="export.h5ad"):
		"""
		Converts a on-disk AnnSQL object to an AnnData object and writes it to a file.

		Parameters:
			filename (str): The name of the file to write the AnnData object to. Defaults to "export.h5ad".
		
		Notes:
			- This method writes the AnnData object to a file in the HDF5 format.
			- The AnnData object contains the 'X', 'obs', 'var', 'PC_scores', 'diff_expression', and 'umap_embeddings' tables.
		
		Example:
			asql = AnnSQL(db='db/pbmc.asql')			
			asql.write_adata(filename='pbmc.h5ad')
		"""

		import anndata as ad
		adata = ad.AnnData(X=self.query("SELECT * EXCLUDE(cell_id) FROM X"), 
							obs=self.query("SELECT * FROM obs"))
		adata.var = self.query("SELECT * FROM var")
		if 'PC_scores' in self.show_tables()['table_name'].tolist():
			adata.obsm["X_pca"] = self.return_pca_scores_matrix().values
		if 'diff_expression' in self.show_tables()['table_name'].tolist():
			adata.uns["diff_expression"] = self.query("SELECT * FROM diff_expression")
		if 'umap_embeddings' in self.show_tables()['table_name'].tolist():
			adata.obsm["X_umap"] = self.query("SELECT * FROM umap_embeddings").values
		adata.write(filename)
		print(f"AnnData object written to {filename}.")