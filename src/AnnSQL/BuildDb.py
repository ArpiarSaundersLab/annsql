import scanpy as sc
import pandas as pd
import polars as pl
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import duckdb
import os 
import json
import time
import gc
import psutil
import warnings
from memory_profiler import profile
from scipy.sparse import issparse
warnings.filterwarnings('ignore')
import sys

class BuildDb:
	

	sql_reserved_keywords = [
		'add', 'all', 'alter', 'and', 'any', 'as', 'asc', 'between', 'by', 'case', 'cast', 'check', 
		'column', 'create', 'cross', 'current_date', 'current_time', 'default', 'delete', 'desc', 
		'distinct', 'drop', 'else', 'exists', 'false', 'for', 'foreign', 'from', 'full', 'group', 
		'having', 'in', 'inner', 'insert', 'interval', 'into', 'is', 'join', 'left', 'like', 'limit', 
		'not', 'null', 'on', 'or', 'order', 'outer', 'primary', 'references', 'right', 'select', 
		'set', 'table', 'then', 'to', 'true', 'union', 'unique', 'update', 'values', 'when', 'where'
	]

	def __init__(self, conn=None, db_path=None, db_name=None, adata=None, create_all_indexes=False, create_basic_indexes=False, convenience_view=True, chunk_size=5000,	make_buffer_file=False, print_output=True, layers=["X", "obs", "var", "var_names", "obsm", "varm", "obsp", "uns"], db_config={}):
		"""
		Initializes the BuildDb object. This object is used to create a database from an AnnData object byway of the MakeDb object.

		Attributes:
			sql_reserved_keywords (list): List of SQL reserved keywords.

		Parameters:
			conn (optional): Connection object to the database.
			db_path (optional): Path to the database file.
			db_name (optional): Name of the database.
			adata (optional): AnnData object.
			create_all_indexes (optional): Flag indicating whether to create all indexes.
			create_basic_indexes (optional): Flag indicating whether to create basic indexes.
			convenience_view (optional): Flag indicating whether to create a convenience view.
			chunk_size (optional): Size of the chunks for processing the data.
			make_buffer_file (optional): Flag indicating whether to create a buffer file.
			print_output (optional): Flag indicating whether to print output.
			layers (optional): List of layers to include in the database.
		Returns:
			None
		"""
		self.adata = adata
		self.conn = conn
		self.db_path = db_path
		self.db_name = db_name
		self.create_all_indexes = create_all_indexes
		self.create_basic_indexes = create_basic_indexes
		self.convenience_view = convenience_view
		self.layers = layers
		self.chunk_size = chunk_size
		self.make_buffer_file = make_buffer_file
		self.print_output = print_output
		self.db_config = db_config
		self.build()
		if "uns" in self.layers: #not recommended for large datasets
			self.build_uns_layer()

	def build(self):
		"""
		Build the database tables for the AnnSQL object.
		This method creates and inserts data into the following tables:
		Tables:
		
			X: Contains cell_id as VARCHAR and var_names_df columns as FLOAT.
			obs: Contains the observation data.
			var_names: Contains the gene names.
			var: Contains the variable data.
			obsm: Contains the observation matrix data.
			varm: Contains the variable matrix data.
			obsp: Contains the observation sparse matrix data.
		
		The method also creates indexes on the tables based on the specified layers.

		Parameters:
			None
		Returns:
			None
		"""
		obs_df = self.adata.obs.reset_index()
		var_names = self.adata.var_names
		var = self.adata.var
		var_names_df = pd.DataFrame(var_names)
		var_names_df.columns = ['gene']
		obs_df.columns = ['cell_id'] + list(obs_df.columns[1:])
		
		#The var_names_make_unique appears doesn't handle case sensitively. SQL requires true unique column names with case sensitivity
		#AnnData would treat something like Gad1 and gad1 as different genes
		var_names_upper = pd.DataFrame(var_names).apply(lambda x: x.str.upper())
		var_names = list(var_names)
		start_time = time.time()
		unique_counter = {}
		for i in range(len(var_names_upper)):
			if var_names_upper.duplicated()[i] == True:
				if var_names_upper.iloc[i][0] in unique_counter:
					unique_counter[var_names_upper.iloc[i][0] ]+=1
					var_names[i] = var_names[i] + f"_{unique_counter[var_names_upper.iloc[i][0] ]}"
				else:
					unique_counter[var_names_upper.iloc[i][0] ] = 1
					var_names[i] = var_names[i] + f"_{unique_counter[var_names_upper.iloc[i][0] ]}"
		end_time = time.time()

		if self.print_output == True:
			print("Time to make var_names unique: ", end_time-start_time)

		#clean the column names for SQL
		var_names_clean = [self.replace_special_chars(col) for col in var_names]

		#Create X with cell_id as varchar and var_names_df columns as float
		#Note: casting as float expecting floating point calculations in future (e.g. normalization)
		#consider making the OG duckdb cast a parameter for users who want to store as int
		start_time = time.time()
		self.conn.execute("CREATE TABLE X (cell_id VARCHAR,	{} )".format(', '.join([f"{self.replace_special_chars(col)} FLOAT" for col in var_names])))
		end_time = time.time()
		if self.print_output == True:
			print("Time to create X table structure: ", end_time-start_time)

		# self.conn.close()
		# self.conn = None

		#handles backed mode or if chunk size <= number of rows
		if self.adata.isbacked or self.chunk_size <= self.adata.shape[0]:
			if "X" in self.layers:
				chunk_size = self.chunk_size 
				if os.path.exists(f"{self.db_path}{self.db_name}_X.parquet"):
					os.remove(f"{self.db_path}{self.db_name}_X.parquet")
				print(f"Starting chunked mode X table data insert. Total rows: {self.adata.shape[0]}")
				writer = None

				for start in range(0, self.adata.shape[0], chunk_size):
					start_time = time.time()
					end = min(start + chunk_size, self.adata.shape[0])

					#reconnect to the database :/
					#self.conn = duckdb.connect(f"{self.db_path}{self.db_name}.asql", config=self.db_config)

					if issparse(self.adata.X) == True:
						X_chunk_df = np.array(self.adata[start:end].X.todense())
					else:
						X_chunk_df = self.adata[start:end].X
					
					X_chunk_df = pl.DataFrame({"cell_id": self.adata.obs.index[start:end],**{name: X_chunk_df[:, idx] for idx, name in enumerate(var_names_clean)}})
					self.conn.register(f"X_chunk_df", X_chunk_df)
					
					if self.make_buffer_file == False:
						self.conn.execute("BEGIN TRANSACTION;")
						self.conn.execute("SET preserve_insertion_order = false;")
						self.conn.execute("INSERT INTO X SELECT * FROM X_chunk_df;")
						self.conn.execute("COMMIT;")
					else:
						table = X_chunk_df.to_arrow()
						if writer is None:
							writer = pq.ParquetWriter(f"{self.db_path}{self.db_name}_X.parquet", table.schema)
						writer.write_table(table)
					
					self.conn.unregister(f"X_chunk_df")

					# self.conn.close()
					# self.conn = None

					del X_chunk_df
					X_chunk_df = None
					gc.collect()
					print(f"Processed chunk {start}-{end-1} in {time.time()-start_time} seconds")

				if writer is not None:
					writer.close()
					
				if self.make_buffer_file == True:
					start_time = time.time()
					print("\nToo close for missiles, switching to guns\nCreating X table from buffer file.\nThis may take a while...")
					self.conn.execute(f"INSERT INTO X SELECT * FROM read_parquet('{self.db_path}{self.db_name}_X.parquet')")
					print(f"Time to create X table from buffer: {time.time()-start_time}")
					os.remove(f"{self.db_path}{self.db_name}_X.parquet")
				if self.print_output == True:
					print(f"Finished inserting X data.")

			else:
				print("Skipping X layer")

		else:
			if "X" in self.layers:
				start_time = time.time()
				#duckdb gives an error when registering sparse polars dataframe, so we need to convert to dense (for now)
				if issparse(self.adata.X):
					if self.print_output == True:
						print("Converting sparse to dense")
					self.adata.X = self.adata.X.todense()

				# #is this an in-memory database?
				if self.db_path == None:
					self.conn.register("X", pl.DataFrame({"cell_id": self.adata.obs.index,**{name: self.adata.X[:, idx] for idx, name in enumerate(var_names_clean)}})) 
				else:
					X_df = self.conn.register("X_df", pl.DataFrame({"cell_id": self.adata.obs.index,**{name: self.adata.X[:, idx] for idx, name in enumerate(var_names_clean)}})) 
					self.conn.execute("BEGIN TRANSACTION;")
					self.conn.execute("SET preserve_insertion_order = false;")
					self.conn.execute("INSERT INTO X SELECT * FROM X_df")
					self.conn.execute("COMMIT;")

				end_time = time.time()
				gc.collect()
				if self.print_output == True:
					print("Time to insert X data: ", end_time-start_time )
			else:
				print("Skipping X layer")



		#these tables usually are not as large as X and can be inserted in one go
		if "obs" in self.layers:
			self.conn.register('obs_df', obs_df)
			self.conn.execute("CREATE OR REPLACE TABLE obs AS SELECT * FROM obs_df")
			self.conn.unregister('obs_df')
			if self.print_output == True:
				print("Finished inserting obs data")
		else:
			print("Skipping obs layer")

		if "var_names" in self.layers:
			self.conn.register('var_names_df', var_names_df)
			self.conn.execute("CREATE OR REPLACE TABLE var_names AS SELECT * FROM var_names_df")
			self.conn.unregister('var_names_df')
			if self.print_output == True:
				print("Finished inserting var_names data")
		else:
			print("Skipping var_names layer")

		if "var" in self.layers:
			var["gene_names_orig"] = var.index
			if 'gene_name' in var.columns:
				var["gene_names"] =  [col for col in var.gene_name]
			else:
				var["gene_names"] = [self.replace_special_chars(col) for col in var_names]

			var = var.reset_index(drop=True)
			self.conn.register('var_df', var)
			self.conn.execute("CREATE OR REPLACE TABLE var AS SELECT * FROM var_df")
			self.conn.unregister('var_df')
			if self.print_output == True:
				print("Finished inserting var data")
		else:
			print("Skipping var layer")

		if "obsm" in self.layers:
			for key in self.adata.obsm.keys():
				obsm_df = pd.DataFrame(self.adata.obsm[key])
				self.conn.register(f'obsm_{key}_df', obsm_df)
				self.conn.execute(f"CREATE OR REPLACE TABLE obsm_{key} AS SELECT * FROM obsm_{key}_df")
				self.conn.unregister(f'obsm_{key}_df')
			if self.print_output == True:
				print("Finished inserting obsm data")
		else:
			print("Skipping obsm layer")


		if "varm" in self.layers:
			for key in self.adata.varm.keys():
				varm_df = pd.DataFrame(self.adata.varm[key])
				self.conn.register(f'varm_{key}_df', varm_df)
				self.conn.execute(f"CREATE OR REPLACE TABLE varm_{key} AS SELECT * FROM varm_{key}_df")
				self.conn.unregister(f'varm_{key}_df')
			if self.print_output == True:
				print("Finished inserting varm data")
		else:
			print("Skipping varm layer")

		if "obsp" in self.layers:
			for key in self.adata.obsp.keys():
				obsp_df = pd.DataFrame(self.adata.obsp[key].toarray())
				self.conn.register(f'obsp_{key}_df', obsp_df)
				self.conn.execute(f"CREATE OR REPLACE TABLE obsp_{key} AS SELECT * FROM obsp_{key}_df")
				self.conn.unregister(f'obsp_{key}_df')
			if self.print_output == True:
				print("Finished inserting obsp data")
		else:
			print("Skipping obsp layer")

		#indexes (Warning: resource intensive. only recommended for small datasets)
		if self.create_all_indexes == True:
			if "X" in self.layers:
				for column in X_df.columns:
					try:
						self.conn.execute(f'CREATE INDEX idx_{column.replace("-", "_").replace(".", "_")}_X ON X ("{column}")')
					except:
						print(f'Could not create index on {column} for X')

			if "obs" in self.layers:
				for column in obs_df.columns:
					try:
						self.conn.execute(f'CREATE INDEX idx_{column.replace("-", "_").replace(".", "_")}_obs ON obs ("{column}")')
					except:
						print(f'Could not create index on {column} for obs')

		#basic indexes
		if self.create_basic_indexes == True:
			if "obs" in self.layers:
				self.conn.execute("CREATE INDEX idx_obs_cell_id ON obs (cell_id)")
			if "X" in self.layers:
				self.conn.execute("CREATE INDEX idx_X_cell_id ON X (cell_id)")

		#view for convenience (not recommended for large datasets)
		if self.convenience_view == True and "X" in self.layers and "obs" in self.layers:
			self.conn.execute("CREATE VIEW adata AS SELECT * FROM obs JOIN X ON obs.cell_id = X.cell_id")

	def make_json_serializable(self,value):
		"""
		Converts a given value into a JSON serializable format.

		Parameters:
			value (any): The value to be converted.
		Returns:
			JSON: The converted value in a JSON serializable format.
		"""
		
		if isinstance(value, np.ndarray):
			return value.tolist()
		elif isinstance(value, (np.int64, np.int32)):
			return int(value)
		elif isinstance(value, (np.float64, np.float32)):
			return float(value)
		elif isinstance(value, dict):
			return {k: self.make_json_serializable(v) for k, v in value.items()}
		elif isinstance(value, list):
			return [self.make_json_serializable(v) for v in value]  
		else:
			return value  

	def build_uns_layer(self):
		"""
		Builds the uns_raw table in the database and inserts the uns data.
		This method creates the uns_raw table in the database if it doesn't exist and inserts the uns data into the table.
		The uns data is a dictionary containing key-value pairs. The keys are stored in the 'key' column, the serialized
		values are stored in the 'value' column, and the data type of each value is stored in the 'data_type' column.
		If an error occurs during table creation or data insertion, the error message is printed to the console.

		Parameters:
			None
		Returns:
			None
		"""

		try:
			self.conn.execute("CREATE TABLE uns_raw (key TEXT, value TEXT, data_type TEXT)")
		except Exception as e:
			print(f"Error creating uns_raw table: {e}")
		for key, value in self.adata.uns.items():
			try:
				serialized_value = self.make_json_serializable(value)
			except TypeError as e:
				print(f"Error serializing key {key}: {e}")
				continue
			if isinstance(value, dict):
				data_type = 'dict'
			elif isinstance(value, list):
				data_type = 'list'
			elif isinstance(value, (int, float, str)):
				data_type = 'scalar'
			elif isinstance(value, np.ndarray):
				data_type = 'array'
				value = value.tolist()
			else:
				data_type = 'unknown'
			try:
				self.conn.execute("INSERT INTO uns_raw VALUES (?, ?, ?)", (key, serialized_value, data_type))
			except Exception as e:
				print(f"Error inserting key {key}: {e}")

	def replace_special_chars(self, string):
		"""
		Replaces special characters in a string with underscores.

		Parameters:
			string (str): The input string to be processed.
		Returns:
			str: The processed string with special characters replaced by underscores.
		"""

		if string.lower() in self.sql_reserved_keywords:
			string = "r_"+string #prefix reserved keywords with r_. The OG can be found in gene_names_orig column
		if string[0].isdigit():
			return 'n'+string.replace("-", "_").replace(".", "_")
		else:
			return string.replace("-", "_").replace(".", "_").replace("(", "_").replace(")", "_").replace(",", "_").replace(" ", "_").replace("/", "_")

	def determine_buffer_status(self):
		"""
		Determines the buffer status based on the available memory.
		If the total available memory is less than or equal to 20GB, a buffer is used
		and building the database will take longer. To disable the buffer, set the
		`make_buffer_file` parameter to `False`.
		
		Parameters:
			None
		Returns:
			None
		"""
		mem = psutil.virtual_memory()
		if mem.total <= 20 * 1024 ** 3:
			print("=========================================================================")
			print("Low memory system detected. (<=20GB)")
			print("Using a buffer. Building the Db will take longer...")
			print("To disable this, explicitly set the parameter to: make_buffer_file=False")
			print("=========================================================================")
			self.make_buffer_file = True
		else:
			self.make_buffer_file = False