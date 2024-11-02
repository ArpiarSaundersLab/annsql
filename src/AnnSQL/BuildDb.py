import scanpy as sc
import pandas as pd
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
warnings.filterwarnings('ignore')

class BuildDb:

	sql_reserved_keywords = [
		'add', 'all', 'alter', 'and', 'any', 'as', 'asc', 'between', 'by', 'case', 'cast', 'check', 
		'column', 'create', 'cross', 'current_date', 'current_time', 'default', 'delete', 'desc', 
		'distinct', 'drop', 'else', 'exists', 'false', 'for', 'foreign', 'from', 'full', 'group', 
		'having', 'in', 'inner', 'insert', 'interval', 'into', 'is', 'join', 'left', 'like', 'limit', 
		'not', 'null', 'on', 'or', 'order', 'outer', 'primary', 'references', 'right', 'select', 
		'set', 'table', 'then', 'to', 'true', 'union', 'unique', 'update', 'values', 'when', 'where'
	]

	def __init__(self, conn=None, db_path=None, db_name=None, adata=None, create_all_indexes=False, create_basic_indexes=False, convenience_view=True, chunk_size=5000,	make_buffer_file=False,	layers=["X", "obs", "var", "var_names", "obsm", "varm", "obsp", "uns"]):
		"""
		Initializes the BuildDb object. This object is used to create a database from an AnnData object byway of the MakeDb object.

		Parameters:
			- conn (optional): Connection object to the database.
			- db_path (optional): Path to the database file.
			- db_name (optional): Name of the database.
			- adata (optional): AnnData object.
			- create_all_indexes (optional): Flag indicating whether to create all indexes.
			- create_basic_indexes (optional): Flag indicating whether to create basic indexes.
			- convenience_view (optional): Flag indicating whether to create a convenience view.
			- chunk_size (optional): Size of the chunks for processing the data.
			- make_buffer_file (optional): Flag indicating whether to create a buffer file.
			- layers (optional): List of layers to include in the database.
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
		self.build()
		if "uns" in self.layers: #not recommended for large datasets
			self.build_uns_layer()
	
	def determine_buffer_status(self):
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


	def print_memory_usage(self,chunk_num, id=1):
		process = psutil.Process()
		mem_info = process.memory_info()
		print(f"Chunk {chunk_num} - id {id} Memory usage: {mem_info.rss / (1024 ** 2):.2f} MB")

	def build(self):
		obs_df = self.adata.obs.reset_index()
		var_names = self.adata.var_names
		var = self.adata.var
		var_names_df = pd.DataFrame(var_names)
		var_names_df.columns = ['gene']
		obs_df.columns = ['cell_id'] + list(obs_df.columns[1:])
		
		#The var_names_make_unique appears doesn't handle case sensitively. SQL requires true unique column names
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
		print("Time to make var_names unique: ", end_time-start_time)

		#clean the column names for SQL
		var_names_clean = [self.replace_special_chars(col) for col in var_names]

		#Create X with cell_id as varchar and var_names_df columns as float
		#Note: casting as float expecting floating point calculations in future (e.g. normalization)
		#consider making the OG duckdb cast a parameter for users who want to store as int
		start_time = time.time()
		self.conn.execute("CREATE TABLE X (cell_id VARCHAR,	{} )".format(', '.join([f"{self.replace_special_chars(col)} FLOAT" for col in var_names])))
		end_time = time.time()
		print("Time to create X table structure: ", end_time-start_time)
		
		#handles backed mode
		if self.adata.isbacked:
			if "X" in self.layers:
				chunk_size = self.chunk_size 
				if os.path.exists(f"{self.db_path}{self.db_name}_X.parquet"):
					os.remove(f"{self.db_path}{self.db_name}_X.parquet")
				print(f"Starting backed mode X table data insert. Total rows: {self.adata.shape[0]}")
				writer = None
				for start in range(0, self.adata.shape[0], chunk_size):
					start_time = time.time()
					end = min(start + chunk_size, self.adata.shape[0])

					X_chunk_df = self.adata[start:end].to_df()
					X_chunk_df = X_chunk_df.reset_index()
					X_chunk_df.columns = ['cell_id'] + list(var_names_clean)

					if self.make_buffer_file == False:
						self.conn.execute("BEGIN TRANSACTION;")
						self.conn.execute("SET preserve_insertion_order = false;")
						self.conn.execute("INSERT INTO X SELECT * FROM X_chunk_df;")
						self.conn.execute("COMMIT;")
					else:
						table = pa.Table.from_pandas(X_chunk_df)
						if writer is None:
							writer = pq.ParquetWriter(f"{self.db_path}{self.db_name}_X.parquet", table.schema)
						writer.write_table(table)

					del X_chunk_df
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

				print(f"Finished inserting X data.")

			else:
				print("Skipping X layer")

		else:
			if "X" in self.layers:
				start_time = time.time()
				X_df = self.adata.to_df()
				X_df = X_df.reset_index()
				X_df.columns = ['cell_id'] + list(var_names_clean)
				self.conn.execute("BEGIN TRANSACTION;")
				self.conn.execute("SET preserve_insertion_order = false;")
				self.conn.execute("INSERT INTO X SELECT * FROM X_df;")
				self.conn.execute("COMMIT;")
				end_time = time.time()
				print("Time to insert X data: ", end_time-start_time )
			else:
				print("Skipping X layer")


		#these tables usually are not as large as X and can be inserted in one go
		if "obs" in self.layers:
			self.conn.register('obs_df', obs_df)
			self.conn.execute("CREATE TABLE obs AS SELECT * FROM obs_df")
			self.conn.unregister('obs_df')
			print("Finished inserting obs data")
		else:
			print("Skipping obs layer")

		if "var_names" in self.layers:
			self.conn.register('var_names_df', var_names_df)
			self.conn.execute("CREATE TABLE var_names AS SELECT * FROM var_names_df")
			self.conn.unregister('var_names_df')
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
			self.conn.execute("CREATE TABLE var AS SELECT * FROM var_df")
			self.conn.unregister('var_df')
			print("Finished inserting var data")
		else:
			print("Skipping var layer")

		if "obsm" in self.layers:
			for key in self.adata.obsm.keys():
				obsm_df = pd.DataFrame(self.adata.obsm[key])
				self.conn.register(f'obsm_{key}_df', obsm_df)
				self.conn.execute(f"CREATE TABLE obsm_{key} AS SELECT * FROM obsm_{key}_df")
				self.conn.unregister(f'obsm_{key}_df')
			print("Finished inserting obsm data")
		else:
			print("Skipping obsm layer")


		if "varm" in self.layers:
			for key in self.adata.varm.keys():
				varm_df = pd.DataFrame(self.adata.varm[key])
				self.conn.register(f'varm_{key}_df', varm_df)
				self.conn.execute(f"CREATE TABLE varm_{key} AS SELECT * FROM varm_{key}_df")
				self.conn.unregister(f'varm_{key}_df')
			print("Finished inserting varm data")
		else:
			print("Skipping varm layer")

		if "obsp" in self.layers:
			for key in self.adata.obsp.keys():
				obsp_df = pd.DataFrame(self.adata.obsp[key].toarray())
				self.conn.register(f'obsp_{key}_df', obsp_df)
				self.conn.execute(f"CREATE TABLE obsp_{key} AS SELECT * FROM obsp_{key}_df")
				self.conn.unregister(f'obsp_{key}_df')
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
		if string.lower() in self.sql_reserved_keywords:
			string = "r_"+string #prefix reserved keywords with r_. The OG can be found in gene_names_orig column
		if string[0].isdigit():
			return 'n'+string.replace("-", "_").replace(".", "_")
		else:
			return string.replace("-", "_").replace(".", "_").replace("(", "_").replace(")", "_").replace(",", "_").replace(" ", "_")