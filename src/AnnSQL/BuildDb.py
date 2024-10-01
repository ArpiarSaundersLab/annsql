import scanpy as sc
import pandas as pd
import numpy as np
import duckdb
import os 
import json
import time
import warnings
warnings.filterwarnings('ignore')

class BuildDb:
	def __init__(self, 
				conn=None, 
				adata=None, 
				create_all_indexes=False, 
				create_basic_indexes=False,
				convenience_view=True,
				layers=["X", "obs", "var", "var_names", "obsm", "varm", "obsp", "uns"]):
		self.adata = adata
		self.conn = conn
		self.create_all_indexes = create_all_indexes
		self.create_basic_indexes = create_basic_indexes
		self.convenience_view = convenience_view
		self.layers = layers
		self.build()
		if "uns" in self.layers: #not recommended for large datasets
			self.build_uns_layer()

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

		#Create X with cell_id as varchar and var_names_df columns as float
		#import to do floating point calculations in future (e.g. normalization)
		start_time = time.time()
		self.conn.execute("""
			CREATE TABLE X (
				cell_id VARCHAR,
				{}
			)
		""".format(', '.join([f"{self.replace_special_chars(col)} FLOAT" for col in var_names])))
		end_time = time.time()
		print("Time to create X table schema: ", end_time-start_time)

		#handle backed mode
		if self.adata.isbacked:
			if "X" in self.layers:
				first_chunk = self.adata.X[:1].toarray() if hasattr(self.adata.X[:1], 'toarray') else self.adata.X[:1]
				X_df = pd.DataFrame(first_chunk, columns=var_names)
				cell_id_df = pd.DataFrame(obs_df['cell_id'][:1]).reset_index(drop=True)
				X_df = pd.concat([cell_id_df, X_df], axis=1)
				X_df.columns = ['cell_id'] + list(X_df.columns[1:])
				chunk_size = 5000 
				print(f"Starting backed mode data insert. Total rows: {self.adata.shape[0]}")
				for start in range(1, self.adata.shape[0], chunk_size):
					start_time = time.time()
					end = min(start + chunk_size, self.adata.shape[0])
					X_chunk = self.adata.X[start:end].toarray() if hasattr(self.adata.X[start:end], 'toarray') else self.adata.X[start:end]
					X_chunk_df = pd.DataFrame(X_chunk, columns=var_names)
					cell_id_chunk_df = pd.DataFrame(obs_df['cell_id'][start:end]).reset_index(drop=True)
					X_chunk_df = pd.concat([cell_id_chunk_df, X_chunk_df], axis=1)
					self.conn.register('X_chunk_df', X_chunk_df)
					self.conn.execute("INSERT INTO X SELECT * FROM X_chunk_df")
					self.conn.unregister('X_chunk_df')
					print(f"Inserted chunk {start}-{end} in {time.time()-start_time} seconds")
				print(f"Finished inserting data in chunks.")
			else:
				print("Skipping X layer")

		else:
			#not backed mode
			if "X" in self.layers:
				start_time = time.time()
				X_df = pd.DataFrame(self.adata.X.toarray() if hasattr(self.adata.X, 'toarray') else self.adata.X,
									columns=var_names)
				cell_id_df = pd.DataFrame(obs_df['cell_id']).reset_index(drop=True)
				X_df = pd.concat([cell_id_df, X_df], axis=1)
				X_df.columns = ['cell_id'] + list(X_df.columns[1:])
				self.conn.register('X_df', X_df)
				self.conn.execute("INSERT INTO X SELECT * FROM X_df")
				self.conn.unregister('X_df')
				end_time = time.time()
				print("Time to insert X data: ", end_time-start_time )
			else:
				print("Skipping X layer")


		#these tables are usually not as large and shouldn't require chunking
		if "obs" in self.layers:
			self.conn.register('obs_df', obs_df)
			self.conn.execute("CREATE TABLE obs AS SELECT * FROM obs_df")
			self.conn.unregister('obs_df')
		else:
			print("Skipping obs layer")

		if "var_names" in self.layers:
			self.conn.register('var_names_df', var_names_df)
			self.conn.execute("CREATE TABLE var_names AS SELECT * FROM var_names_df")
			self.conn.unregister('var_names_df')
		else:
			print("Skipping var_names layer")

		if "var" in self.layers:
			var["gene_names_orig"] = var.index
			var["gene_names"] = [self.replace_special_chars(col) for col in var_names]
			var = var.reset_index(drop=True)
			self.conn.register('var_df', var)
			self.conn.execute("CREATE TABLE var AS SELECT * FROM var_df")
			self.conn.unregister('var_df')
		else:
			print("Skipping var layer")

		if "obsm" in self.layers:
			for key in self.adata.obsm.keys():
				obsm_df = pd.DataFrame(self.adata.obsm[key])
				self.conn.register(f'obsm_{key}_df', obsm_df)
				self.conn.execute(f"CREATE TABLE obsm_{key} AS SELECT * FROM obsm_{key}_df")
				self.conn.unregister(f'obsm_{key}_df')
		else:
			print("Skipping obsm layer")


		if "varm" in self.layers:
			for key in self.adata.varm.keys():
				varm_df = pd.DataFrame(self.adata.varm[key])
				self.conn.register(f'varm_{key}_df', varm_df)
				self.conn.execute(f"CREATE TABLE varm_{key} AS SELECT * FROM varm_{key}_df")
				self.conn.unregister(f'varm_{key}_df')
		else:
			print("Skipping varm layer")

		if "obsp" in self.layers:
			for key in self.adata.obsp.keys():
				obsp_df = pd.DataFrame(self.adata.obsp[key].toarray())
				self.conn.register(f'obsp_{key}_df', obsp_df)
				self.conn.execute(f"CREATE TABLE obsp_{key} AS SELECT * FROM obsp_{key}_df")
				self.conn.unregister(f'obsp_{key}_df')
		else:
			print("Skipping obsp layer")

		#indexes (resource intensive. only recommended for small datasets)
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
		reserved_words = ['cell_id','CAST','SELECT','UPDATE']
		if string in reserved_words:
			string = "r_"+string
		if string[0].isdigit():
			return 'n'+string.replace("-", "_").replace(".", "_")
		else:
			return string.replace("-", "_").replace(".", "_")



