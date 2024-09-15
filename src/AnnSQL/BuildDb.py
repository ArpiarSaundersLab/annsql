import scanpy as sc
import pandas as pd
import numpy as np
import duckdb
import os 
import json

class BuildDb:
	def __init__(self, conn=None, adata=None, create_all_indexes=False, add_uns=True, convenience_view=True):
		self.adata = adata
		self.conn = conn
		self.create_all_indexes = create_all_indexes
		self.convenience_view = convenience_view
		self.build()
		if add_uns == True: #not recommended for large datasets
			self.build_uns_layer()

	def build(self):
		obs_df = self.adata.obs.reset_index()
		var_names = self.adata.var_names
		var = self.adata.var
		var_names_df = pd.DataFrame(var_names)
		var_names_df.columns = ['gene']
		obs_df.columns = ['cell_id'] + list(obs_df.columns[1:])
		
		#handle backed mode
		if self.adata.isbacked:
			first_chunk = self.adata.X[:1].toarray() if hasattr(self.adata.X[:1], 'toarray') else self.adata.X[:1]
			X_df = pd.DataFrame(first_chunk, columns=var_names)
			cell_id_df = pd.DataFrame(obs_df['cell_id'][:1]).reset_index(drop=True)
			X_df = pd.concat([cell_id_df, X_df], axis=1)
			X_df.columns = ['cell_id'] + list(X_df.columns[1:])

			self.conn.register('X_df', X_df)
			self.conn.execute("CREATE TABLE X AS SELECT * FROM X_df")
			self.conn.unregister('X_df')

			chunk_size = 2000 
			print(f"Starting backed mode db creation. Total rows: {self.adata.shape[0]}")
			for start in range(1, self.adata.shape[0], chunk_size):
				end = min(start + chunk_size, self.adata.shape[0])
				X_chunk = self.adata.X[start:end].toarray() if hasattr(self.adata.X[start:end], 'toarray') else self.adata.X[start:end]
				X_chunk_df = pd.DataFrame(X_chunk, columns=var_names)
				cell_id_chunk_df = pd.DataFrame(obs_df['cell_id'][start:end]).reset_index(drop=True)
				X_chunk_df = pd.concat([cell_id_chunk_df, X_chunk_df], axis=1)
				self.conn.register('X_chunk_df', X_chunk_df)
				self.conn.execute("INSERT INTO X SELECT * FROM X_chunk_df")
				self.conn.unregister('X_chunk_df')
				print(f"Inserted chunk {start}-{end}")

			#insert obs, var, var_names, obsm, varm, obsp

			print(f"Finished inserting data in chunks.")

		else:
			
			#not backed mode
			X_df = pd.DataFrame(self.adata.X.toarray() if hasattr(self.adata.X, 'toarray') else self.adata.X,
								columns=var_names)
			cell_id_df = pd.DataFrame(obs_df['cell_id']).reset_index(drop=True)
			X_df = pd.concat([cell_id_df, X_df], axis=1)
			X_df.columns = ['cell_id'] + list(X_df.columns[1:])
			self.conn.register('X_df', X_df)
			self.conn.execute("CREATE TABLE X AS SELECT * FROM X_df")
			self.conn.unregister('X_df')

		#these tables are usually not as large and don't require chunking
		self.conn.register('obs_df', obs_df)
		self.conn.execute("CREATE TABLE obs AS SELECT * FROM obs_df")
		self.conn.unregister('obs_df')
		self.conn.register('var_names_df', var_names_df)
		self.conn.execute("CREATE TABLE var_names AS SELECT * FROM var_names_df")
		self.conn.unregister('var_names_df')
		self.conn.register('var_df', var)
		self.conn.execute("CREATE TABLE var AS SELECT * FROM var_df")
		self.conn.unregister('var_df')

		for key in self.adata.obsm.keys():
			obsm_df = pd.DataFrame(self.adata.obsm[key])
			self.conn.register(f'obsm_{key}_df', obsm_df)
			self.conn.execute(f"CREATE TABLE obsm_{key} AS SELECT * FROM obsm_{key}_df")
			self.conn.unregister(f'obsm_{key}_df')

		for key in self.adata.varm.keys():
			varm_df = pd.DataFrame(self.adata.varm[key])
			self.conn.register(f'varm_{key}_df', varm_df)
			self.conn.execute(f"CREATE TABLE varm_{key} AS SELECT * FROM varm_{key}_df")
			self.conn.unregister(f'varm_{key}_df')

		for key in self.adata.obsp.keys():
			obsp_df = pd.DataFrame(self.adata.obsp[key].toarray())
			self.conn.register(f'obsp_{key}_df', obsp_df)
			self.conn.execute(f"CREATE TABLE obsp_{key} AS SELECT * FROM obsp_{key}_df")
			self.conn.unregister(f'obsp_{key}_df')

		#indexes (resource intensive)
		if self.create_all_indexes == True:
			for column in X_df.columns:
				try:
					self.conn.execute(f'CREATE INDEX idx_{column.replace("-", "_").replace(".", "_")}_X ON X ("{column}")')
				except:
					print(f'Could not create index on {column} for X')

			for column in obs_df.columns:
				try:
					self.conn.execute(f'CREATE INDEX idx_{column.replace("-", "_").replace(".", "_")}_obs ON obs ("{column}")')
				except:
					print(f'Could not create index on {column} for obs')

		#basic indexes
		self.conn.execute("CREATE INDEX idx_obs_cell_id ON obs (cell_id)")
		self.conn.execute("CREATE INDEX idx_X_cell_id ON X (cell_id)")

		#view for convenience (not recommended for large datasets)
		if self.convenience_view == True:
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
			