from BuildDb import BuildDb
import scanpy as sc
import pandas as pd
import duckdb
import warnings
import logging
import os 

class AnnSQL:
	def __init__(self, adata=None, db=None, create_all_indexes=False, layers=["X", "obs", "var", "var_names", "obsm", "varm", "obsp", "uns"]):
		self.adata = adata
		self.db = db
		self.create_all_indexes = create_all_indexes
		self.layers = layers
		self.validate_params()
		if self.db is None:
			self.build_db()
		else:
			self.open_db()

	def validate_params(self):
		if self.adata is None and self.db is None:
			raise ValueError('Both adata and db parameters not defined. Select an option')
		if self.adata is not None:
			if not isinstance(self.adata, sc.AnnData):
				raise ValueError('adata must be a scanpy AnnData object')
		if self.db is not None:
			if not os.path.exists(self.db):
				raise ValueError('The db provided doesn\'t exist. Please check the path')
		if self.db is not None and self.adata is None:
			warnings.warn('Warning: No adata object provided. return_type="adata" is disabled.')

	def open_db(self):
		if self.db is not None:
			self.conn = duckdb.connect(self.db)

	def close_db(self):
		if self.db is not None:
			self.conn.close()

	def build_db(self):
		self.conn = duckdb.connect(':memory:')
		db = BuildDb(adata=self.adata, conn=self.conn, create_all_indexes=self.create_all_indexes, layers=self.layers)
		self.conn = db.conn

	def query(self, query, return_type='pandas'):
		if return_type not in ['pandas', 'adata', 'parquet']:
			raise ValueError('return_type must be either pandas, parquet or adata')
		if 'UPDATE' in query.upper() or 'DELETE' in query.upper() or 'INSERT' in query.upper():
			raise ValueError('UPDATE, DELETE, and INSERT detected. Please use update_query() instead')

		self.open_db()
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
			return self.adata[result_df["cell_id"]]

	def query_raw(self, query):
		self.open_db()
		result = self.conn.execute(query)
		self.close_db()
		return result

	def update_query(self, query):
		if 'SELECT' in query.upper():
			raise ValueError('SELECT detected. Please use query() instead')
		try:
			self.open_db()
			self.conn.execute(query)
			self.close_db()
			print("Query Successful")
		except Exception as e:
			print("Update Query Error:", e)

	def show_tables(self):
		self.open_db()
		result = self.conn.execute("SELECT table_name FROM information_schema.tables  WHERE table_schema='main'").df()
		self.close_db()
		return result

	def show_settings(self):
		self.open_db()
		result = self.conn.execute("SELECT * FROM duckdb_settings()").df()
		self.close_db()
		return result

	def export_parquet(self):
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