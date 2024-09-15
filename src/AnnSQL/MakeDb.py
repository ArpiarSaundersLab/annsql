from BuildDb import BuildDb
import scanpy as sc
import pandas as pd
import duckdb
import os 

class MakeDb:
	def __init__(self, adata=None, db_name=None, db_path="db/", create_all_indexes=False, add_uns=True, convenience_view=True):
		self.adata = adata
		self.db_name = db_name
		self.db_path = db_path
		self.create_all_indexes = create_all_indexes
		self.add_uns = add_uns
		self.convenience_view = convenience_view
		self.validate_params()
		self.build_db()

	def validate_params(self):
		if self.db_name is None:
			raise ValueError('db_name is required and must be a string')
		if self.db_path is None:
			raise ValueError('db_path is required and must be a valid system path')
		if self.adata is not None:
			if not isinstance(self.adata, sc.AnnData):
				raise ValueError('adata must be a scanpy AnnData object')

	def create_db(self):
		if os.path.exists(self.db_path+self.db_name+'.asql'):
			raise ValueError('The database'+ self.db_path+self.db_name+'  exists already.')
		else:
			if not os.path.exists(self.db_path):
				os.makedirs(self.db_path)
			self.conn = duckdb.connect(self.db_path+self.db_name+'.asql')

	def build_db(self):
		self.create_db()
		BuildDb(adata=self.adata, conn=self.conn, create_all_indexes=self.create_all_indexes, add_uns=self.add_uns, convenience_view=self.convenience_view)
		self.conn.close()
