from AnnSQL.BuildDb import BuildDb
import scanpy as sc
import pandas as pd
import duckdb
import os 

class MakeDb:
	def __init__(self, adata=None, 	db_name=None, db_path="db/", create_all_indexes=False, create_basic_indexes=False, convenience_view=True, chunk_size=5000,make_buffer_file=False, layers=["X", "obs", "var", "var_names", "obsm", "varm", "obsp", "uns"]):
		"""
		Initializes the MakeDb object. This object is used to create a database from an AnnData object.

		Parameters:
			- adata: AnnData object, (optional): The AnnData object to be used for creating the database.
			- db_name: str, (optional): The name of the database.
			- db_path: str, (optional): The path where the database will be created. Must have a trailing slash.
			- create_all_indexes: bool, (optional): Whether to create indexes for all layers in the database.
			- create_basic_indexes: bool, (optional): Whether to create indexes for basic layers in the database.
			- convenience_view: bool, (optional): Whether to create a convenience view for the database.
			- chunk_size: int, (optional):	The number of cells to be processed in each chunk.
			- make_buffer_file: bool, (optional): Whether to create a buffer file for storing intermediate data. Necessary for low memory systems (<=12Gb).
			- layers: list of str, (optional): The layers to be included in the database.
		Returns:
			None
		"""
		self.adata = adata
		self.db_name = db_name
		self.db_path = db_path
		if not self.db_path.endswith('/'): #add trailing slash
			self.db_path += '/'
		self.layers = layers
		self.create_all_indexes = create_all_indexes
		self.create_basic_indexes = create_basic_indexes
		self.convenience_view = convenience_view
		self.chunk_size = chunk_size
		self.make_buffer_file = make_buffer_file
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
		BuildDb(adata=self.adata, conn=self.conn, create_all_indexes=self.create_all_indexes, create_basic_indexes=self.create_basic_indexes, convenience_view=self.convenience_view, layers=self.layers, chunk_size=self.chunk_size, db_path=self.db_path, db_name=self.db_name, make_buffer_file=self.make_buffer_file)
		self.conn.close()
