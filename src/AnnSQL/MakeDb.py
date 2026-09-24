from AnnSQL.BuildDb import BuildDb
import scanpy as sc
import pandas as pd
import duckdb
import os 

class MakeDb:
	def __init__(self, adata=None, 	db_name=None, db_path="db/", create_all_indexes=False, create_basic_indexes=False, convenience_view=True, chunk_size=10000,make_buffer_file=False, layers=["X", "obs", "var", "var_names", "obsm", "varm", "obsp", "uns"], print_output=True, db_config={}, delete_existing_db=False, row_group_size=None, block_size=None):
		"""
		Initializes the MakeDb object. This object is used to create a database from an AnnData object by using the BuildDb method.

		Args:
			adata (AnnData, optional): The AnnData object to be used for creating the database.
			db_name (str, optional): The name of the database.
			db_path (str, optional): The path where the database will be created. Must have a trailing slash.
			create_all_indexes (bool, optional): Whether to create indexes for all layers in the database.
			create_basic_indexes (bool, optional): Whether to create indexes for basic layers in the database.
			convenience_view (bool, optional): Whether to create a convenience view for the database.
			chunk_size (int, optional): The number of cells to be processed in each chunk.
			make_buffer_file (bool, optional): Whether to create a buffer file for storing intermediate data. Necessary for low memory systems (<=12Gb).
			layers (list of str, optional): The layers to be included in the database.
			print_output (bool, optional): Whether to print output messages.
			row_group_size (int, optional): DuckDB row group size for the database file. Smaller values (e.g. 8192) bound how many
				rows of X are buffered in memory before being written to disk. Defaults to None (DuckDB default of 122880).
			block_size (int, optional): DuckDB storage block size in bytes (power of two, e.g. 16384). DuckDB reserves about one block
				per column while inserting, so with ~30k genes the default 262144 costs ~8-15GB. Defaults to None (DuckDB default).
				For wide datasets use together with row_group_size and a memory limit, e.g. db_config={'memory_limit': '4GB'}.

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
		self.print_output = print_output
		self.db_config = db_config
		self.delete_existing_db = delete_existing_db
		self.row_group_size = row_group_size
		self.block_size = block_size
		self.validate_params()
		self.build_db()

	def validate_params(self):
		"""
		Validates the parameters required for creating a database.

		Raises:
			ValueError: If `db_name` is not provided or is not a string.
			ValueError: If `db_path` is not provided or is not a valid system path.
			ValueError: If `adata` is provided but is not an instance of `scanpy.AnnData`.
		"""

		if self.db_name is None:
			raise ValueError('db_name is required and must be a string')
		if self.db_path is None:
			raise ValueError('db_path is required and must be a valid system path')
		if self.adata is not None:
			if not isinstance(self.adata, sc.AnnData):
				raise ValueError('adata must be a scanpy AnnData object')
		if self.row_group_size is not None:
			#duckdb requires the row group size to be a positive multiple of its vector size (2048)
			if not isinstance(self.row_group_size, int) or isinstance(self.row_group_size, bool) or self.row_group_size <= 0:
				raise ValueError('row_group_size must be a positive integer')
			if self.row_group_size % 2048 != 0:
				raise ValueError(f'row_group_size must be a multiple of 2048 (duckdb vector size), got {self.row_group_size}')
		if self.block_size is not None:
			#duckdb accepts a power of two between its minimum and maximum block size
			if not isinstance(self.block_size, int) or isinstance(self.block_size, bool) or self.block_size <= 0:
				raise ValueError('block_size must be a positive integer')
			if self.block_size & (self.block_size - 1) != 0:
				raise ValueError(f'block_size must be a power of two, got {self.block_size}')
			if not 16384 <= self.block_size <= 262144:
				raise ValueError(f'block_size must be between 16384 and 262144, got {self.block_size}')

	def create_db(self):
		"""
		Creates a new database if it does not already exist.

		Raises:
			ValueError: If the database already exists.
		Returns:
			None
		"""
		if self.delete_existing_db == True:
			if os.path.exists(self.db_path+self.db_name+'.asql'):
				os.remove(self.db_path+self.db_name+'.asql')
				if self.print_output:
					print('Deleted existing database: '+self.db_path+self.db_name+'.asql')

		if os.path.exists(self.db_path+self.db_name+'.asql'):
			raise ValueError('The database'+ self.db_path+self.db_name+'  exists already.')
		else:
			if not os.path.exists(self.db_path):
				os.makedirs(self.db_path)
			attach_options = []
			if self.block_size is not None:
				attach_options.append(f"BLOCK_SIZE {self.block_size}")
			if self.row_group_size is not None:
				attach_options.append(f"ROW_GROUP_SIZE {self.row_group_size}")

			if len(attach_options) == 0:
				self.conn = duckdb.connect(self.db_path+self.db_name+'.asql', config=self.db_config)
			else:
				#block and row group sizes can only be set when attaching a database file
				db_file = (self.db_path+self.db_name+'.asql').replace("'", "''")
				self.conn = duckdb.connect(config=self.db_config)
				self.conn.execute(f"ATTACH '{db_file}' AS asql_db ({', '.join(attach_options)})")
				self.conn.execute("USE asql_db")

	def build_db(self):
		"""
		Builds the database by creating it, executing the necessary SQL statements, and closing the connection.

		Parameters:
			None
		Returns:
			None
		"""
		self.create_db()
		BuildDb(adata=self.adata, conn=self.conn, create_all_indexes=self.create_all_indexes, create_basic_indexes=self.create_basic_indexes, convenience_view=self.convenience_view, layers=self.layers, chunk_size=self.chunk_size, db_path=self.db_path, db_name=self.db_name, make_buffer_file=self.make_buffer_file, print_output=self.print_output, db_config=self.db_config)
		self.conn.close()
