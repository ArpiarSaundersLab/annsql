import unittest
import AnnSQL as AnnSQL
from AnnSQL.MakeDb import MakeDb
import scanpy as sc
import os
import time
import warnings
warnings.filterwarnings('ignore')


class TestDatabase(unittest.TestCase):
	def setUp(self):
		self.adata = sc.datasets.pbmc68k_reduced()
		self.db_path = "tests/db/"
		self.db_name = "pbmc68k_reduced"
		self.db_file = os.path.join(self.db_path, f"{self.db_name}.asql")

	def make_wide_adata(self, genes, cells=100):
		#pbmc68k_reduced has only 765 genes, so the wide-data path needs a synthetic matrix
		import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
		X = sp.random(cells, genes, density=0.05, format="csr", dtype=np.float32, random_state=0)
		return ad.AnnData(X, obs=pd.DataFrame(index=[f"c{i}" for i in range(cells)]),
						  var=pd.DataFrame(index=[f"G{i}" for i in range(genes)]))

	def block_size_of(self, db_file):
		import duckdb
		conn = duckdb.connect(db_file, read_only=True)
		size = conn.execute("SELECT max(block_size) FROM pragma_database_size() WHERE block_size > 0").fetchone()[0]
		conn.close()
		return size

	def test_wide_data_block_size_default(self):
		#a matrix wider than the threshold should pick duckdb's minimum block size automatically,
		#otherwise duckdb reserves ~one 262144 byte block per column and runs out of memory
		wide = self.make_wide_adata(MakeDb.WIDE_GENE_THRESHOLD + 1000)
		if os.path.exists(self.db_file):
			os.remove(self.db_file)
		MakeDb(adata=wide, db_name=self.db_name, db_path=self.db_path, print_output=False,
			   db_config={"memory_limit": "4GB"})
		self.assertEqual(self.block_size_of(self.db_file), MakeDb.WIDE_BLOCK_SIZE)
		os.remove(self.db_file)

	def test_narrow_data_keeps_duckdb_default(self):
		#below the threshold nothing should change, so existing databases keep their current format
		narrow = self.make_wide_adata(MakeDb.WIDE_GENE_THRESHOLD - 1000)
		if os.path.exists(self.db_file):
			os.remove(self.db_file)
		MakeDb(adata=narrow, db_name=self.db_name, db_path=self.db_path, print_output=False,
			   db_config={"memory_limit": "4GB"})
		self.assertNotEqual(self.block_size_of(self.db_file), MakeDb.WIDE_BLOCK_SIZE)
		os.remove(self.db_file)

	def test_explicit_block_size_overrides_default(self):
		#an explicit block_size must win even when the matrix is wide
		wide = self.make_wide_adata(MakeDb.WIDE_GENE_THRESHOLD + 1000)
		if os.path.exists(self.db_file):
			os.remove(self.db_file)
		MakeDb(adata=wide, db_name=self.db_name, db_path=self.db_path, print_output=False,
			   block_size=262144, db_config={"memory_limit": "4GB"})
		self.assertEqual(self.block_size_of(self.db_file), 262144)
		os.remove(self.db_file)

	def test_build_database(self):
		if os.path.exists(self.db_file): #tearDown 
			os.remove(self.db_file)
		MakeDb(adata=self.adata, db_name=self.db_name, db_path=self.db_path, print_output=False)
		self.assertTrue(os.path.exists(self.db_file))

	def test_query_database(self):
		# Build a fresh database for this test (don't rely on test_build_database running first)
		if os.path.exists(self.db_file):
			os.remove(self.db_file)
		MakeDb(adata=self.adata, db_name=self.db_name, db_path=self.db_path, print_output=False)
		adata_sql = AnnSQL.AnnSQL(db=self.db_file, print_output=False)
		result = adata_sql.query("SELECT * FROM X")
		if os.path.exists(self.db_file): #tearDown
			os.remove(self.db_file)
		self.assertEqual(len(result), self.adata.shape[0])

	def test_backed_mode(self):
		import warnings
		warnings.filterwarnings('ignore')
		self.adata = sc.datasets.pbmc3k_processed()
		self.adata = sc.read_h5ad("data/pbmc3k_processed.h5ad", backed="r")
		MakeDb(adata=self.adata, db_name=self.db_name, db_path=self.db_path, print_output=False)
		adata_sql = AnnSQL.AnnSQL(db=self.db_file)
		result = adata_sql.query("SELECT * FROM X")
		if os.path.exists("data"): #tearDown here. 
			os.remove("data/pbmc3k_processed.h5ad")
			os.rmdir("data")
		if os.path.exists(self.db_file): #tearDown here. 
			os.remove(self.db_file)
		self.assertEqual(len(result), self.adata.shape[0])

	def test_backed_mode_buffer_file(self):
		import warnings
		warnings.filterwarnings('ignore')
		self.adata = sc.datasets.pbmc3k_processed()
		self.adata = sc.read_h5ad("data/pbmc3k_processed.h5ad", backed="r")
		MakeDb(adata=self.adata, db_name=self.db_name, db_path=self.db_path, chunk_size=500, make_buffer_file=True, print_output=False)
		adata_sql = AnnSQL.AnnSQL(db=self.db_file)
		result = adata_sql.query("SELECT * FROM X")
		if os.path.exists("data"): #tearDown here. 
			os.remove("data/pbmc3k_processed.h5ad")
			os.rmdir("data")
		if os.path.exists(self.db_file): #tearDown here. 
			os.remove(self.db_file)
		self.assertEqual(len(result), self.adata.shape[0])

	def test_export_X_to_csv(self):
		# Build a fresh database for this test
		if os.path.exists(self.db_file):
			os.remove(self.db_file)
		MakeDb(adata=self.adata, db_name=self.db_name, db_path=self.db_path, print_output=False)
		
		# Test the export_X_to_csv method
		adata_sql = AnnSQL.AnnSQL(db=self.db_file, print_output=False)
		test_csv_file = "test_X_export.csv"
		
		# Export X table to CSV
		adata_sql.export_X_to_csv(test_csv_file)
		
		# Verify the CSV file was created
		self.assertTrue(os.path.exists(test_csv_file))
		
		# Verify the content by reading the CSV and comparing with query result
		import pandas as pd
		csv_data = pd.read_csv(test_csv_file)
		query_result = adata_sql.query("SELECT * FROM X")
		
		# Check if the shapes match
		self.assertEqual(csv_data.shape, query_result.shape)
		
		# Cleanup
		if os.path.exists(test_csv_file):
			os.remove(test_csv_file)
		if os.path.exists(self.db_file):
			os.remove(self.db_file)

if __name__ == "__main__":
	unittest.main()