import unittest
import AnnSQL as AnnSQL
from MakeDb import MakeDb
import scanpy as sc
import os
import time

class TestDatabase(unittest.TestCase):
	def setUp(self):
		self.adata = sc.datasets.pbmc68k_reduced()
		self.db_path = "tests/db/"
		self.db_name = "pbmc68k_reduced"
		self.db_file = os.path.join(self.db_path, f"{self.db_name}.asql")

	def test_build_database(self):
		if os.path.exists(self.db_file): #tearDown here. We need this file
			os.remove(self.db_file)
		MakeDb(adata=self.adata, db_name=self.db_name, db_path=self.db_path)
		self.assertTrue(os.path.exists(self.db_file))

	def test_query_database(self):
		adata_sql = AnnSQL.AnnSQL(db=self.db_file)
		result = adata_sql.query("SELECT * FROM X")
		if os.path.exists(self.db_file): #tearDown here. We need this file
			os.remove(self.db_file)
		self.assertEqual(len(result), self.adata.shape[0])

if __name__ == "__main__":
	unittest.main()
