import unittest
import AnnSQL as AnnSQL
import scanpy as sc

class TestQuery(unittest.TestCase):
	def setUp(self):
		self.adata = sc.datasets.pbmc68k_reduced()
		self.adata_sql = AnnSQL.AnnSQL(adata=self.adata)

	def test_select_query(self):
		result = self.adata_sql.query("SELECT * FROM X")
		self.assertEqual(len(result), self.adata.shape[0])

if __name__ == '__main__':
	unittest.main()