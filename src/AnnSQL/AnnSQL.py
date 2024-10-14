from BuildDb import BuildDb
import scanpy as sc
import pandas as pd
import duckdb
import warnings
import logging
import os 

class AnnSQL:
	def __init__(self, adata=None, db=None, create_all_indexes=False, create_basic_indexes=False, layers=["X", "obs", "var", "var_names", "obsm", "varm", "obsp", "uns"]):
		self.adata = self.open_anndata(adata)
		self.db = db
		self.create_basic_indexes = create_basic_indexes
		self.create_all_indexes = create_all_indexes
		self.layers = layers
		self.validate_params()
		self.is_open = False
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

	def open_anndata(self,adata):
		if not isinstance(adata, sc.AnnData) and isinstance(adata, str):	
			return sc.read_h5ad(adata)
		else:
			return adata
			
	def open_db(self):
		if self.db is not None:
			self.conn = duckdb.connect(self.db)
			self.is_open = True

	def close_db(self):
		if self.db is not None:
			self.conn.close()
			self.is_open = False

	def build_db(self):
		self.conn = duckdb.connect(':memory:')
		db = BuildDb(adata=self.adata, conn=self.conn, create_all_indexes=self.create_all_indexes, create_basic_indexes=self.create_basic_indexes, layers=self.layers)
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
			if self.db is not None and self.adata is None:
				print('Warning: No adata object provided. return_type="adata" is disabled.')
				return result_df
			return self.adata[result_df["cell_id"]]

	def query_raw(self, query):
		self.open_db()
		result = self.conn.execute(query)
		self.close_db()
		return result

	def update_query(self, query, suppress_message=False):
		if 'SELECT' in query.upper() or 'DELETE' in query.upper():
			raise ValueError('SELECT detected. Please use query() instead')
		try:
			self.open_db()
			self.conn.execute(query)
			self.close_db()
			if suppress_message == False:
				print("Query Successful")
		except Exception as e:
			print("Update Query Error:", e)

	def delete_query(self, query, suppress_message=False):
		if 'SELECT' in query.upper() or 'UPDATE' in query.upper():
			raise ValueError('SELECT detected. Please use query() instead')
		try:
			self.open_db()
			self.conn.execute(query)
			self.close_db()
			if suppress_message == False:
				print("Delete Query Successful")
		except Exception as e:
			print("Delete Query Error:", e)

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
	
	def replace_special_chars(self, string):
		return string.replace("-", "_").replace(".", "_")

	def expression_normalize(self, total_counts_per_cell=1e4, chunk_size=200, print_progress=False):
		self.check_chunk_size(chunk_size)
		if 'total_counts' not in self.query("SELECT * FROM X LIMIT 1").columns:
			print("Total counts not found...")
			self.calculate_total_counts(chunk_size=chunk_size,print_progress=print_progress)
		
		print("Expression Normalization Started")
		gene_names = self.query(f"Describe X")['column_name'][1:].values
		if 'total_counts' in gene_names:
			gene_names = gene_names[:-1]
		for i in range(0, len(gene_names), chunk_size):
			updates = []
			chunk = gene_names[i:i + chunk_size]
			for gene in chunk:
				if gene == 'total_counts':
					continue
				updates.append(f"{gene} = (({gene} / total_counts) * {total_counts_per_cell})")
			update_query = f"UPDATE X SET {', '.join(updates)}"
			self.update_query(update_query, suppress_message=True)
			if print_progress == True:
				print(f"Processed chunk {i // chunk_size + 1}")
		print("Expression Normalization Complete")

	def expression_log(self, log_type="LN", chunk_size=200, print_progress=False):
		#log_type can be LN, LOG (LOG2 alias), LOG2, LOG10
		self.check_chunk_size(chunk_size)
		gene_names = self.query(f"Describe X")['column_name'][1:].values
		if 'total_counts' in gene_names:
			gene_names = gene_names[:-1]
		
		print("Log Transform Started")
		for i in range(0, len(gene_names), chunk_size):
			updates = []
			chunk = gene_names[i:i + chunk_size]
			for gene in chunk:
				if gene == 'total_counts':
					continue
				updates.append(f"{gene} = {log_type}({gene}+1)")
				#handle zero values
				#updates.append(f"{gene} = CASE WHEN {gene} = 0 OR {gene} = 0.0 THEN 0.0 ELSE {log_type}({gene}) END")
			update_query = f"UPDATE X SET {', '.join(updates)}"
			self.update_query(update_query, suppress_message=True)
			if print_progress == True:
				print(f"Processed chunk {i // chunk_size + 1}")
		print("Log Transform Complete")

	
	def calculate_total_counts(self, chunk_size=250, print_progress=False):
		self.check_chunk_size(chunk_size)
		gene_names = self.query(f"Describe X")['column_name'][1:].values
		
		if "total_counts" in gene_names:
			self.update_query(f"UPDATE X SET total_counts = 0;")
			gene_names = gene_names[:-1] 
		else:
			self.query(f"ALTER TABLE X ADD COLUMN total_counts FLOAT DEFAULT 0;")
		
		print("Total Counts Calculation Started")
		for i in range(0, len(gene_names), chunk_size):
			chunk = gene_names[i:i+chunk_size]
			chunk = " + ".join(chunk) + " + total_counts"
			self.update_query(f"UPDATE X SET total_counts = ({chunk});", suppress_message=True)
			if print_progress == True:
				print(f"Processed chunk {i // chunk_size + 1}")

		#set obs total_counts
		if 'total_counts' not in self.query("SELECT * FROM obs LIMIT 1").columns:
			self.query_raw("ALTER TABLE obs ADD COLUMN total_counts FLOAT DEFAULT 0;")
		self.query_raw("UPDATE obs SET total_counts = (SELECT total_counts FROM X WHERE obs.cell_id = X.cell_id)")
		print("Total Counts Calculation Complete")

	def calculate_gene_counts(self, chunk_size=200, print_progress=False, gene_field="gene_names"):
		self.check_chunk_size(chunk_size)
		gene_names_df = self.query(f"SELECT {gene_field} FROM var")
		gene_names_df["gene_counts"] = 0.0		
		gene_names_df = gene_names_df.reset_index(drop=True)
		
		var_table = self.query("SELECT * FROM var LIMIT 1")
		if var_table.shape[0] == 0:
			print("Creating Var Table")
			self.open_db()
			self.conn.register("gene_names_df", gene_names_df)
			self.conn.execute(f"CREATE TABLE var AS SELECT * FROM gene_names_df")
		else:
			print("Updating Var Table")
			if "gene_counts" not in var_table.columns:
				self.update_query("ALTER TABLE var ADD COLUMN gene_counts FLOAT DEFAULT 0;", suppress_message=True)
			else:
				self.update_query("UPDATE var SET gene_counts = 0.0;", suppress_message=True)	

		print("Gene Counts Calculation Started")
		gene_counts = []
		for i in range(0, len(gene_names_df), chunk_size):
			chunk = gene_names_df[gene_field][i:i+chunk_size]
			query = f"SELECT {', '.join([f'SUM({gene}) as {gene}' for gene in chunk])} FROM X;"
			counts_chunk = self.query(query)
			gene_counts.extend(counts_chunk.values.flatten())
			if print_progress == True:
				print(f"Processed chunk {i // chunk_size + 1}")

		#insert these values into the var table matching on the index.
		gene_counts_df = pd.DataFrame({"gene_counts": gene_counts})
		gene_counts_df[gene_field] = gene_names_df[gene_field]

		#update the var table with the gene_counts values
		self.open_db()
		self.conn.execute("DROP TABLE IF EXISTS gene_counts_df")
		self.conn.register("gene_counts_df", gene_counts_df)
		self.conn.execute(f"CREATE TABLE gene_counts_df AS SELECT * FROM gene_counts_df")
		self.conn.execute(f"UPDATE var SET gene_counts = (SELECT gene_counts FROM gene_counts_df WHERE var.{gene_field} = gene_counts_df.{gene_field})")
		self.conn.execute("DROP VIEW IF EXISTS gene_counts_df")
		print("Gene Counts Calculation Complete")

	def calculate_variable_genes(self, chunk_size=100, print_progress=False, gene_field="gene_names"):
		self.check_chunk_size(chunk_size)
		gene_names_df = self.query(f"SELECT {gene_field} FROM var")
		gene_names_df["variance"] = 0.0		
		gene_names_df = gene_names_df.reset_index(drop=True)
		
		var_table = self.query("SELECT * FROM var LIMIT 1")
		if var_table.shape[0] == 0:
			print("Creating Var Table")
			self.open_db()
			self.conn.register("gene_names_df", gene_names_df)
			self.conn.execute(f"CREATE TABLE var AS SELECT * FROM gene_names_df")
		else:
			print("Updating Var Table")
			if "variance" not in var_table.columns:
				self.update_query("ALTER TABLE var ADD COLUMN variance FLOAT DEFAULT 0;", suppress_message=True)
			else:
				self.update_query("UPDATE var SET variance = 0.0;", suppress_message=True)	

		variance_values = []
		for i in range(0, len(gene_names_df), chunk_size):
			chunk = gene_names_df[gene_field][i:i+chunk_size]
			query = f"SELECT {', '.join([f'VARIANCE({gene}) as {gene}' for gene in chunk])} FROM X;"
			variance_chunk = self.query(query)
			variance_values.extend(variance_chunk.values.flatten())
			if print_progress == True:
				print(f"Processed chunk {i // chunk_size + 1}")

		#insert these values into the var table matching on the index.
		variance_df = pd.DataFrame({"variance": variance_values})
		variance_df[gene_field] = gene_names_df[gene_field]

		#update the var table with the variance values
		self.open_db()
		self.conn.execute("DROP TABLE IF EXISTS variance_df")
		self.conn.register("variance_df", variance_df)
		self.conn.execute(f"CREATE TABLE variance_df AS SELECT * FROM variance_df")
		self.conn.execute(f"UPDATE var SET variance = (SELECT variance FROM variance_df WHERE var.{gene_field} = variance_df.{gene_field})")
		self.conn.execute("DROP VIEW IF EXISTS variance_df")
		print("Variance Calculation Complete")




	def check_chunk_size(self, chunk_size):
		if chunk_size > 999:
			raise ValueError('chunk_size must be less than 1000. DuckDb limitation')