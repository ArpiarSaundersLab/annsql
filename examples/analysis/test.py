#import libraries
import scanpy as sc
from memory_profiler import memory_usage
import time
from AnnSQL import AnnSQL
from AnnSQL.MakeDb import MakeDb
import pandas as pd
import numpy as np


#open the database
asql = AnnSQL(db="../db/PCCM.asql")
asql.show_tables()

# Get top 2000 genes
genes_df = asql.query("SELECT gene_names_orig FROM var ORDER BY variance DESC")
genes = genes_df["gene_names_orig"].tolist()[:2000]

#temp buffer table to reduce memory usage
asql.query_raw("DROP TABLE IF EXISTS PC_scores_temp;")
asql.query_raw("CREATE TABLE PC_scores_temp (cell_id STRING, pc INT, partial_score DOUBLE);")

start_time = time.time()
dot_product_chunk_query_2 = """
WITH A AS (
    SELECT pc, loading
    FROM PC_loadings
),
B AS (
	SELECT
    X.cell_id,
    A.pc,
    X.Zeb1 * A.loading AS partial_score
	FROM X_standard_wide X
	CROSS JOIN A
)
SELECT cell_id, pc, SUM(partial_score) FROM B GROUP BY cell_id, pc ORDER BY cell_id, pc;
"""
asql.query(dot_product_chunk_query_2)
print("--- %s seconds dot product ---" % (time.time() - start_time))

asql.query('SELECT * FROM PC_scores')


#process genes in chunks
for i in range(0, len(genes), chunk_size):
	chunk_genes = genes[i:i+chunk_size]
	genes_clause = ", ".join([f"'{gene}'" for gene in chunk_genes])
	dot_product_chunk_query = f"""
	SELECT 
		X.cell_id,
		P.pc,
		SUM(X.value * P.loading) AS partial_score
	FROM X_standard X
	JOIN PC_loadings P ON X.gene = P.gene
	WHERE X.gene IN ({genes_clause})
	GROUP BY X.cell_id, P.pc
	ORDER BY X.cell_id, P.pc;
	"""
	#add partial results into the temp
	asql.query_raw(f"INSERT INTO PC_scores_temp {dot_product_chunk_query}")
	print(f"PCs Chunk {i} of {len(genes)}")

#pull it all together
asql.query_raw("DROP TABLE IF EXISTS PC_scores;")
asql.query_raw("""
	CREATE TABLE PC_scores AS
		SELECT cell_id, pc, SUM(partial_score) AS pc_score FROM PC_scores_temp
	GROUP BY cell_id, pc
	ORDER BY cell_id, pc;
	""")

#drop temp table
asql.query_raw("DROP TABLE PC_scores_temp;")
