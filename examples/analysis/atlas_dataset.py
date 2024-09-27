#import libraries
import scanpy as sc
from MakeDb import MakeDb
from AnnSQL import AnnSQL
import time
import pandas as pd

# #same as above, but using scanpy's read function
#adata = sc.read_h5ad("/home/kenny/Documents/OHSU/Projects/TAP/data/celltypist_models/chunked_approach/Macosko_Mouse_Atlas_Single_Nuclei.Use_Backed.h5ad", backed="r")

# #build a database to query later
#MakeDb(adata=adata, 
# 		db_name="Macosko_Mouse_Atlas_4mil", 
# 		db_path="../db/", 
# 		create_basic_indexes=True,
# 		layers=["X", "obs","var"])

#call the AnnSQL class
adata_sql = AnnSQL(db="../db/Macosko_Mouse_Atlas.asql")

###################################################
# attempt to normalize the data
###################################################
#get the gene names
gene_names = adata_sql.query(f"Describe X")['column_name'][1:].values

#removes total_counts from the array & sets to 0
if "total_counts" in gene_names:
	adata_sql.update_query(f"UPDATE X SET total_counts = 0;")
	gene_names = gene_names[:-1] 
else:
	adata_sql.query(f"ALTER TABLE X ADD COLUMN total_counts FLOAT DEFAULT 0;")

#iterates gene_names in chunks
start_time = time.time()
chunk_size = 990 #Ddb limits arguments to 1k
for i in range(0, len(gene_names), chunk_size):
	chunk = gene_names[i:i+chunk_size]
	chunk = " + ".join(chunk) + " + total_counts"
	adata_sql.update_query(f"UPDATE X SET total_counts = ({chunk});")
end_time = time.time()
print("Total Counts Time: ", end_time-start_time)

#normalize to 10k and log2
start_time = time.time()
chunk_size = 200  #reduces db memory usage
for i in range(0, len(gene_names), chunk_size):
	updates = []
	interval_time = time.time()
	chunk = gene_names[i:i + chunk_size]
	for gene in chunk:
		updates.append(f"{gene} = LOG2(({gene} / total_counts) * 1e4 + 1e-5)")
	update_query = f"UPDATE X SET {', '.join(updates)}"
	adata_sql.update_query(update_query)
	print(f"Processed chunk {i // chunk_size + 1} in {time.time() - interval_time} seconds")
end_time = time.time()
print("Normalize & Log2 Time: ", end_time - start_time)



###################################################
#get the top highly variable genes
###################################################
adata.var.to_dict()
df = pd.DataFrame(adata.var.to_dict(), columns=["esn_name","gene_name"], index=adata.var_names)
df["esn_name"] = df["index"]
df = df.drop(columns=["index"])
df = df[["gene_name", "esn_name"]]
df["variance"] = 0
df

adata_sql.open_db()
adata_sql.conn.register("df", df)
adata_sql.conn.execute("DROP TABLE IF EXISTS var")
adata_sql.conn.execute(f"CREATE TABLE var AS SELECT * FROM df")
adata_sql.query("SELECT * FROM var")

adata_sql.query(f"SELECT VARIANCE(ENSMUSG00000051951) FROM X;")
adata_sql.query_raw(f"UPDATE var SET ENSMUSG00000051951 = VARIANCE(SELECT ENSMUSG00000051951 FROM X)")

adata_sql.query_raw(f"""
    WITH variance_computation AS (
        SELECT VARIANCE(ENSMUSG00000051951) AS computed_variance
        FROM X
    )
    UPDATE var
    SET ENSMUSG00000051951 = (SELECT computed_variance FROM variance_computation);
""")