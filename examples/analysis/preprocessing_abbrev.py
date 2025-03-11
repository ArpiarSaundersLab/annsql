#preprocessing, condensed version.
import os
import scanpy as sc
from AnnSQL import AnnSQL
from AnnSQL.MakeDb import MakeDb

#make the db
filepath = "../data/PCCM.h5ad"
adata = sc.read_h5ad(filepath)
adata = adata.raw.to_adata()
if os.path.exists("db/PCCM2.asql"):
	os.remove("db/PCCM2.asql")
MakeDb(adata=adata, db_name="PCCM2", db_path="db")

#preprocess the data
asql = AnnSQL(db="db/PCCM2.asql")
asql.save_raw()
asql.calculate_total_counts(chunk_size=750)
asql.calculate_gene_counts(chunk_size=750) 
asql.filter_by_cell_counts(min_cell_count=2500, max_cell_count=40000)
asql.filter_by_gene_counts(min_gene_counts=100, max_gene_counts=10000)
asql.expression_normalize(total_counts_per_cell=10000, chunk_size=750)
asql.expression_log(log_type="LN", chunk_size=750)
asql.calculate_variable_genes(chunk_size=750, save_var_names=False)
asql.save_highly_variable_genes(top_variable_genes=2500)
asql.calculate_pca(n_pcs=50, 
					top_variable_genes=2500, 
					chunk_size=750, 
					max_cells_memory_threshold=5000, 
					print_progress=True)
asql.calculate_umap()
asql.calculate_leiden_clusters(resolution=1.0, n_neighbors=5)
asql.plot_umap()
asql.plot_umap(color_by="leiden_clusters", annotate=True)
asql.plot_umap(color_by="Gad2")
asql.plot_umap(color_by="Gad1", counts_table="X_raw")