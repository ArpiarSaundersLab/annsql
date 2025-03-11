library(Seurat)
library(SeuratDisk)
library(anndata)

#Note: Open rstudio with sceasy conda environment

# Define the file increments
n_cells <- c("1000", "5000", "10000", "15000", "20000", "25000", "30000", "35000", "40000", "45000", "50000", "75000", "100000", "250000")

#Iterate splatter data and convert to Seurat
for (i in n_cells) {
	print(i)
	h5ad_file <- paste0("/home/kenny/Documents/OHSU/Projects/AnnSql/examples/data/splatter/data_", i, ".h5ad")
	seurat_file <- paste0("/home/kenny/Documents/OHSU/Projects/AnnSql/examples/rdata/splatter/data_", i, ".rds")
	Convert(h5ad_file, dest = "h5seurat", overwrite = TRUE)
	seurat_obj <- LoadH5Seurat(sub(".h5ad$", ".h5seurat", h5ad_file))
	saveRDS(seurat_obj, file = seurat_file)
	file.remove(sub(".h5ad$", ".h5seurat", h5ad_file))
	seurat_obj <- readRDS(seurat_file)
	s4_seurat <- UpdateSeuratObject(object = seurat_obj)
	saveRDS(s4_seurat, seurat_file)
	rm(seurat_obj)
	rm(s4_seurat)
}


# Iterate random data and convert to Seurat
for (i in n_cells) {
	print(i)
	h5ad_file <- paste0("/home/kenny/Documents/OHSU/Projects/AnnSql/examples/data/random/data_", i, ".h5ad")
	seurat_file <- paste0("/home/kenny/Documents/OHSU/Projects/AnnSql/examples/rdata/random/data_", i, ".rds")
	Convert(h5ad_file, dest = "h5seurat", overwrite = TRUE)
	seurat_obj <- LoadH5Seurat(sub(".h5ad$", ".h5seurat", h5ad_file), meta.data = FALSE, misc = FALSE)
	adata <- read_h5ad(h5ad_file, backed = TRUE)
	seurat_obj <- AddMetaData(seurat_obj, adata$obs)
	saveRDS(seurat_obj, file = seurat_file)
	file.remove(sub(".h5ad$", ".h5seurat", h5ad_file))
	seurat_obj <- readRDS(seurat_file)
	s4_seurat <- UpdateSeuratObject(object = seurat_obj)
	saveRDS(s4_seurat, seurat_file)
	rm(seurat_obj)
	rm(s4_seurat)
}
