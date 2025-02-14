library(Seurat)

#the different file increments in an array (must be in quotes R, converts to 1e5 at 100k if not in quotes)
n_cells <- c("1000", "5000", "10000", "15000", "20000", "25000", "30000", "35000", "40000", "45000", "50000", "75000", "100000")
files <- c(1,2,3,4,5,6)

for (file in files){


	#iterate through the array
	for (i in n_cells) {

		print(paste("File", file))
		print(paste("Starting", i))

		#open the seurat_obj file
		seurat_obj <- readRDS( paste0("/home/kenny/Documents/OHSU/Projects/AnnSql/examples/rdata/random/data_", i, ".rds"))

		#an array with size,method,runtime,filter,runtime_log,runtime_log10,size_log,size_log10 columns
		df <- data.frame(matrix(ncol = 8))
		colnames(df) <- c("size", "method", "runtime", "filter", "runtime_log", "runtime_log10", "size_log", "size_log10")

		#filter 1
		elapsed_time_1 <- system.time({
		filtered_seurat = (subset(seurat_obj, cells = colnames(subset(seurat_obj, subset = `gene-1` > 0.5, return.null=TRUE))))
		})

		#insert into the dataframe
		df <- rbind(df, c(i, "Seurat", as.numeric(elapsed_time_1["elapsed"][1]), "Filter1", log(as.numeric(elapsed_time_1["elapsed"][1])), log10(as.numeric(elapsed_time_1["elapsed"][1])), log(1000), log10(1000)))

		#filter 2
		elapsed_time_2 <- system.time({
		filtered_seurat = (subset(seurat_obj, cells = colnames(subset(seurat_obj, subset = `gene-1` > 0.5 & `gene-2` < 0.5, return.null=TRUE))))
		})

		#insert into the dataframe
		df <- rbind(df, c(i, "Seurat", as.numeric(elapsed_time_2["elapsed"][1]), "Filter2", log(as.numeric(elapsed_time_2["elapsed"][1])), log10(as.numeric(elapsed_time_2["elapsed"][1])), log(1000), log10(1000)))

		#filter 3
		elapsed_time_3 <- system.time({
		filtered_seurat = (subset(seurat_obj, cells = colnames(subset(seurat_obj, subset = `gene-1` > 0.5 & `gene-2` < 0.5 & `gene-3` > 0.5, return.null=TRUE))))
		})

		#insert into the dataframe
		df <- rbind(df, c(i, "Seurat", as.numeric(elapsed_time_3["elapsed"][1]), "Filter3", log(as.numeric(elapsed_time_3["elapsed"][1])), log10(as.numeric(elapsed_time_3["elapsed"][1])), log(1000), log10(1000)))

		#filter 4
		elapsed_time_4 <- system.time({
		filtered_seurat = (subset(seurat_obj, cells = colnames(subset(seurat_obj, subset = `gene-1` > 0.5 & `gene-2` < 0.5 & `gene-3` > 0.5 & `gene-4` < 0.5, return.null=TRUE))))
		})

		#insert into the dataframe
		df <- rbind(df, c(i, "Seurat", as.numeric(elapsed_time_4["elapsed"][1]), "Filter4", log(as.numeric(elapsed_time_4["elapsed"][1])), log10(as.numeric(elapsed_time_4["elapsed"][1])), log(1000), log10(1000)))

		#filter 5
		elapsed_time_5 <- system.time({
		seurat_obj$gene_3_div_gene_4 <- seurat_obj[["RNA"]]@data["gene-3", ] / seurat_obj[["RNA"]]@data["gene-4", ]
		filtered_seurat <- subset(seurat_obj, subset = (`gene-1` > 0.5 & `gene-2` < 0.3) & (gene_3_div_gene_4 > 0 & gene_3_div_gene_4 < 0.5))
		})

		#insert into the dataframe
		df <- rbind(df, c(i, "Seurat", as.numeric(elapsed_time_5["elapsed"][1]), "Filter5", log(as.numeric(elapsed_time_5["elapsed"][1])), log10(as.numeric(elapsed_time_5["elapsed"][1])), log(1000), log10(1000)))

		#filter 6
		elapsed_time_6 <- system.time({
		cell_types <- seurat_obj$cell_type
		avg_gene_1 <- aggregate(seurat_obj[["RNA"]]@data["gene-1", ] ~ cell_types, FUN = mean)
		colnames(avg_gene_1) <- c("cell_type", "avg_gene_1")
		})

		#insert into the dataframe
		df <- rbind(df, c(i, "Seurat", as.numeric(elapsed_time_6["elapsed"][1]), "Filter6", log(as.numeric(elapsed_time_6["elapsed"][1])), log10(as.numeric(elapsed_time_6["elapsed"][1])), log(1000), log10(1000)))

		# print the dataframe
		df <- df[-1, ]

		#append the dataframe to the csv file
		write.table(df, file = paste0("/home/kenny/Documents/OHSU/Projects/AnnSql/examples/results/comparisons_random_",file,".csv"), 
						sep = ",", 
						col.names = !file.exists(paste0("/home/kenny/Documents/OHSU/Projects/AnnSql/examples/results/comparisons_random_",file,".csv")), 
						row.names = FALSE, 
						append = TRUE,
						quote = FALSE
						)
	}

}