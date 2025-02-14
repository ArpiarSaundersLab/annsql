library(Seurat)

#the different file increments in an array (must be in quotes R, converts to 1e5 at 100k if not in quotes)
n_cells <- c("1000", "5000", "10000", "15000", "20000", "25000", "30000", "35000", "40000", "45000", "50000", "75000", "100000")

#the random values for gene_ids generated from the python script for consistency
gene_ids <- c(2732, 9845, 3264, 4859, 9225, 7891, 4373, 5874, 6744, 3468, 705, 2599, 2222, 7768, 2897, 9893, 537, 6216, 6921, 6036)
gene_ids_2 <- c(2163, 5072, 4851, 7877, 2046, 1871, 7599, 2496, 8291, 755, 797, 659, 3219, 8615, 7456, 3337, 2745, 4735, 8736, 6687)

#function to confirm the result already exists. This helps avoid running same code twice.
check_combination <- function(gene_id, size, method, filter) {
	existing_data <- read.csv("/home/kenny/Documents/OHSU/Projects/AnnSql/examples/results/comparisons_splatter_R.csv")
	combination_exists <- any(existing_data$gene_id == gene_id & existing_data$size == size & existing_data$method == method & existing_data$filter == filter)
	return(combination_exists)
}

for (i in n_cells){

	#open the seurat_obj file
	seurat_obj <- readRDS( paste0("/home/kenny/Documents/OHSU/Projects/AnnSql/examples/rdata/splatter/data_", i, ".rds"))


	#iterate through the array
	for (gene_id in gene_ids) {

		print(paste("File", i, "Gene", gene_id))

		#array with size,method,runtime,filter,runtime_log,runtime_log10,size_log,size_log10 columns
		df <- data.frame(matrix(ncol = 9))
		colnames(df) <- c("size", "method", "runtime", "filter", "gene_id", "runtime_log", "runtime_log10", "size_log", "size_log10")
		
		#if combination returnss false, then run the code
		if (!check_combination(gene_id, i, "Seurat", "Filter1")) {
			print("Running Filter 1")
			elapsed_time_1 <- system.time({
				gene_column <- paste0("gene-", gene_id) 
				filtered_cells <- colnames(seurat_obj)[seurat_obj@assays$RNA@data[gene_column, ] > 0.5]
				filtered_seurat <- subset(seurat_obj, cells = filtered_cells, return.null=TRUE)
			})
			#insert into the dataframe
			df <- rbind(df, c(i, "Seurat", as.numeric(elapsed_time_1["elapsed"][1]), "Filter1", gene_id, log(as.numeric(elapsed_time_1["elapsed"][1])), log10(as.numeric(elapsed_time_1["elapsed"][1])), log(as.numeric(i)), log10(as.numeric(i))))
		}

		#Filter 2
		if (!check_combination(gene_id, i, "Seurat", "Filter2")) {
			print("Running Filter 2")
			elapsed_time_2 <- system.time({
				gene_column <- paste0("gene-", gene_id) 
				gene_column_2 <- paste0("gene-", sample(gene_ids_2, 1))
				filtered_cells <- colnames(seurat_obj)[seurat_obj@assays$RNA@data[gene_column, ] > 0.5 & seurat_obj@assays$RNA@data[gene_column_2, ] < 0.5]
				filtered_seurat <- subset(seurat_obj, cells = filtered_cells, return.null=TRUE)
			})
			#insert into the dataframe
			df <- rbind(df, c(i, "Seurat", as.numeric(elapsed_time_2["elapsed"][1]), "Filter2", gene_id, log(as.numeric(elapsed_time_2["elapsed"][1])), log10(as.numeric(elapsed_time_2["elapsed"][1])), log(as.numeric(i)), log10(as.numeric(i))))
		}


		#Filter 3
		if (!check_combination(gene_id, i, "Seurat", "Filter3")) {
			print("Running Filter 3")
			elapsed_time_3 <- system.time({
				gene_column <- paste0("gene-", gene_id) 
				gene_column_2 <- paste0("gene-", sample(gene_ids_2, 1))
				gene_column_3 <- paste0("gene-", sample(gene_ids_2, 1))
				filtered_cells <- colnames(seurat_obj)[seurat_obj@assays$RNA@data[gene_column, ] > 0.5 & seurat_obj@assays$RNA@data[gene_column_2, ] < 0.5 & seurat_obj@assays$RNA@data[gene_column_3, ] > 0.5]
				filtered_seurat <- subset(seurat_obj, cells = filtered_cells, return.null=TRUE)
			})
			#insert into the dataframe
			df <- rbind(df, c(i, "Seurat", as.numeric(elapsed_time_3["elapsed"][1]), "Filter3", gene_id, log(as.numeric(elapsed_time_3["elapsed"][1])), log10(as.numeric(elapsed_time_3["elapsed"][1])), log(as.numeric(i)), log10(as.numeric(i))))
		}


		#Filter 4
		if (!check_combination(gene_id, i, "Seurat", "Filter4")) {
			print("Running Filter 4")
			elapsed_time_4 <- system.time({
				gene_column <- paste0("gene-", gene_id) 
				gene_column_2 <- paste0("gene-", sample(gene_ids_2, 1))
				gene_column_3 <- paste0("gene-", sample(gene_ids_2, 1))
				gene_column_4 <- paste0("gene-", sample(gene_ids_2, 1))
				filtered_cells <- colnames(seurat_obj)[seurat_obj@assays$RNA@data[gene_column, ] > 0.5 & seurat_obj@assays$RNA@data[gene_column_2, ] < 0.5 & seurat_obj@assays$RNA@data[gene_column_3, ] > 0.5 & seurat_obj@assays$RNA@data[gene_column_4, ] < 0.5]
				filtered_seurat <- subset(seurat_obj, cells = filtered_cells, return.null=TRUE)
			})
			#insert into the dataframe
			df <- rbind(df, c(i, "Seurat", as.numeric(elapsed_time_4["elapsed"][1]), "Filter4", gene_id, log(as.numeric(elapsed_time_4["elapsed"][1])), log10(as.numeric(elapsed_time_4["elapsed"][1])), log(as.numeric(i)), log10(as.numeric(i))))
		}

		#Filter 5
		if (!check_combination(gene_id, i, "Seurat", "Filter5")) {
			print("Running Filter 5")
			elapsed_time_5 <- system.time({
				gene_column <- paste0("gene-", gene_id) 
				gene_column_2 <- paste0("gene-", sample(gene_ids_2, 1))
				#seurat_obj$gene_3_div_gene_4 <- seurat_obj[["RNA"]]@data[gene_column, ] / seurat_obj[["RNA"]]@data[gene_column_2, ]
				filtered_cells <- colnames(seurat_obj)[seurat_obj@assays$RNA@data[gene_column, ] > 0 & seurat_obj[["RNA"]]@data[gene_column, ] / seurat_obj[["RNA"]]@data[gene_column_2, ] < 0.5]
				if (length(filtered_cells) > 0) {
					filtered_seurat <- subset(seurat_obj, cells = filtered_cells, return.null = TRUE)
				} else {
					print("No cells left after filtering.")
				}								
			})
			#insert into the dataframe
			df <- rbind(df, c(i, "Seurat", as.numeric(elapsed_time_5["elapsed"][1]), "Filter5", gene_id, log(as.numeric(elapsed_time_5["elapsed"][1])), log10(as.numeric(elapsed_time_5["elapsed"][1])), log(as.numeric(i)), log10(as.numeric(i))))
		}


		#Filter 6
		if (!check_combination(gene_id, i, "Seurat", "Filter6")) {
			print("Running Filter 6")
			elapsed_time_6 <- system.time({
				cell_types <- seurat_obj$cell_type
				gene_column <- paste0("gene-", gene_id) 
				avg_gene_1 <- aggregate(seurat_obj[["RNA"]]@data[gene_column, ] ~ cell_types, FUN = mean)
				colnames(avg_gene_1) <- c("cell_type", "avg_gene_1")
			})
			#insert into the dataframe
			df <- rbind(df, c(i, "Seurat", as.numeric(elapsed_time_6["elapsed"][1]), "Filter6", gene_id, log(as.numeric(elapsed_time_6["elapsed"][1])), log10(as.numeric(elapsed_time_6["elapsed"][1])), log(as.numeric(i)), log10(as.numeric(i))))
		}

		#adjust dataframe (removes top row)
		df <- df[-1, ]

		#append the dataframe to the csv file
		write.table(df, file = paste0("/home/kenny/Documents/OHSU/Projects/AnnSql/examples/results/comparisons_splatter_R.csv"), 
						sep = ",", 
						col.names = !file.exists(paste0("/home/kenny/Documents/OHSU/Projects/AnnSql/examples/results/comparisons_splatter_R.csv")), 
						row.names = FALSE, 
						append = TRUE,
						quote = FALSE
						)
	}

}