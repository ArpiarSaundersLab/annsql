library(Seurat)
library(profmem)

profile_memory_usage <- function(code_block) {
  mem_profile <- profmem(code_block)
  total_mem_used_bytes <- sum(mem_profile$bytes, na.rm = TRUE)
  total_mem_used_mb <- total_mem_used_bytes / (1024^2)
  total_mem_used_mb <- round(total_mem_used_mb, digits = 2)
  return(total_mem_used_mb)
}

#the different file increments in an array (must be in quotes R, converts to 1e5 at 100k if not in quotes)
n_cells <- c("1000", "5000", "10000", "15000", "20000", "25000", "30000", "35000", "40000", "45000", "50000", "75000", "100000")

#the random values for gene_ids generated from the python script for consistency
gene_ids <- c(2732, 9845, 3264, 4859, 9225, 7891, 4373, 5874, 6744, 3468, 705, 2599, 2222, 7768, 2897, 9893, 537, 6216, 6921, 6036)
gene_ids_2 <- c(2163, 5072, 4851, 7877, 2046, 1871, 7599, 2496, 8291, 755, 797, 659, 3219, 8615, 7456, 3337, 2745, 4735, 8736, 6687)

#function to confirm the result already exists. This helps avoid running same code twice.
check_combination <- function(gene_id, size, method, filter) {
	existing_data <- read.csv("/home/kenny/Documents/OHSU/Projects/AnnSql/examples/results/comparisons_memory_splatter_R.csv")
	combination_exists <- any(existing_data$gene_id == gene_id & existing_data$size == size & existing_data$method == method & existing_data$filter_type == filter)
	return(combination_exists)
}

#iterate through the cell sizes
for (i in n_cells) {

	#iterate through the rand generated genes
	for (gene_id in gene_ids){

		print(paste("File", i, "Gene", gene_id))

		#array with the same columns in the csv
		df <- data.frame(matrix(ncol = 9))
		colnames(df) <- c("size", "method", "memory", "filter_type", "gene_id", "memory_log", "memory_log10", "size_log", "size_log10")
		
		#open the seurat_obj
		total_loading_mb <- profile_memory_usage({
			seurat_obj <- readRDS(paste0("/home/kenny/Documents/OHSU/Projects/AnnSql/examples/rdata/splatter/data_", i, ".rds"))
		})

   		print(paste("Total loading memory:", total_loading_mb, "MB"))

		if (!check_combination(gene_id, i, "Seurat", "Filter1")) {
			memory_filter_1 <- profile_memory_usage({
				print("Running Filter 1")
				gene_column <- paste0("gene-", gene_id) 
				filtered_cells <- colnames(seurat_obj)[seurat_obj@assays$RNA@data[gene_column, ] > 0.5]
				filtered_seurat <- subset(seurat_obj, cells = filtered_cells, return.null=TRUE)
			})

			total_filter_1_mb <- total_loading_mb + memory_filter_1
			print(paste("Memory of filter 1 Only:", memory_filter_1, "MB"))
			print(paste("Memory of filter 1 Total:", total_filter_1_mb, "MB"))
			df <- rbind(df, c(i, "Seurat", total_filter_1_mb, "Filter1", gene_id, log(total_filter_1_mb), log10(total_filter_1_mb), log(as.numeric(i)), log10(as.numeric(i))))
		}


		if (!check_combination(gene_id, i, "Seurat", "Filter2")) {
			memory_filter_2 <- profile_memory_usage({
				print("Running Filter 2")
				gene_column <- paste0("gene-", gene_id)
				gene_column_2 <- paste0("gene-", sample(gene_ids_2, 1))
				filtered_cells <- colnames(seurat_obj)[seurat_obj@assays$RNA@data[gene_column, ] > 0.5 & seurat_obj@assays$RNA@data[gene_column_2, ] < 0.5]
				filtered_seurat <- subset(seurat_obj, cells = filtered_cells, return.null=TRUE)
			})

			total_filter_2_mb <- total_loading_mb + memory_filter_2
			print(paste("Memory of filter 2 Only:", memory_filter_2, "MB"))
			print(paste("Memory of filter 2 Total:", total_filter_2_mb, "MB"))
			df <- rbind(df, c(i, "Seurat", total_filter_2_mb, "Filter2", gene_id, log(total_filter_2_mb), log10(total_filter_2_mb), log(as.numeric(i)), log10(as.numeric(i))))
		}

		if (!check_combination(gene_id, i, "Seurat", "Filter3")) {
			memory_filter_3 <- profile_memory_usage({
				print("Running Filter 3")
				gene_column <- paste0("gene-", gene_id)
				gene_column_2 <- paste0("gene-", sample(gene_ids_2, 1))
				gene_column_3 <- paste0("gene-", sample(gene_ids_2, 1))
				filtered_cells <- colnames(seurat_obj)[seurat_obj@assays$RNA@data[gene_column, ] > 0.5 & seurat_obj@assays$RNA@data[gene_column_2, ] < 0.5 & seurat_obj@assays$RNA@data[gene_column_3, ] > 0.5]
				filtered_seurat <- subset(seurat_obj, cells = filtered_cells, return.null=TRUE)
			})

			total_filter_3_mb <- total_loading_mb + memory_filter_3
			print(paste("Memory of filter 3 Only:", memory_filter_3, "MB"))
			print(paste("Memory of filter 3 Total:", total_filter_3_mb, "MB"))
			df <- rbind(df, c(i, "Seurat", total_filter_3_mb, "Filter3", gene_id, log(total_filter_3_mb), log10(total_filter_3_mb), log(as.numeric(i)), log10(as.numeric(i))))
		}

		if (!check_combination(gene_id, i, "Seurat", "Filter4")) {
			memory_filter_4 <- profile_memory_usage({
				print("Running Filter 4")
				gene_column <- paste0("gene-", gene_id)
				gene_column_2 <- paste0("gene-", sample(gene_ids_2, 1))
				gene_column_3 <- paste0("gene-", sample(gene_ids_2, 1))
				gene_column_4 <- paste0("gene-", sample(gene_ids_2, 1))
				filtered_cells <- colnames(seurat_obj)[seurat_obj@assays$RNA@data[gene_column, ] > 0.5 & seurat_obj@assays$RNA@data[gene_column_2, ] < 0.5 & seurat_obj@assays$RNA@data[gene_column_3, ] > 0.5 & seurat_obj@assays$RNA@data[gene_column_4, ] < 0.5]
				filtered_seurat <- subset(seurat_obj, cells = filtered_cells, return.null=TRUE)
			})

			total_filter_4_mb <- total_loading_mb + memory_filter_4
			print(paste("Memory of filter 4 Only:", memory_filter_4, "MB"))
			print(paste("Memory of filter 4 Total:", total_filter_4_mb, "MB"))
			df <- rbind(df, c(i, "Seurat", total_filter_4_mb, "Filter4", gene_id, log(total_filter_4_mb), log10(total_filter_4_mb), log(as.numeric(i)), log10(as.numeric(i))))
		}


		if (!check_combination(gene_id, i, "Seurat", "Filter5")) {
			memory_filter_5 <- profile_memory_usage({
				print("Running Filter 5")
				gene_column <- paste0("gene-", gene_id)
				gene_column_2 <- paste0("gene-", sample(gene_ids_2, 1))				
				seurat_obj$gene_3_div_gene_4 <- seurat_obj[["RNA"]]@data[gene_column_3, ] / seurat_obj[["RNA"]]@data[gene_column_4, ]
				filtered_cells <- colnames(seurat_obj)[seurat_obj@assays$RNA@data[gene_column, ] > 0.5 & seurat_obj@assays$RNA@data[gene_column_2, ] < 0.3 & seurat_obj@assays$RNA@data[gene_3_div_gene_4, ] > 0 & seurat_obj@assays$RNA@data[gene_3_div_gene_4, ] < 0.5]
				filtered_seurat <- subset(seurat_obj, cells = filtered_cells, return.null=TRUE)
			})

			total_filter_5_mb <- total_loading_mb + memory_filter_5
			print(paste("Memory of filter 5 Only:", memory_filter_5, "MB"))
			print(paste("Memory of filter 5 Total:", total_filter_5_mb, "MB"))
			df <- rbind(df, c(i, "Seurat", total_filter_5_mb, "Filter5", gene_id, log(total_filter_5_mb), log10(total_filter_5_mb), log(as.numeric(i)), log10(as.numeric(i))))

		}

		

		# same as above but for filter 6
		if (!check_combination(gene_id, i, "Seurat", "Filter6")) {
			memory_filter_6 <- profile_memory_usage({
				print("Running Filter 6")
				gene_column <- paste0("gene-", gene_id)
				cell_types <- seurat_obj$cell_type
				avg_gene_1 <- aggregate(seurat_obj[["RNA"]]@data[gene_column, ] ~ cell_types, FUN = mean)
				colnames(avg_gene_1) <- c("cell_type", "avg_gene_1")
			})

			total_filter_6_mb <- total_loading_mb + memory_filter_6
			print(paste("Memory of filter 6 Only:", memory_filter_6, "MB"))
			print(paste("Memory of filter 6 Total:", total_filter_6_mb, "MB"))
			df <- rbind(df, c(i, "Seurat", total_filter_6_mb, "Filter6", gene_id, log(total_filter_6_mb), log10(total_filter_6_mb), log(as.numeric(i)), log10(as.numeric(i))))
		}
		
		# print the dataframe
		df <- df[-1, ]

		#append the dataframe to the csv file
		write.table(df, file = paste0("/home/kenny/Documents/OHSU/Projects/AnnSql/examples/results/comparisons_memory_splatter_R.csv"), 
						sep = ",", 
						col.names = !file.exists(paste0("/home/kenny/Documents/OHSU/Projects/AnnSql/examples/results/comparisons_memory_splatter_R.csv")), 
						row.names = FALSE, 
						append = TRUE,
						quote = FALSE
						)

	}

}