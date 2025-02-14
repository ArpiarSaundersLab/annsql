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
files <- c(1,2,3,4,5,6)

for (file in files){

	#iterate through the array
	for (i in n_cells) {

		print(paste("File", file))
		print(paste("Starting", i))

		df <- data.frame(matrix(ncol = 8))
		colnames(df) <- c("size", "method", "memory", "filter", "memory_log", "memory_log10", "size_log", "size_log10")

		#open the seurat_obj file
		total_loading_mb <- profile_memory_usage({
			seurat_obj <- readRDS(paste0("/home/kenny/Documents/OHSU/Projects/AnnSql/examples/rdata/random/data_", i, ".rds"))
		})

   		print(paste("Total loading memory:", total_loading_mb, "MB"))

		# memory of the filter 1
		memory_filter_1 <- profile_memory_usage({
			filtered_seurat = (subset(seurat_obj, cells = colnames(subset(seurat_obj, subset = `gene-1` > 0.5, return.null=TRUE))))
		})

		total_filter_1_mb <- total_loading_mb + memory_filter_1
		print(paste("Memory of filter 1 Only:", memory_filter_1, "MB"))
		print(paste("Memory of filter 1 Total:", total_filter_1_mb, "MB"))

		#insert into the dataframe
   		df <- rbind(df, c(i, "Seurat", total_filter_1_mb, "Filter1", log(total_filter_1_mb), log10(total_filter_1_mb), log(as.numeric(i)), log10(as.numeric(i))))


		# memory of the filter 2
		memory_filter_2 <- profile_memory_usage({
			filtered_seurat = (subset(seurat_obj, cells = colnames(subset(seurat_obj, subset = `gene-1` > 0.5 & `gene-2` < 0.5, return.null=TRUE))))
		})

		total_filter_2_mb <- total_loading_mb + memory_filter_2
		print(paste("Memory of filter 2 Only:", memory_filter_2, "MB"))
		print(paste("Memory of filter 2 Total:", total_filter_2_mb, "MB"))

		#insert into the dataframe
   		df <- rbind(df, c(i, "Seurat", total_filter_2_mb, "Filter2", log(total_filter_2_mb), log10(total_filter_2_mb), log(as.numeric(i)), log10(as.numeric(i))))

		# memory of the filter 3
		memory_filter_3 <- profile_memory_usage({
			filtered_seurat = (subset(seurat_obj, cells = colnames(subset(seurat_obj, subset = `gene-1` > 0.5 & `gene-2` < 0.5 & `gene-3` > 0.5, return.null=TRUE))))
		})

		total_filter_3_mb <- total_loading_mb + memory_filter_3
		print(paste("Memory of filter 3 Only:", memory_filter_3, "MB"))
		print(paste("Memory of filter 3 Total:", total_filter_3_mb, "MB"))

		#insert into the dataframe
   		df <- rbind(df, c(i, "Seurat", total_filter_3_mb, "Filter3", log(total_filter_3_mb), log10(total_filter_3_mb), log(as.numeric(i)), log10(as.numeric(i))))

		# memory of the filter 4
		memory_filter_4 <- profile_memory_usage({
			filtered_seurat = (subset(seurat_obj, cells = colnames(subset(seurat_obj, subset = `gene-1` > 0.5 & `gene-2` < 0.5 & `gene-3` > 0.5 & `gene-4` < 0.5, return.null=TRUE))))
		})

		total_filter_4_mb <- total_loading_mb + memory_filter_4
		print(paste("Memory of filter 4 Only:", memory_filter_4, "MB"))
		print(paste("Memory of filter 4 Total:", total_filter_4_mb, "MB"))

		#insert into the dataframe
   		df <- rbind(df, c(i, "Seurat", total_filter_4_mb, "Filter4", log(total_filter_4_mb), log10(total_filter_4_mb), log(as.numeric(i)), log10(as.numeric(i))))

		# memory of the filter 5
		memory_filter_5 <- profile_memory_usage({
			seurat_obj$gene_3_div_gene_4 <- seurat_obj[["RNA"]]@data["gene-3", ] / seurat_obj[["RNA"]]@data["gene-4", ]
			filtered_seurat <- subset(seurat_obj, subset = (`gene-1` > 0.5 & `gene-2` < 0.3) & (gene_3_div_gene_4 > 0 & gene_3_div_gene_4 < 0.5))
		})

		total_filter_5_mb <- total_loading_mb + memory_filter_5
		print(paste("Memory of filter 5 Only:", memory_filter_5, "MB"))
		print(paste("Memory of filter 5 Total:", total_filter_5_mb, "MB"))

		#insert into the dataframe
   		df <- rbind(df, c(i, "Seurat", total_filter_5_mb, "Filter5", log(total_filter_5_mb), log10(total_filter_5_mb), log(as.numeric(i)), log10(as.numeric(i))))

		# memory of the filter 6
		memory_filter_6 <- profile_memory_usage({
			cell_types <- seurat_obj$cell_type
			avg_gene_1 <- aggregate(seurat_obj[["RNA"]]@data["gene-1", ] ~ cell_types, FUN = mean)
			colnames(avg_gene_1) <- c("cell_type", "avg_gene_1")
		})

		total_filter_6_mb <- total_loading_mb + memory_filter_6
		print(paste("Memory of filter 6 Only:", memory_filter_6, "MB"))
		print(paste("Memory of filter 6 Total:", total_filter_6_mb, "MB"))

		#insert into the dataframe
   		df <- rbind(df, c(i, "Seurat", total_filter_6_mb, "Filter6", log(total_filter_6_mb), log10(total_filter_6_mb), log(as.numeric(i)), log10(as.numeric(i))))

		# print the dataframe
		df <- df[-1, ]

		#append the dataframe to the csv file
		write.table(df, file = paste0("/home/kenny/Documents/OHSU/Projects/AnnSql/examples/results/comparisons_random_memory_",file,".csv"), 
						sep = ",", 
						col.names = !file.exists(paste0("/home/kenny/Documents/OHSU/Projects/AnnSql/examples/results/comparisons_random_memory_",file,".csv")), 
						row.names = FALSE, 
						append = TRUE,
						quote = FALSE
						)

	}

}