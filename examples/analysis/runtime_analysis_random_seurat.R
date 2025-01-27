library(Seurat)

#the different file increments in an array (must be in quotes R, converts to 1e5 at 100k if not in quotes)
n_cells <- c("1000", "5000", "10000", "15000", "20000", "25000", "30000", "35000", "40000", "45000", "50000", "75000", "100000", "250000")

num_iterations <- 30
elapsed_times <- numeric(num_iterations)
random_min <- 1
random_max <- 999

#filter 1
for (i in seq_len(num_iterations)) {
  set.seed(as.numeric(Sys.time()+i))
  random_gene_number <- sample(random_min:random_max, 1)
  random_gene <- paste0("gene-", random_gene_number)
  elapsed_times[i] <- system.time({
    filtered_seurat <- subset(
      seurat_obj,
      subset = !!rlang::sym(random_gene) > 0.5, return.null=TRUE
    )
  })["elapsed"]
  print(paste("Iteration", i, "gene",random_gene, "Elapsed time:", elapsed_times[i]))
}
par(mar = c(2, 2, 2, 2))  # c(bottom, left, top, right)
boxplot(elapsed_times, main = "Elapsed Times", ylab = "Time (seconds)")
points(jitter(rep(1, length(elapsed_times))), elapsed_times, pch = 19)

#filter 2
system.time({
  filtered_seurat = (subset(seurat_obj, cells = colnames(subset(seurat_obj, subset = `gene-1` > 0.5 & `gene-2` < 0.5, return.null=TRUE))))
})["elapsed"]

#filter 3
system.time({
  filtered_seurat = (subset(seurat_obj, cells = colnames(subset(seurat_obj, subset = `gene-1` > 0.5 & `gene-2` < 0.5 & `gene-3` > 0.5, return.null=TRUE))))
})["elapsed"]

#filter 4
system.time({
  filtered_seurat = (subset(seurat_obj, cells = colnames(subset(seurat_obj, subset = `gene-1` > 0.5 & `gene-2` < 0.5 & `gene-3` > 0.5 & `gene-4` < 0.5, return.null=TRUE))))
})["elapsed"]

#filter 5
system.time({
  seurat_obj$gene_3_div_gene_4 <- seurat_obj[["RNA"]]@data["gene-3", ] / seurat_obj[["RNA"]]@data["gene-4", ]
  filtered_seurat <- subset(seurat_obj, subset = (`gene-1` > 0.5 & `gene-2` < 0.3) & (gene_3_div_gene_4 > 0 & gene_3_div_gene_4 < 0.5))
})["elapsed"]

#filter 6
system.time({
  cell_types <- seurat_obj$cell_type
  avg_gene_1 <- aggregate(seurat_obj[["RNA"]]@data["gene-1", ] ~ cell_types, FUN = mean)
  colnames(avg_gene_1) <- c("cell_type", "avg_gene_1")
})["elapsed"]



