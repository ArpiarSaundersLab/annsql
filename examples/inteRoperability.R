#install the duckdb api dependency and load the library
if (!requireNamespace("duckdb", quietly = TRUE)) {
	install.packages("duckdb")
}

library("duckdb")

#connect to the AnnSQL db using the duckdb driver
con <- dbConnect(duckdb::duckdb(), dbdir = "db/pbmc3k.asql")

#run any query your heart desires
results <- dbGetQuery(con, "SELECT * FROM adata LIMIT 5")

#take a look at the results
print(results)

#don't forget to close the door
dbDisconnect(con)




