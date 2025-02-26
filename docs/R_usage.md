# Opening AnnSQL Databases in R
Annsql uses the duckdb api to build the database and perform operations. You can access the .asql file in R by simply installing the duckdb package and loading the library and database. Once you've loaded the db you can run any query your heart desires.
<br>
<br>

```
#install the duckdb api dependency and load the library
if (!requireNamespace("duckdb", quietly = TRUE)) {
	install.packages("duckdb")
}

#import duckdb
library("duckdb")

#connect to the AnnSQL db using the duckdb driver
con <- dbConnect(duckdb::duckdb(), dbdir = "db/pbmc3k.asql")

#run any query your heart desires
results <- dbGetQuery(con, "SELECT * FROM adata LIMIT 5")

#take a look at the results
print(results)

#don't forget to close the door
dbDisconnect(con)
```
<br>

While extended functionality in our Python based package is not currently implemented in R. AnnSQL databases can be opened, queried, and explored in R. 