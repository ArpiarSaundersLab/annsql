# Interoperability

Behind the scenes, **AnnSQL** uses the **DuckDB** in-process analytical database engine to create the `asql` on-disk database. The database may be accessed and queried using any of the languages below, without any need for conversion. Simply open the database as shown in the example below, then craft any SQL query to access your data with blazing fast speed. 

## C
```c
#include <stdio.h>
#include "duckdb.h"

int main() {
    duckdb_database db;
    duckdb_connection conn;
    duckdb_result result;

    if (duckdb_open("db/pbmc.asql", &db) == DuckDBError) {
        printf("Failed to open database\n");
        return 1;
    }
    duckdb_connect(db, &conn);

    if (duckdb_query(conn, "SELECT * FROM X LIMIT 5", &result) == DuckDBError) {
        printf("Query failed\n");
    }
    
    duckdb_destroy_result(&result);
    duckdb_disconnect(&conn);
    duckdb_close(&db);
    return 0;
}
```

## Command Line
```sh
duckdb db/pbmc.asql -c "SELECT * FROM X LIMIT 5"
```

## Go
```go
package main

import (
	"fmt"
	"github.com/marcboeker/go-duckdb"
)

func main() {
	db, _ := duckdb.Open("db/pbmc.asql")
	defer db.Close()

	rows, _ := db.Query("SELECT * FROM X LIMIT 5")
	defer rows.Close()

	for rows.Next() {
		var data string
		rows.Scan(&data)
		fmt.Println(data)
	}
}
```

## Java
```java
import java.sql.*;
import org.duckdb.DuckDBConnection;

public class DuckDBExample {
    public static void main(String[] args) throws Exception {
        Connection conn = (DuckDBConnection) DriverManager.getConnection("jdbc:duckdb:db/pbmc.asql");
        Statement stmt = conn.createStatement();
        ResultSet rs = stmt.executeQuery("SELECT * FROM X LIMIT 5");

        while (rs.next()) {
            System.out.println(rs.getString(1));
        }

        rs.close();
        stmt.close();
        conn.close();
    }
}
```

## Julia
```julia
using DuckDB
conn = DBInterface.connect(DuckDB.DB, "db/pbmc.asql")
result = DBInterface.execute(conn, "SELECT * FROM X LIMIT 5")
for row in result
    println(row)
end
DBInterface.close!(conn)
```

## Node
```javascript
const duckdb = require('duckdb');

const db = new duckdb.Database('db/pbmc.asql');
db.all("SELECT * FROM X LIMIT 5", (err, rows) => {
    if (err) throw err;
    console.log(rows);
    db.close();
});
```

## Python
```python
from AnnSQL import AnnSQL

asql = AnnSQL(db="db/pbmc3k.asql")

result = asql.query("SELECT * FROM X LIMIT 5")

print(result)
conn.close()
```

## R
```r
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

## Rust
```rust
use duckdb::{Connection, Result};

fn main() -> Result<()> {
    let conn = Connection::open("db/pbmc.asql")?;
    let mut stmt = conn.prepare("SELECT * FROM X LIMIT 5")?;
    let rows = stmt.query_map([], |row| row.get::<_, String>(0))?;
    for row in rows {
        println!("{}", row?);
    }
    Ok(())
}
```
