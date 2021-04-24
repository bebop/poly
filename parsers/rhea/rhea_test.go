package rhea

import (
	"database/sql"
	_ "github.com/mattn/go-sqlite3"
	"log"
	"testing"
)

var db *sql.DB

func TestInsertRheaSqlite(t *testing.T) {
	// Connect to in-memory sqlite database
	var err error
	db, err = sql.Open("sqlite3", ":memory:")
	if err != nil {
		t.Fatalf("Could not connect to in-memory sqlite database")
	}

	err = CreateRheaSqlite(db)
	if err != nil {
		log.Fatalf("CreateDatabase() Failed: %s", err)
	}

	err = InsertRheaSqlite(db, "data/rhea.rdf.gz")
	if err != nil {
		log.Fatalf("InsertRhea Failed: %s", err)
	}
}
