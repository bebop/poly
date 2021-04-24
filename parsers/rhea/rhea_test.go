package rhea

import (
	"database/sql"
	"fmt"
	_ "github.com/mattn/go-sqlite3"
	"log"
	"os"
	"testing"
)

var db *sql.DB
var rhea Rhea

func TestMain(m *testing.M) {
	// Read the Compressed Rhea XML to bytes
	rheaBytes, err := ReadRhea("data/rhea.rdf.gz")
	if err != nil {
		log.Fatalf("Failed to read rhea")
	}

	// Parse the Rhea bytes into the rhea struct
	rhea, err = ParseRhea(rheaBytes)
	if err != nil {
		log.Fatalf("Failed to parse rhea")
	}

	// Start running tests
	code := m.Run()
	os.Exit(code)
}

func ExampleExportRhea() {
	// Convert rhea to JSON
	rheaJson, _ := rhea.Export()

	fmt.Println(string(rheaJson)[:100])
	// Output: {"reactionParticipants":[{"reactionside":"http://rdf.rhea-db.org/10000_L","contains":1,"containsn":f
}

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

	err = InsertRheaSqlite(db, rhea)
	if err != nil {
		log.Fatalf("InsertRhea Failed: %s", err)
	}
}
