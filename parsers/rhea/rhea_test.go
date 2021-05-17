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
	rheaBytes, err := ReadGzippedXml("data/rhea.rdf.gz")
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

	err = RheaSqlite(db, rhea)
	if err != nil {
		log.Fatalf("InsertRhea Failed: %s", err)
	}
}

func TestReadRhea2Uniprot(t *testing.T) {
	lines := make(chan Rhea2Uniprot, 100)
	go ReadRhea2UniprotTrembl("data/rhea2uniprot_sprot_minimized.tsv.gz", lines)

	var line Rhea2Uniprot
	for l := range lines {
		line = l
	}

	if line.UniprotID != "P06106" {
		log.Fatalf("Got wrong uniprotId. Expected P06106, got %s", line.UniprotID)
	}
}

func ExampleReadRhea2UniprotSprot() {
	lines := make(chan Rhea2Uniprot, 100)
	go ReadRhea2UniprotSprot("data/rhea2uniprot_sprot_minimized.tsv", lines)

	var line Rhea2Uniprot
	for l := range lines {
		line = l
	}

	fmt.Println(line)
	// Output: {10048 UN 10048 P06106}
}
