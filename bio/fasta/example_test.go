package fasta_test

import (
	_ "embed"
	"errors"
	"fmt"
	"io"
	"os"
	"strings"

	"github.com/TimothyStiles/poly/bio/fasta"
)

// This example shows how to open a file with the fasta parser. The sequences
// within that file can then be analyzed further with different software.
func Example_basic() {
	fastaString := ">testing\nATGC\n"
	file := strings.NewReader(fastaString)   // this "file" can be replaced with os.Open
	parser := fasta.NewParser(file, 32*1024) // Initialize a parser for this file
	var records []fasta.Record
	for {
		record, err := parser.Next()
		if err != nil {
			if errors.Is(err, io.EOF) {
				err = nil // EOF not treated as parsing error.
			}
			if err != nil {
				fmt.Printf("Got error: %s\n", err) // check if there was a parse error
			}
			break
		}
		records = append(records, *record)
	}
	fmt.Println(records[0].Sequence)
	// Output: ATGC
}

// ExampleWrite shows basic usage of the  writer.
func ExampleRecord_WriteTo() {
	// Get a fasta record
	fastaString := ">testing\nATGC\n"
	fakeFile := strings.NewReader(fastaString)
	parser := fasta.NewParser(fakeFile, 32*1024)
	record, _ := parser.Next()

	// Write it to a file
	file, _ := os.Create("data/test.fasta")
	_, _ = record.WriteTo(file)
	file.Close()

	// Read that file
	file2, _ := os.Open("data/test.fasta")
	parser = fasta.NewParser(file2, 32*1024)
	newRecord, _ := parser.Next()
	file2.Close()

	// Remove the test file
	os.Remove("data/test.fasta")

	fmt.Println(newRecord.Identifier)
	// Output: testing
}
