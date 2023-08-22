package fasta2_test

import (
	"fmt"

	"github.com/TimothyStiles/poly/io/fasta2"
)

// ExampleReadFile shows basic usage for ReadFile
func ExampleReadFile() {
	fastas, _ := fasta2.ReadFile("testdata/base.fasta")

	fmt.Println(fastas[0].Header)
	// Output: gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]
}
