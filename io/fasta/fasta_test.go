package fasta

import (
	"bytes"
	"fmt"
	"os"
)

// ExampleRead shows basic usage for Read.
func ExampleRead() {
	fastas, _ := Read("data/base.fasta")
	fmt.Println(fastas[0].Name)
	// Output: gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]
}

// ExampleParse shows basic usage for Parse.
func ExampleParse() {
	file, _ := os.Open("data/base.fasta")
	fastas, _ := Parse(file)

	fmt.Println(fastas[0].Name)
	// Output: gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]
}

// ExampleBuild shows basic usage for Build
func ExampleBuild() {
	fastas, _ := Read("data/base.fasta") // get example data
	fasta, _ := Build(fastas)            // build a fasta byte array
	firstLine := string(bytes.Split(fasta, []byte("\n"))[0])

	fmt.Println(firstLine)
	// Output: >gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]
}

// ExampleWrite shows basic usage of the  writer.
func ExampleWrite() {
	fastas, _ := Read("data/base.fasta")       // get example data
	_ = Write(fastas, "data/test.fasta")       // write it out again
	testSequence, _ := Read("data/test.fasta") // read it in again

	os.Remove("data/test.fasta") // getting rid of test file

	fmt.Println(testSequence[0].Name)
	// Output: gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]
}

// ExampleReadGz shows basic usage for ReadGz on a gzip'd file.
func ExampleReadGz() {
	fastas, _ := ReadGz("data/uniprot_1mb_test.fasta.gz")
	var name string
	for _, fasta := range fastas {
		name = fasta.Name
	}

	fmt.Println(name)
	// Output: sp|P86857|AGP_MYTCA Alanine and glycine-rich protein (Fragment) OS=Mytilus californianus OX=6549 PE=1 SV=1
}

// ExampleReadGzConcurrent shows how to use the concurrent  parser for larger files.
func ExampleReadGzConcurrent() {
	fastas := make(chan Fasta, 1000)
	go ReadGzConcurrent("data/uniprot_1mb_test.fasta.gz", fastas)
	var name string
	for fasta := range fastas {
		name = fasta.Name
	}

	fmt.Println(name)
	// Output: sp|P86857|AGP_MYTCA Alanine and glycine-rich protein (Fragment) OS=Mytilus californianus OX=6549 PE=1 SV=1
}

// ExampleReadConcurrent shows how to use the concurrent  parser for decompressed fasta files.
func ExampleReadConcurrent() {
	fastas := make(chan Fasta, 100)
	go ReadConcurrent("data/base.fasta", fastas)
	var name string
	for fasta := range fastas {
		name = fasta.Name
	}

	fmt.Println(name)
	// Output: MCHU - Calmodulin - Human, rabbit, bovine, rat, and chicken
}
