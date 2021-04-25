package fasta

import (
	"bytes"
	"fmt"
	"os"
)

// ExampleReadFASTA shows basic usage for ReadFASTA.
func ExampleReadFASTA() {
	fastas := ReadFASTA("data/base.fasta")
	fmt.Println(fastas[0].Name)
	// Output: gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]
}

// ExampleParseFASTA shows basic usage for ParseFASTA.
func ExampleParseFASTA() {
	file, _ := os.Open("data/base.fasta")
	fastas := ParseFASTA(file)

	fmt.Println(fastas[0].Name)
	// Output: gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]
}

// ExampleBuildFASTA shows basic usage for BuildFASTA
func ExampleBuildFASTA() {
	fastas := ReadFASTA("data/base.fasta") // get example data
	fasta := BuildFASTA(fastas)            // build a fasta byte array
	firstLine := string(bytes.Split(fasta, []byte("\n"))[0])

	fmt.Println(firstLine)
	// Output: >gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]
}

// ExampleWriteFASTA shows basic usage of the FASTA writer.
func ExampleWriteFASTA() {
	fastas := ReadFASTA("data/base.fasta")       // get example data
	WriteFASTA(fastas, "data/test.fasta")        // write it out again
	testSequence := ReadFASTA("data/test.fasta") // read it in again

	os.Remove("data/test.fasta") // getting rid of test file

	fmt.Println(testSequence[0].Name)
	// Output: gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]
}

// ExampleReadFASTAGz shows basic usage for ReadFASTAGz on a gzip'd file.
func ExampleReadFASTAGz() {
	fastas := ReadFASTAGz("data/uniprot_1mb_test.fasta.gz")
	var name string
	for _, fasta := range fastas {
		name = fasta.Name
	}

	fmt.Println(name)
	// Output: sp|P86857|AGP_MYTCA Alanine and glycine-rich protein (Fragment) OS=Mytilus californianus OX=6549 PE=1 SV=1
}

// ExampleReadFASTAConcurrent shows how to use the concurrent FASTA parser for larger files.
func ExampleReadFASTAConcurrent() {
	fastas := make(chan Fasta, 1000)
	go ReadFASTAGzConcurrent("data/uniprot_1mb_test.fasta.gz", fastas)
	var name string
	for fasta := range fastas {
		name = fasta.Name
	}

	fmt.Println(name)
	// Output: sp|P86857|AGP_MYTCA Alanine and glycine-rich protein (Fragment) OS=Mytilus californianus OX=6549 PE=1 SV=1
}
