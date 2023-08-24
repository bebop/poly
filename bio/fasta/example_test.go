package fasta_test

import (
	"bytes"
	_ "embed"
	"fmt"
	"os"
	"strings"

	"github.com/TimothyStiles/poly/bio/fasta"
)

//go:embed data/base.fasta
var baseFasta string

// This example shows how to open a file with the fasta parser. The sequences
// within that file can then be analyzed further with different software.
func Example_basic() {
	fastas, _ := fasta.Read("data/base.fasta")
	fmt.Println(fastas[1].Sequence)
	// Output: ADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNGTIDFPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLTDEEVDEMIREADIDGDGQVNYEEFVQMMTAK*
}

// ExampleRead shows basic usage for Read.
func ExampleRead() {
	fastas, _ := fasta.Read("data/base.fasta")
	fmt.Println(fastas[0].Name)
	// Output: gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]
}

// ExampleParse shows basic usage for Parse.
func ExampleParse() {
	file, _ := os.Open("data/base.fasta")
	fastas, _ := fasta.Parse(file)

	fmt.Println(fastas[0].Name)
	// Output: gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]
}

// ExampleBuild shows basic usage for Build
func ExampleBuild() {
	fastas, _ := fasta.Read("data/base.fasta") // get example data
	var buffer bytes.Buffer                    // Initialize a buffer to write fastas into
	for _, fasta := range fastas {
		_ = fasta.Write(&buffer) // build a fasta byte array

	}
	firstLine := string(bytes.Split(buffer.Bytes(), []byte("\n"))[0])

	fmt.Println(firstLine)
	// Output: >gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]
}

// ExampleWrite shows basic usage of the  writer.
func ExampleWrite() {
	fastas, _ := fasta.Read("data/base.fasta")       // get example data
	_ = fasta.WriteFile(fastas, "data/test.fasta")   // write it out again
	testSequence, _ := fasta.Read("data/test.fasta") // read it in again

	os.Remove("data/test.fasta") // getting rid of test file

	fmt.Println(testSequence[0].Name)
	// Output: gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]
}

// ExampleReadGz shows basic usage for ReadGz on a gzip'd file.
func ExampleReadGz() {
	fastas, _ := fasta.ReadGz("data/uniprot_1mb_test.fasta.gz")
	var name string
	for _, fasta := range fastas {
		name = fasta.Name
	}

	fmt.Println(name)
	// Output: sp|P86857|AGP_MYTCA Alanine and glycine-rich protein (Fragment) OS=Mytilus californianus OX=6549 PE=1 SV=1
}

// ExampleReadGzConcurrent shows how to use the concurrent  parser for larger files.
func ExampleReadGzConcurrent() {
	fastas := make(chan fasta.Fasta, 1000)
	go fasta.ReadGzConcurrent("data/uniprot_1mb_test.fasta.gz", fastas)
	var name string
	for fasta := range fastas {
		name = fasta.Name
	}

	fmt.Println(name)
	// Output: sp|P86857|AGP_MYTCA Alanine and glycine-rich protein (Fragment) OS=Mytilus californianus OX=6549 PE=1 SV=1
}

// ExampleReadConcurrent shows how to use the concurrent  parser for decompressed fasta files.
func ExampleReadConcurrent() {
	fastas := make(chan fasta.Fasta, 100)
	go fasta.ReadConcurrent("data/base.fasta", fastas)
	var name string
	for fasta := range fastas {
		name = fasta.Name
	}

	fmt.Println(name)
	// Output: MCHU - Calmodulin - Human, rabbit, bovine, rat, and chicken
}

func ExampleParser() {
	parser := fasta.NewParser(strings.NewReader(baseFasta), 256)
	for {
		fasta, _, err := parser.ParseNext()
		if err != nil {
			fmt.Println(err)
			break
		}
		fmt.Println(fasta.Name)
	}
	//Output:
	// gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]
	// MCHU - Calmodulin - Human, rabbit, bovine, rat, and chicken
	// EOF
}
