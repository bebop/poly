package main

import (
	"bytes"
	"io/ioutil"
	"os"
	"strings"
	"testing"

	"github.com/TimothyStiles/poly"
	"github.com/google/go-cmp/cmp"
	"github.com/google/go-cmp/cmp/cmpopts"
)

/******************************************************************************
Oct, 15, 2020

Testing command line utilities via subroutines can be annoying so
if you're doing it from the commandline be sure to compile first.
From the project's root directory often use:

go build && go install && go test -v ./...

To accurately test your commands you MUST make sure to rebuild and reinstall
before you run your tests. Otherwise your system version will be out of
date and will give you results using an older build.

Happy hacking,
Tim


TODO:

write subtest to check for empty output before merge
******************************************************************************/

func TestConvertPipe(t *testing.T) {

	var writeBuffer bytes.Buffer
	app := application()
	app.Writer = &writeBuffer

	args := os.Args[0:1]                                // Name of the program.
	args = append(args, "c", "-i", "gbk", "-o", "json") // Append a flag

	file, _ := ioutil.ReadFile("../data/puc19.gbk")
	app.Reader = bytes.NewReader(file)

	err := app.Run(args)

	if err != nil {
		t.Fatalf("Run error: %s", err)
	}

	// getting test sequence from non-pipe io to compare against io to stdout
	baseTestSequence := poly.ReadGbk("../data/puc19.gbk")

	pipeOutputTestSequence := poly.ParseJSON(writeBuffer.Bytes())

	if diff := cmp.Diff(baseTestSequence, pipeOutputTestSequence, cmpopts.IgnoreFields(poly.Feature{}, "ParentAnnotatedSequence")); diff != "" {
		t.Errorf(" mismatch from convert pipe input test (-want +got):\n%s", diff)
	}

}

func TestConvertFile(t *testing.T) {

	app := application()

	args := os.Args[0:1]                                                                // Name of the program.
	args = append(args, "c", "-o", "json", "../data/puc19.gbk", "../data/t4_intron.gb") // Append a flag

	err := app.Run(args)

	if err != nil {
		t.Fatalf("Run error: %s", err)
	}

	puc19InputTestSequence := poly.ReadGbk("../data/puc19.gbk")
	puc19OutputTestSequence := poly.ReadJSON("../data/puc19.json")

	//clearing test data.
	os.Remove("../data/puc19.json")

	// compared input gff from resulting output json. Fail test and print diff if error.
	if diff := cmp.Diff(puc19InputTestSequence, puc19OutputTestSequence, cmpopts.IgnoreFields(poly.Feature{}, "ParentAnnotatedSequence")); diff != "" {
		t.Errorf(" mismatch from concurrent gbk input test (-want +got):\n%s", diff)
	}

	t4InputTestSequence := poly.ReadGbk("../data/t4_intron.gb")
	t4OutputTestSequence := poly.ReadJSON("../data/t4_intron.json")

	// clearing test data.
	os.Remove("../data/t4_intron.json")

	// compared input gbk from resulting output json. Fail test and print diff if error.
	if diff := cmp.Diff(t4InputTestSequence, t4OutputTestSequence, cmpopts.IgnoreFields(poly.Feature{}, "ParentAnnotatedSequence")); diff != "" {
		t.Errorf(" mismatch from concurrent gbk input test (-want +got):\n%s", diff)
	}
}
func TestHashFile(t *testing.T) {

	puc19GbkBlake3Hash := "4b0616d1b3fc632e42d78521deb38b44fba95cca9fde159e01cd567fa996ceb9"
	var writeBuffer bytes.Buffer

	app := application()
	app.Writer = &writeBuffer

	// testing file matching hash
	args := os.Args[0:1]                             // Name of the program.
	args = append(args, "hash", "../data/puc19.gbk") // Append a flag

	err := app.Run(args)

	if err != nil {
		t.Fatalf("Run error: %s", err)
	}

	if writeBuffer.Len() == 0 {
		t.Error("TestHash did not write output to desired writer.")
	}

	hashOutputString := strings.TrimSpace(writeBuffer.String())
	if hashOutputString != puc19GbkBlake3Hash {
		t.Errorf("TestHashFile has failed. Returned %q, want %q", hashOutputString, puc19GbkBlake3Hash)
	}

}

func TestHashPipe(t *testing.T) {

	puc19GbkBlake3Hash := "4b0616d1b3fc632e42d78521deb38b44fba95cca9fde159e01cd567fa996ceb9"
	var writeBuffer bytes.Buffer

	// create a mock application
	app := application()
	app.Writer = &writeBuffer
	file, _ := ioutil.ReadFile("../data/puc19.gbk")
	app.Reader = bytes.NewReader(file)

	args := os.Args[0:1]                     // Name of the program.
	args = append(args, "hash", "-i", "gbk") // Append a flag

	err := app.Run(args)

	if err != nil {
		t.Fatalf("Run error: %s", err)
	}

	hashOutputString := strings.TrimSpace(writeBuffer.String())

	if hashOutputString != puc19GbkBlake3Hash {
		t.Errorf("TestHashPipe has failed. Returned %q, want %q", hashOutputString, puc19GbkBlake3Hash)
	}

}

func TestHashJSON(t *testing.T) {
	// testing json write output

	puc19GbkBlake3Hash := "4b0616d1b3fc632e42d78521deb38b44fba95cca9fde159e01cd567fa996ceb9"

	app := application()

	args := os.Args[0:1]                                           // Name of the program.
	args = append(args, "hash", "-o", "json", "../data/puc19.gbk") // Append a flag
	err := app.Run(args)

	if err != nil {
		t.Fatalf("Run error: %s", err)
	}

	hashOutputString := poly.ReadJSON("../data/puc19.json").SequenceHash
	os.Remove("../data/puc19.json")

	if hashOutputString != puc19GbkBlake3Hash {
		t.Errorf("TestHashJSON has failed. Returned %q, want %q", hashOutputString, puc19GbkBlake3Hash)
	}

}

func TestOptimizeString(t *testing.T) {

	gfpTranslation := "MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK*"
	var writeBuffer bytes.Buffer

	app := application()
	app.Writer = &writeBuffer
	app.Reader = bytes.NewBufferString(gfpTranslation)

	args := os.Args[0:1]                                      // Name of the program.
	args = append(args, "op", "-aa", "-wt", "data/puc19.gbk") // Append a flag
	err := app.Run(args)

	if err != nil {
		t.Fatalf("Run error: %s", err)
	}
	app.Reader = os.Stdin

	// should return codon optimized sequence
	optimizeOutputString := strings.TrimSpace(writeBuffer.String())

	translation := poly.Translate(optimizeOutputString, poly.DefaultCodonTablesByNumber[1])

	if translation != gfpTranslation {
		t.Errorf("TestOptimizeCommand for string output has failed. Returned %q, want %q", translation, gfpTranslation)
	}

}

func TestTranslationString(t *testing.T) {
	gfpDnaSequence := "ATGGCTAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCTACATACGGAAAGCTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCATATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAACGCACTATATCTTTCAAAGATGACGGGAACTACAAGACGCGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATCGTATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTCGGACACAAACTCGAGTACAACTATAACTCACACAATGTATACATCACGGCAGACAAACAAAAGAATGGAATCAAAGCTAACTTCAAAATTCGCCACAACATTGAAGATGGATCCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCGACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGCGTGACCACATGGTCCTTCTTGAGTTTGTAACTGCTGCTGGGATTACACATGGCATGGATGAGCTCTACAAATAA"
	gfpTranslation := "MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK*"

	var writeBuffer bytes.Buffer

	app := application()
	app.Writer = &writeBuffer
	app.Reader = bytes.NewBufferString(gfpDnaSequence)

	args := os.Args[0:1]                   // Name of the program.
	args = append(args, "tr", "-ct", "11") // Append a flag
	err := app.Run(args)

	if err != nil {
		t.Fatalf("Run error: %s", err)
	}
	app.Reader = os.Stdin

	translation := strings.TrimSpace(writeBuffer.String())

	if translation != gfpTranslation {
		t.Errorf("TestTranslationString for string output has failed. Returned %q, want %q", translation, gfpTranslation)
	}

}
