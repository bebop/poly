package genbank_test

import (
	"bytes"
	"fmt"
	"os"
	"path/filepath"

	"github.com/TimothyStiles/poly/io/genbank"
)

// This example shows how to open a genbank file and search for a gene given
// its name. After finding it, notes about the particular gene are read.
func Example_basic() {
	sequences, _ := genbank.Read("../../data/puc19.gbk")
	for _, feature := range sequences.Features {
		if feature.Attributes["gene"] == "bla" {
			fmt.Println(feature.Attributes["note"])
		}
	}
	// Output: confers resistance to ampicillin, carbenicillin, andrelated antibiotics
}

func ExampleRead() {
	sequence, _ := genbank.Read("../../data/puc19.gbk")
	fmt.Println(sequence.Meta.Locus.ModificationDate)
	// Output: 22-OCT-2019
}

func ExampleWrite() {
	tmpDataDir, err := os.MkdirTemp("", "data-*")
	if err != nil {
		fmt.Println(err.Error())
	}
	defer os.RemoveAll(tmpDataDir)

	sequences, _ := genbank.Read("../../data/puc19.gbk")

	tmpGbkFilePath := filepath.Join(tmpDataDir, "puc19.gbk")
	_ = genbank.Write(sequences, tmpGbkFilePath)

	testSequence, _ := genbank.Read(tmpGbkFilePath)

	fmt.Println(testSequence.Meta.Locus.ModificationDate)
	// Output: 22-OCT-2019
}

func ExampleBuild() {
	sequences, _ := genbank.Read("../../data/puc19.gbk")
	gbkBytes, _ := genbank.Build(sequences)
	testSequence, _ := genbank.Parse(bytes.NewReader(gbkBytes))

	fmt.Println(testSequence.Meta.Locus.ModificationDate)
	// Output: 22-OCT-2019
}

func ExampleParse() {
	file, _ := os.Open("../../data/puc19.gbk")
	sequence, _ := genbank.Parse(file)

	fmt.Println(sequence.Meta.Locus.ModificationDate)
	// Output: 22-OCT-2019
}

func ExampleReadMulti() {
	sequences, err := genbank.ReadMulti("../../data/multiGbk_test.seq")
	if err != nil {
		fmt.Println(err.Error())
	}

	fmt.Println(sequences[1].Meta.Locus.ModificationDate)
	// Output: 05-FEB-1999
}

func ExampleWriteMulti() {
	tmpDataDir, err := os.MkdirTemp("", "data-*")
	if err != nil {
		fmt.Println(err.Error())
	}

	sequences, _ := genbank.ReadMulti("../../data/multiGbk_test.seq")
	tmpGbkFilePath := filepath.Join(tmpDataDir, "multiGbk_test.seq")

	err = genbank.WriteMulti(sequences, tmpGbkFilePath)

	if err != nil {
		fmt.Println(err.Error())
	}

	testSequences, _ := genbank.ReadMulti(tmpGbkFilePath)
	isEqual := sequences[1].Meta.Locus.ModificationDate == testSequences[1].Meta.Locus.ModificationDate
	fmt.Println(isEqual)
	// Output: true
}

func ExampleBuildMulti() {
	sequences, _ := genbank.ReadMulti("../../data/multiGbk_test.seq")
	gbkBytes, _ := genbank.BuildMulti(sequences)
	testSequences, _ := genbank.ParseMulti(bytes.NewReader(gbkBytes))

	isEqual := sequences[1].Meta.Locus.ModificationDate == testSequences[1].Meta.Locus.ModificationDate
	fmt.Println(isEqual)
	// Output: true
}

func ExampleParseMulti() {
	file, _ := os.Open("../../data/multiGbk_test.seq")
	sequences, _ := genbank.ParseMulti(file)
	fmt.Println(sequences[1].Meta.Locus.ModificationDate)
	// Output: 05-FEB-1999
}

func ExampleGenbank_AddFeature() {
	// Sequence for greenflourescent protein (GFP) that we're using as test data for this example.
	gfpSequence := "ATGGCTAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCTACATACGGAAAGCTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCATATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAACGCACTATATCTTTCAAAGATGACGGGAACTACAAGACGCGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATCGTATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTCGGACACAAACTCGAGTACAACTATAACTCACACAATGTATACATCACGGCAGACAAACAAAAGAATGGAATCAAAGCTAACTTCAAAATTCGCCACAACATTGAAGATGGATCCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCGACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGCGTGACCACATGGTCCTTCTTGAGTTTGTAACTGCTGCTGGGATTACACATGGCATGGATGAGCTCTACAAATAA"

	// initialize sequence and feature structs.
	var sequence genbank.Genbank
	var feature genbank.Feature

	// set the initialized sequence struct's sequence.
	sequence.Sequence = gfpSequence

	// Set the initialized feature name and sequence location.
	feature.Description = "Green Fluorescent Protein"
	feature.Location = genbank.Location{}
	feature.Location.Start = 0
	feature.Location.End = len(sequence.Sequence)

	// Add the GFP feature to the sequence struct.
	_ = sequence.AddFeature(&feature)

	// get the GFP feature sequence string from the sequence struct.
	featureSequence, _ := feature.GetSequence()

	// check to see if the feature was inserted properly into the sequence.
	fmt.Println(gfpSequence == featureSequence)

	// Output: true
}

func ExampleFeature_GetSequence() {
	// Sequence for greenflourescent protein (GFP) that we're using as test data for this example.
	gfpSequence := "ATGGCTAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCTACATACGGAAAGCTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCATATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAACGCACTATATCTTTCAAAGATGACGGGAACTACAAGACGCGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATCGTATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTCGGACACAAACTCGAGTACAACTATAACTCACACAATGTATACATCACGGCAGACAAACAAAAGAATGGAATCAAAGCTAACTTCAAAATTCGCCACAACATTGAAGATGGATCCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCGACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGCGTGACCACATGGTCCTTCTTGAGTTTGTAACTGCTGCTGGGATTACACATGGCATGGATGAGCTCTACAAATAA"

	// initialize sequence and feature structs.
	var sequence genbank.Genbank
	var feature genbank.Feature

	// set the initialized sequence struct's sequence.
	sequence.Sequence = gfpSequence

	// Set the initialized feature name and sequence location.
	feature.Description = "Green Fluorescent Protein"
	feature.Location.Start = 0
	feature.Location.End = len(sequence.Sequence)

	// Add the GFP feature to the sequence struct.
	_ = sequence.AddFeature(&feature)

	// get the GFP feature sequence string from the sequence struct.
	featureSequence, _ := feature.GetSequence()

	// check to see if the feature was inserted properly into the sequence.
	fmt.Println(gfpSequence == featureSequence)

	// Output: true
}
