package gff_test

import (
	"bytes"
	"fmt"
	"os"
	"path/filepath"
	"testing"

	"github.com/TimothyStiles/poly/bio/gff"
	"github.com/TimothyStiles/poly/transform"
)

// This example shows how to open a gff file and search for a gene given its
// locus tag. We then display the EC number of that particular gene.
func Example_basic() {
	sequence, _ := gff.Read("../../data/ecoli-mg1655-short.gff")
	for _, feature := range sequence.Features {
		if feature.Attributes["locus_tag"] == "b0003" {
			fmt.Println(feature.Attributes["EC_number"])
		}
	}
	// Output: 2.7.1.39
}

func ExampleRead() {
	sequence, _ := gff.Read("../../data/ecoli-mg1655-short.gff")
	fmt.Println(sequence.Meta.Name)
	// Output: U00096.3
}

func ExampleParse() {
	file, _ := os.Open("../../data/ecoli-mg1655-short.gff")
	sequence, _ := gff.Parse(file)

	fmt.Println(sequence.Meta.Name)
	// Output: U00096.3
}

func ExampleBuild() {
	sequence, _ := gff.Read("../../data/ecoli-mg1655-short.gff")
	gffBytes, _ := gff.Build(sequence)
	gffReader := bytes.NewReader(gffBytes)
	reparsedSequence, _ := gff.Parse(gffReader)

	fmt.Println(reparsedSequence.Meta.Name)
	// Output: U00096.3
}

func ExampleWrite() {
	tmpDataDir, err := os.MkdirTemp("", "data-*")
	if err != nil {
		fmt.Println(err.Error())
	}
	defer os.RemoveAll(tmpDataDir)

	sequence, _ := gff.Read("../../data/ecoli-mg1655-short.gff")

	tmpGffFilePath := filepath.Join(tmpDataDir, "ecoli-mg1655-short.gff")
	_ = gff.Write(sequence, tmpGffFilePath)

	testSequence, _ := gff.Read(tmpGffFilePath)

	fmt.Println(testSequence.Meta.Name)
	// Output: U00096.3
}

func ExampleGff_AddFeature() {
	// Sequence for greenflourescent protein (GFP) that we're using as test data for this example.
	gfpSequence := "ATGGCTAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCTACATACGGAAAGCTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCATATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAACGCACTATATCTTTCAAAGATGACGGGAACTACAAGACGCGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATCGTATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTCGGACACAAACTCGAGTACAACTATAACTCACACAATGTATACATCACGGCAGACAAACAAAAGAATGGAATCAAAGCTAACTTCAAAATTCGCCACAACATTGAAGATGGATCCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCGACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGCGTGACCACATGGTCCTTCTTGAGTTTGTAACTGCTGCTGGGATTACACATGGCATGGATGAGCTCTACAAATAA"

	// initialize sequence and feature structs.
	var sequence gff.Gff
	var feature gff.Feature

	// set the initialized sequence struct's sequence.
	sequence.Sequence = gfpSequence

	// Set the initialized feature name and sequence location.
	feature.Location = gff.Location{}
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
	var sequence gff.Gff
	var feature gff.Feature

	// set the initialized sequence struct's sequence.
	sequence.Sequence = gfpSequence

	// Set the initialized feature name and sequence location.
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

func TestFeature_GetSequence(t *testing.T) {
	// This test is a little too complex and contrived for an example function.
	// Essentially, it's testing GetSequence()'s ability to parse and retrieve sequences from complex location structures.
	// This was originally covered in the old package system  it was not covered in the new package system so I decided to include it here.

	// Sequence for greenflourescent protein (GFP) that we're using as test data for this example.
	gfpSequence := "ATGGCTAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCTACATACGGAAAGCTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCATATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAACGCACTATATCTTTCAAAGATGACGGGAACTACAAGACGCGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATCGTATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTCGGACACAAACTCGAGTACAACTATAACTCACACAATGTATACATCACGGCAGACAAACAAAAGAATGGAATCAAAGCTAACTTCAAAATTCGCCACAACATTGAAGATGGATCCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCGACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGCGTGACCACATGGTCCTTCTTGAGTTTGTAACTGCTGCTGGGATTACACATGGCATGGATGAGCTCTACAAATAA"

	sequenceLength := len(gfpSequence)

	// Splitting the sequence into two parts to make a multi-location feature.
	sequenceFirstHalf := gfpSequence[:sequenceLength/2]
	sequenceSecondHalf := transform.ReverseComplement(gfpSequence[sequenceLength/2:]) // This feature is reverse complemented.

	// rejoining the two halves into a single string where the second half of the sequence is reverse complemented.
	gfpSequenceModified := sequenceFirstHalf + sequenceSecondHalf

	// initialize sequence and feature structs.
	var sequence gff.Gff
	var feature gff.Feature

	// set the initialized sequence struct's sequence.
	sequence.Sequence = gfpSequenceModified
	// initialize sublocations to be usedin the feature.

	var subLocation gff.Location
	var subLocationReverseComplemented gff.Location

	subLocation.Start = 0
	subLocation.End = sequenceLength / 2

	subLocationReverseComplemented.Start = sequenceLength / 2
	subLocationReverseComplemented.End = sequenceLength
	subLocationReverseComplemented.Complement = true // According to genbank complement means reverse complement. What a country.

	feature.Location.SubLocations = []gff.Location{subLocation, subLocationReverseComplemented}

	// Add the GFP feature to the sequence struct.
	_ = sequence.AddFeature(&feature)

	// get the GFP feature sequence string from the sequence struct.
	featureSequence, _ := feature.GetSequence()

	// check to see if the feature was inserted properly into the sequence.
	if gfpSequence != featureSequence {
		t.Error("Feature sequence was not properly retrieved.")
	}
}
