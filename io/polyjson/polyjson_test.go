package polyjson

import (
	"errors"
	"fmt"
	"io/ioutil"
	"os"
	"path/filepath"
	"testing"
	"time"

	"github.com/TimothyStiles/poly/seqhash"
	"github.com/TimothyStiles/poly/transform"
	"github.com/stretchr/testify/assert"
)

/******************************************************************************

JSON related tests begin here.

******************************************************************************/

func Example() {

	// this example also is run by the poly's suite so this just sets up a temporary directory for writing files
	tmpDataDir, err := ioutil.TempDir("", "data-*")
	if err != nil {
		fmt.Println(err.Error())
	}

	defer os.RemoveAll(tmpDataDir)

	// initiate a new polyjson sequence struct
	var sequence Poly

	// define the meta section of our sequence.
	sequence.Meta.Name = "Cat DNA"
	sequence.Meta.Description = "Synthetic Cat DNA for testing purposes."
	sequence.Meta.CreatedBy = "Catz (all you basepair are belong to us)"
	sequence.Meta.CreatedOn = time.Now()
	sequence.Meta.URL = "www.allyourbasepair.com/catz"
	sequence.Meta.CreatedWith = "Poly - The world's most modern, open source software library for engineering organisms."

	// add our sequence string and its hash to use as a unique identifier.
	sequence.Sequence = "CATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCAT"
	sequence.Meta.Hash, _ = seqhash.Hash(sequence.Sequence, "DNA", false, true)

	// add our sequence features
	catFeature := Feature{}
	catFeature.Name = "Cat coding region."
	catFeature.Description = "a cat coding region at the beginning of our sequence."
	catFeature.Type = "CDS"
	catFeature.Location.Start = 0
	catFeature.Location.End = 8
	catFeature.Tags = map[string]string{"product": "cat protein"}

	_ = sequence.AddFeature(&catFeature) // add the feature annotation to our sequence

	// write our sequence to a JSON file
	tmpJSONFilePath := filepath.Join(tmpDataDir, "sample.json")
	_ = Write(sequence, tmpJSONFilePath)

	exportedSequence, _ := Read(tmpJSONFilePath)

	// print our struct DNA sequence
	fmt.Println(exportedSequence.Sequence)
	// Output: CATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCAT
}

func ExampleRead() {
	sequence, _ := Read("../../data/cat.json")

	fmt.Println(sequence.Sequence)
	// Output: CATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCAT
}

func ExampleParse() {
	file, _ := ioutil.ReadFile("../../data/cat.json")
	sequence, _ := Parse(file)

	fmt.Println(sequence.Sequence)
	// Output: CATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCAT
}

func ExampleWrite() {
	tmpDataDir, err := ioutil.TempDir("", "data-*")
	if err != nil {
		fmt.Println(err.Error())
	}
	defer os.RemoveAll(tmpDataDir)

	sequence, _ := Read("../../data/cat.json")

	tmpJSONFilePath := filepath.Join(tmpDataDir, "sample.json")
	_ = Write(sequence, tmpJSONFilePath)

	testSequence, _ := Read(tmpJSONFilePath)

	fmt.Println(testSequence.Sequence)
	// Output: CATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCAT
}

/******************************************************************************

JSON related tests end here.

******************************************************************************/

func ExamplePoly_AddFeature() {

	// Sequence for greenflourescent protein (GFP) that we're using as test data for this example.
	gfpSequence := "ATGGCTAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCTACATACGGAAAGCTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCATATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAACGCACTATATCTTTCAAAGATGACGGGAACTACAAGACGCGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATCGTATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTCGGACACAAACTCGAGTACAACTATAACTCACACAATGTATACATCACGGCAGACAAACAAAAGAATGGAATCAAAGCTAACTTCAAAATTCGCCACAACATTGAAGATGGATCCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCGACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGCGTGACCACATGGTCCTTCTTGAGTTTGTAACTGCTGCTGGGATTACACATGGCATGGATGAGCTCTACAAATAA"

	// initialize sequence and feature structs.
	var sequence Poly
	var feature Feature

	// set the initialized sequence struct's sequence.
	sequence.Sequence = gfpSequence

	// Set the initialized feature name and sequence location.
	feature.Description = "Green Flourescent Protein"
	feature.Location = Location{}
	feature.Location.Start = 0
	feature.Location.End = len(sequence.Sequence)

	// Add the GFP feature to the sequence struct.
	_ = sequence.AddFeature(&feature)

	// get the GFP feature sequence string from the sequence struct.
	featureSequence := feature.GetSequence()

	// check to see if the feature was inserted properly into the sequence.
	fmt.Println(gfpSequence == featureSequence)

	// Output: true
}

func ExampleFeature_GetSequence() {

	// Sequence for greenflourescent protein (GFP) that we're using as test data for this example.
	gfpSequence := "ATGGCTAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCTACATACGGAAAGCTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCATATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAACGCACTATATCTTTCAAAGATGACGGGAACTACAAGACGCGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATCGTATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTCGGACACAAACTCGAGTACAACTATAACTCACACAATGTATACATCACGGCAGACAAACAAAAGAATGGAATCAAAGCTAACTTCAAAATTCGCCACAACATTGAAGATGGATCCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCGACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGCGTGACCACATGGTCCTTCTTGAGTTTGTAACTGCTGCTGGGATTACACATGGCATGGATGAGCTCTACAAATAA"

	// initialize sequence and feature structs.
	var sequence Poly
	var feature Feature

	// set the initialized sequence struct's sequence.
	sequence.Sequence = gfpSequence

	// Set the initialized feature name and sequence location.
	feature.Description = "Green Flourescent Protein"
	feature.Location.Start = 0
	feature.Location.End = len(sequence.Sequence)

	// Add the GFP feature to the sequence struct.
	_ = sequence.AddFeature(&feature)

	// get the GFP feature sequence string from the sequence struct.
	featureSequence := feature.GetSequence()

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
	var sequence Poly
	var feature Feature

	// set the initialized sequence struct's sequence.
	sequence.Sequence = gfpSequenceModified
	// initialize sublocations to be usedin the feature.

	var subLocation Location
	var subLocationReverseComplemented Location

	subLocation.Start = 0
	subLocation.End = sequenceLength / 2

	subLocationReverseComplemented.Start = sequenceLength / 2
	subLocationReverseComplemented.End = sequenceLength
	subLocationReverseComplemented.Complement = true // According to genbank complement means reverse complement. What a country.

	feature.Description = "Green Flourescent Protein"
	feature.Location.SubLocations = []Location{subLocation, subLocationReverseComplemented}

	// Add the GFP feature to the sequence struct.
	_ = sequence.AddFeature(&feature)

	// get the GFP feature sequence string from the sequence struct.
	featureSequence := feature.GetSequence()

	// check to see if the feature was inserted properly into the sequence.
	if gfpSequence != featureSequence {
		t.Error("Feature sequence was not properly retrieved.")
	}

}

func TestParse_error(t *testing.T) {
	unmarshalErr := errors.New("unmarshal error")
	oldUnmarshalFn := unmarshalFn
	unmarshalFn = func(data []byte, v interface{}) error {
		return unmarshalErr
	}
	defer func() {
		unmarshalFn = oldUnmarshalFn
	}()
	_, err := Parse([]byte{})
	assert.EqualError(t, err, unmarshalErr.Error())
}

func TestRead_error(t *testing.T) {
	readErr := errors.New("read error")
	oldReadFileFn := readFileFn
	readFileFn = func(filename string) ([]byte, error) {
		return nil, readErr
	}
	defer func() {
		readFileFn = oldReadFileFn
	}()
	_, err := Read("/tmp/file")
	assert.EqualError(t, err, readErr.Error())
}

func TestWrite_error(t *testing.T) {
	marshalIndentErr := errors.New("marshal indent error")
	oldMarshalIndentFn := marshalIndentFn
	marshalIndentFn = func(v interface{}, prefix, indent string) ([]byte, error) {
		return nil, marshalIndentErr
	}
	defer func() {
		marshalIndentFn = oldMarshalIndentFn
	}()
	err := Write(Poly{}, "/tmp/file")
	assert.EqualError(t, err, marshalIndentErr.Error())
}
