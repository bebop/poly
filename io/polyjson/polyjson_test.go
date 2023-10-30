package polyjson

import (
	"errors"
	"os"
	"strings"
	"testing"

	"github.com/TimothyStiles/poly/transform"
	"github.com/stretchr/testify/assert"
)

/******************************************************************************

JSON related tests begin here.

******************************************************************************/

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

	feature.Description = "Green Fluorescent Protein"
	feature.Location.SubLocations = []Location{subLocation, subLocationReverseComplemented}

	// Add the GFP feature to the sequence struct.
	_ = sequence.AddFeature(&feature)

	// get the GFP feature sequence string from the sequence struct.
	featureSequence, _ := feature.GetSequence()

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
	_, err := Parse(strings.NewReader(""))
	assert.EqualError(t, err, unmarshalErr.Error())
}

func TestRead_error(t *testing.T) {
	readErr := errors.New("read error")
	oldReadFileFn := readFileFn
	readFileFn = func(filename string) (*os.File, error) {
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

func TestConvert(t *testing.T) {
	sequence, err := Read("../../data/cat.json")
	assert.NoError(t, err)
	sequence.ToGenbank()
}
