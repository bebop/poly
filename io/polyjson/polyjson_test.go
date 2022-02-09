package polyjson_test

import (
	"fmt"
	"io/ioutil"
	"os"
	"path/filepath"
	"time"

	"github.com/TimothyStiles/poly/io/polyjson"
	"github.com/TimothyStiles/poly/seqhash"
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
	var sequence polyjson.Poly

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
	catFeature := polyjson.Feature{}
	catFeature.Name = "Cat coding region."
	catFeature.Description = "a cat coding region at the beginning of our sequence."
	catFeature.Type = "CDS"
	catFeature.Location.Start = 0
	catFeature.Location.End = 8
	catFeature.Tags = map[string]string{"product": "cat protein"}

	sequence.AddFeature(&catFeature) // add the feature annotation to our sequence

	// write our sequence to a JSON file
	tmpJSONFilePath := filepath.Join(tmpDataDir, "sample.json")
	polyjson.Write(sequence, tmpJSONFilePath)

	exportedSequence, _ := polyjson.Read(tmpJSONFilePath)

	// print our struct DNA sequence
	fmt.Println(exportedSequence.Sequence)
	// Output: CATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCAT
}

func ExampleRead() {
	sequence, _ := polyjson.Read("../../data/cat.json")

	fmt.Println(sequence.Sequence)
	// Output: CATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCAT
}

func ExampleParse() {
	file, _ := ioutil.ReadFile("../../data/cat.json")
	sequence, _ := polyjson.Parse(file)

	fmt.Println(sequence.Sequence)
	// Output: CATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCAT
}

func ExampleWrite() {
	tmpDataDir, err := ioutil.TempDir("", "data-*")
	if err != nil {
		fmt.Println(err.Error())
	}
	defer os.RemoveAll(tmpDataDir)

	sequence, _ := polyjson.Read("../../data/cat.json")

	tmpJSONFilePath := filepath.Join(tmpDataDir, "sample.json")
	polyjson.Write(sequence, tmpJSONFilePath)

	testSequence, _ := polyjson.Read(tmpJSONFilePath)

	fmt.Println(testSequence.Sequence)
	// Output: CATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCAT
}

/******************************************************************************

JSON related tests end here.

******************************************************************************/
