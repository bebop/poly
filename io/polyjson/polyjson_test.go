package polyjson_test

import (
	"fmt"
	"io/ioutil"
	"os"
	"path/filepath"
	"testing"
	"time"

	"github.com/TimothyStiles/poly/io/genbank"
	"github.com/TimothyStiles/poly/io/gff"
	"github.com/TimothyStiles/poly/io/poly"
	"github.com/TimothyStiles/poly/io/polyjson"
	"github.com/TimothyStiles/poly/seqhash"

	"github.com/google/go-cmp/cmp"
	"github.com/google/go-cmp/cmp/cmpopts"
)

/******************************************************************************

JSON related tests begin here.

******************************************************************************/

func Example() {
	//initiate a new polyjson sequence struct
	var sequence polyjson.Poly

	// Define the meta section of our sequence.
	sequence.Meta.Name = "Cat DNA"
	sequence.Meta.Description = "Synthetic Cat DNA for testing purposes."
	sequence.Meta.CreatedBy = "Catz (all you basepair are belong to us)"
	sequence.Meta.CreatedOn = time.Now()
	sequence.Meta.URL.Host = "www.catz.com"
	sequence.Meta.URL.Path = "/catz"
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

}

func ExampleRead() {
	sequence := polyjson.Read("../../data/sample.json")

	fmt.Println(sequence.Meta)
	//output: Saccharomyces cerevisiae (baker's yeast)
}

func ExampleParse() {
	file, _ := ioutil.ReadFile("../../data/sample.json")
	sequence := polyjson.Parse(file)

	fmt.Println(sequence.Meta.Source)
	//output: Saccharomyces cerevisiae (baker's yeast)
}

func ExampleWrite() {
	tmpDataDir, err := ioutil.TempDir("", "data-*")
	if err != nil {
		fmt.Println(err.Error())
	}
	defer os.RemoveAll(tmpDataDir)

	sequence := polyjson.Read("../../data/sample.json")

	tmpJSONFilePath := filepath.Join(tmpDataDir, "sample.json")
	polyjson.Write(sequence, tmpJSONFilePath)

	testSequence := polyjson.Read(tmpJSONFilePath)

	fmt.Println(testSequence.Meta.Source)
	//output: Saccharomyces cerevisiae (baker's yeast)
}

func TestGbkToJSON(t *testing.T) {
	tmpDataDir, err := ioutil.TempDir("", "data-*")
	if err != nil {
		t.Error(err)
	}
	defer os.RemoveAll(tmpDataDir)

	testSequence := genbank.Read("../../data/puc19.gbk")

	tmpJSONFilePath := filepath.Join(tmpDataDir, "puc19.json")
	polyjson.Write(testSequence, tmpJSONFilePath)

	readTestSequence := polyjson.Read(tmpJSONFilePath)

	if diff := cmp.Diff(testSequence, readTestSequence, cmpopts.IgnoreFields(poly.Feature{}, "ParentSequence")); diff != "" {
		t.Errorf(" mismatch (-want +got):\n%s", diff)
	}
}

func TestGffToJSON(t *testing.T) {
	tmpDataDir, err := ioutil.TempDir("", "data-*")
	if err != nil {
		t.Error(err)
	}
	defer os.RemoveAll(tmpDataDir)

	gffTestSequence := gff.Read("../../data/ecoli-mg1655-short.gff")

	tmpJSONFilePath := filepath.Join(tmpDataDir, "ecoli-mg1655-short.json")
	polyjson.Write(gffTestSequence, tmpJSONFilePath)

	gffReadTestSequence := polyjson.Read(tmpJSONFilePath)

	if diff := cmp.Diff(gffTestSequence, gffReadTestSequence, cmpopts.IgnoreFields(poly.Feature{}, "ParentSequence")); diff != "" {
		t.Errorf(" mismatch (-want +got):\n%s", diff)
	}

}

/******************************************************************************

JSON related tests end here.

******************************************************************************/
