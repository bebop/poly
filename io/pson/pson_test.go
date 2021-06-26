package pson

import (
	"fmt"
	"io/ioutil"
	"os"
	"path/filepath"
	"testing"

	"github.com/TimothyStiles/poly"
	"github.com/TimothyStiles/poly/io/genbank"
	"github.com/TimothyStiles/poly/io/gff"

	"github.com/google/go-cmp/cmp"
	"github.com/google/go-cmp/cmp/cmpopts"
)

/******************************************************************************

JSON related tests begin here.

******************************************************************************/

func ExampleRead() {
	sequence := Read("../../data/sample.json")

	fmt.Println(sequence.Meta.Source)
	//output: Saccharomyces cerevisiae (baker's yeast)
}

func ExampleParse() {
	file, _ := ioutil.ReadFile("../../data/sample.json")
	sequence := Parse(file)

	fmt.Println(sequence.Meta.Source)
	//output: Saccharomyces cerevisiae (baker's yeast)
}

func ExampleWrite() {
	tmpDataDir, err := ioutil.TempDir("", "data-*")
	if err != nil {
		fmt.Println(err.Error())
	}
	defer os.RemoveAll(tmpDataDir)

	sequence := Read("../../data/sample.json")

	tmpJSONFilePath := filepath.Join(tmpDataDir, "sample.json")
	Write(sequence, tmpJSONFilePath)

	testSequence := Read(tmpJSONFilePath)

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
	Write(testSequence, tmpJSONFilePath)

	readTestSequence := Read(tmpJSONFilePath)

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
	Write(gffTestSequence, tmpJSONFilePath)

	gffReadTestSequence := Read(tmpJSONFilePath)

	if diff := cmp.Diff(gffTestSequence, gffReadTestSequence, cmpopts.IgnoreFields(poly.Feature{}, "ParentSequence")); diff != "" {
		t.Errorf(" mismatch (-want +got):\n%s", diff)
	}

}

/******************************************************************************

JSON related tests end here.

******************************************************************************/
