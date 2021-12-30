package gff

import (
	"fmt"
	"io/ioutil"
	"os"
	"path/filepath"
	"testing"

	"github.com/TimothyStiles/poly/io/poly"
	"github.com/google/go-cmp/cmp"
	"github.com/google/go-cmp/cmp/cmpopts"
	"github.com/pmezard/go-difflib/difflib"
)

/******************************************************************************

Gff related tests and benchmarks begin here.

******************************************************************************/

func ExampleRead() {

	sequence := Read("../../data/ecoli-mg1655-short.gff")
	fmt.Println(sequence.Meta.Name)
	// Output: U00096.3
}

func ExampleParse() {
	file, _ := ioutil.ReadFile("../../data/ecoli-mg1655-short.gff")
	sequence := Parse(file)

	fmt.Println(sequence.Meta.Name)
	// Output: U00096.3
}

func ExampleBuild() {

	sequence := Read("../../data/ecoli-mg1655-short.gff")
	gffBytes := Build(sequence)
	reparsedSequence := Parse(gffBytes)

	fmt.Println(reparsedSequence.Meta.Name)
	// Output: U00096.3

}

func ExampleWrite() {
	tmpDataDir, err := ioutil.TempDir("", "data-*")
	if err != nil {
		fmt.Println(err.Error())
	}
	defer os.RemoveAll(tmpDataDir)

	sequence := Read("../../data/ecoli-mg1655-short.gff")

	tmpGffFilePath := filepath.Join(tmpDataDir, "ecoli-mg1655-short.gff")
	Write(sequence, tmpGffFilePath)

	testSequence := Read(tmpGffFilePath)

	fmt.Println(testSequence.Meta.Name)
	// Output: U00096.3
}

// TODO should delete output files.

func TestGffIO(t *testing.T) {
	tmpDataDir, err := ioutil.TempDir("", "data-*")
	if err != nil {
		t.Error(err)
	}
	defer os.RemoveAll(tmpDataDir)

	testInputPath := "../../data/ecoli-mg1655-short.gff"
	tmpGffFilePath := filepath.Join(tmpDataDir, "ecoli-mg1655-short.gff")

	testSequence := Read(testInputPath)
	Write(testSequence, tmpGffFilePath)

	readTestSequence := Read(tmpGffFilePath)

	if diff := cmp.Diff(testSequence, readTestSequence, cmpopts.IgnoreFields(poly.Feature{}, "ParentSequence")); diff != "" {
		t.Errorf("Parsing the output of Build() does not produce the same output as parsing the original file read with ReadGff(). Got this diff:\n%s", diff)
	}

	original, _ := ioutil.ReadFile(testInputPath)
	builtOutput, _ := ioutil.ReadFile(tmpGffFilePath)
	gffDiff := difflib.UnifiedDiff{
		A:        difflib.SplitLines(string(original)),
		B:        difflib.SplitLines(string(builtOutput)),
		FromFile: testInputPath,
		ToFile:   tmpGffFilePath,
		Context:  3,
	}

	gffDiffText, _ := difflib.GetUnifiedDiffString(gffDiff)

	if gffDiffText != "" {
		t.Errorf("Build() does not output the same file as was input through ReadGff(). Got this diff:\n%s", gffDiffText)
	}

}

func BenchmarkReadGff(b *testing.B) {
	for i := 0; i < b.N; i++ {
		Read("../../data/ecoli-mg1655-short.gff")
	}
}

func BenchmarkReadGff1(b *testing.B)     { BenchmarkReadGff(b) }
func BenchmarkReadGff10(b *testing.B)    { BenchmarkReadGff(b) }
func BenchmarkReadGff100(b *testing.B)   { BenchmarkReadGff(b) }
func BenchmarkReadGff1000(b *testing.B)  { BenchmarkReadGff(b) }
func BenchmarkReadGff10000(b *testing.B) { BenchmarkReadGff(b) }

/******************************************************************************

Gff related tests and benchmarks end here.

******************************************************************************/
