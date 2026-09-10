package gff

import (
	"bytes"
	"errors"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/google/go-cmp/cmp"
	"github.com/google/go-cmp/cmp/cmpopts"
	"github.com/pmezard/go-difflib/difflib"
)

/******************************************************************************

Gff related tests and benchmarks begin here.

******************************************************************************/

// TODO should delete output files.

func TestGffIO(t *testing.T) {
	tmpDataDir, err := os.MkdirTemp("", "data-*")
	if err != nil {
		t.Error(err)
	}
	defer os.RemoveAll(tmpDataDir)

	testInputPath := "../../data/ecoli-mg1655-short.gff"
	tmpGffFilePath := filepath.Join(tmpDataDir, "ecoli-mg1655-short.gff")

	testSequence, _ := Read(testInputPath)
	_ = Write(testSequence, tmpGffFilePath)

	readTestSequence, _ := Read(tmpGffFilePath)

	if diff := cmp.Diff(testSequence, readTestSequence, cmpopts.IgnoreFields(Feature{}, "ParentSequence")); diff != "" {
		t.Errorf("Parsing the output of Build() does not produce the same output as parsing the original file read with ReadGff(). Got this diff:\n%s", diff)
	}

	original, _ := os.ReadFile(testInputPath)
	builtOutput, _ := os.ReadFile(tmpGffFilePath)
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

// testing that readAllFn() returns an error.
func TestParseReader_error(t *testing.T) {
	parseErr := errors.New("parse error")
	oldReadAllFn := readAllFn
	readAllFn = func(r io.Reader) ([]byte, error) {
		return nil, parseErr
	}
	defer func() {
		readAllFn = oldReadAllFn
	}()
	_, err := Parse(strings.NewReader(""))
	if err != parseErr {
		t.Errorf("Parse() did not return the expected error. Got %v, expected %v", err, parseErr)
	}
}

// testing that all Atoi() calls return an error.
func TestParseAtoi_error(t *testing.T) {
	file, _ := openFn("../../data/ecoli-mg1655-short.gff")
	fileBytes, _ := readAllFn(file)
	fileString := string(fileBytes)
	recievedAtoiInputs := []string{}

	parseErr := errors.New("parse error")
	oldAtoiFn := atoiFn
	// helper function to see if the input string for atoiFn is in the recievedAtoiInputs slice.
	contains := func(stringSlice []string, expression string) bool {
		for _, input := range stringSlice {
			if input == expression {
				return true
			}
		}
		return false
	}
	atoiFn = func(gffString string) (int, error) {
		if !contains(recievedAtoiInputs, gffString) {
			recievedAtoiInputs = append(recievedAtoiInputs, gffString)
			// length = len(recievedAtoiInputs)
			return 0, parseErr
		}
		return oldAtoiFn(gffString)
	}

	defer func() {
		atoiFn = oldAtoiFn
	}()
	for index := 0; index <= 10; index++ {
		var gffBuffer bytes.Buffer
		gffBuffer.WriteString(fileString)
		_, _ = Parse(&gffBuffer)
	}
}

// testing that Read can return an appropriate error.
func TestRead_error(t *testing.T) {
	readErr := errors.New("open : no such file or directory")
	oldOpenFn := openFn
	openFn = func(filepath string) (*os.File, error) {
		return nil, readErr
	}
	defer func() {
		openFn = oldOpenFn
	}()
	_, err := Read("../../data/ecoli-mg1655-short.gff")
	if err != readErr {
		t.Errorf("Read() did not return the expected error. Got %v, expected %v", err, readErr)
	}
}

// TestAddFeature_doesNotShareState verifies that AddFeature stores an
// independent copy of the given feature: mutating the caller's Feature
// (its Attributes map or Location.SubLocations slice) after the call must
// not affect the feature stored on the Gff, and vice versa.
func TestAddFeature_doesNotShareState(t *testing.T) {
	var sequence Gff

	feature := Feature{
		Name: "testFeature",
		Attributes: map[string]string{
			"gene_id": "original",
		},
		Location: Location{
			Start: 0,
			End:   10,
			SubLocations: []Location{
				{Start: 0, End: 5},
			},
		},
	}

	if err := sequence.AddFeature(&feature); err != nil {
		t.Fatalf("AddFeature() returned an unexpected error: %v", err)
	}

	// Mutating the caller's copy after AddFeature must not leak into the
	// feature stored on sequence.Features.
	feature.Attributes["gene_id"] = "mutated"
	feature.Location.SubLocations[0].Start = 999

	storedFeature := sequence.Features[0]
	if storedFeature.Attributes["gene_id"] != "original" {
		t.Errorf("mutating the caller's Attributes map after AddFeature changed the stored feature's Attributes: got %q, want %q", storedFeature.Attributes["gene_id"], "original")
	}
	if storedFeature.Location.SubLocations[0].Start != 0 {
		t.Errorf("mutating the caller's SubLocations after AddFeature changed the stored feature's SubLocations: got %d, want %d", storedFeature.Location.SubLocations[0].Start, 0)
	}

	// Mutating the stored feature must not leak back into the caller's
	// feature either.
	sequence.Features[0].Attributes["gene_id"] = "mutated again"
	sequence.Features[0].Location.SubLocations[0].Start = -1

	if feature.Attributes["gene_id"] != "mutated" {
		t.Errorf("mutating the stored feature's Attributes changed the caller's Attributes: got %q, want %q", feature.Attributes["gene_id"], "mutated")
	}
	if feature.Location.SubLocations[0].Start != 999 {
		t.Errorf("mutating the stored feature's SubLocations changed the caller's SubLocations: got %d, want %d", feature.Location.SubLocations[0].Start, 999)
	}

	// The stored feature's ParentSequence should point back at sequence.
	if storedFeature.ParentSequence != &sequence {
		t.Errorf("stored feature's ParentSequence = %p, want %p", storedFeature.ParentSequence, &sequence)
	}
}

func BenchmarkReadGff(b *testing.B) {
	for i := 0; i < b.N; i++ {
		_, _ = Read("../../data/ecoli-mg1655-short.gff")
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
