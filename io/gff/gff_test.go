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
