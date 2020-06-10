package main

import (
	"io/ioutil"
	"testing"

	"github.com/google/go-cmp/cmp"
	"github.com/pmezard/go-difflib/difflib"
)

func BenchmarkReadGff(b *testing.B) {
	for i := 0; i < b.N; i++ {
		ParseGff("data/ecoli-mg1655.gff")
	}
}

func TestGffIO(t *testing.T) {
	testSequence := ReadGff("data/ecoli-mg1655.gff")
	WriteGff(testSequence, "data/test.gff")
	readTestSequence := ReadGff("data/test.gff")
	if diff := cmp.Diff(testSequence, readTestSequence); diff != "" {
		t.Errorf("Parsing the output of BuildGff() does not produce the same output as parsing the original file read with ReadGff(). Got this diff:\n%s", diff)
	}

	original, _ := ioutil.ReadFile("data/ecoli-mg1655.gff")
	builtOutput, _ := ioutil.ReadFile("data/test.gff")
	gffDiff := difflib.UnifiedDiff{
		A:        difflib.SplitLines(string(original)),
		B:        difflib.SplitLines(string(builtOutput)),
		FromFile: "data/ecoli-mg1655.gff",
		ToFile:   "data/test.gff",
		Context:  3,
	}

	gffDiffText, _ := difflib.GetUnifiedDiffString(gffDiff)
	if gffDiffText != "" {
		t.Errorf("BuildGff() does not output the same file as was input through ReadGff(). Got this diff:\n%s", gffDiffText)
	}
	


}

func BenchmarkReadGff1(b *testing.B)     { BenchmarkReadGff(b) }
func BenchmarkReadGff10(b *testing.B)    { BenchmarkReadGff(b) }
func BenchmarkReadGff100(b *testing.B)   { BenchmarkReadGff(b) }
func BenchmarkReadGff1000(b *testing.B)  { BenchmarkReadGff(b) }
func BenchmarkReadGff10000(b *testing.B) { BenchmarkReadGff(b) }
