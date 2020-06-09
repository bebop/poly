package main

import (
	"testing"

	"github.com/google/go-cmp/cmp"
)

func BenchmarkReadGff(b *testing.B) {
	for i := 0; i < b.N; i++ {
		ParseGff("data/ecoli-mg1655.gff")
	}
}

func TestGffIO(t *testing.T) {
	testSequence := ReadGff("data/ecoli-mg1655.gff")
	WriteGff(testSequence, "data/test.gff")
	WriteJSON(testSequence, "data/testgff.json")
	readTestSequence := ReadGff("data/test.gff")
	if diff := cmp.Diff(testSequence, readTestSequence); diff != "" {
		t.Errorf("MakeGatewayInfo() mismatch (-want +got):\n%s", diff)
	}
}

func BenchmarkReadGff1(b *testing.B)     { BenchmarkReadGff(b) }
func BenchmarkReadGff10(b *testing.B)    { BenchmarkReadGff(b) }
func BenchmarkReadGff100(b *testing.B)   { BenchmarkReadGff(b) }
func BenchmarkReadGff1000(b *testing.B)  { BenchmarkReadGff(b) }
func BenchmarkReadGff10000(b *testing.B) { BenchmarkReadGff(b) }
