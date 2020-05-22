package main

import "testing"

func BenchmarkParseGff(b *testing.B) {
	for i := 0; i < b.N; i++ {
		parseGff("data/ecoli-mg1655.gff")
	}
}

func BenchmarkParseGff1(b *testing.B)     { BenchmarkParseGbk(b) }
func BenchmarkParseGff10(b *testing.B)    { BenchmarkParseGbk(b) }
func BenchmarkParseGff100(b *testing.B)   { BenchmarkParseGbk(b) }
func BenchmarkParseGff1000(b *testing.B)  { BenchmarkParseGbk(b) }
func BenchmarkParseGff10000(b *testing.B) { BenchmarkParseGbk(b) }
