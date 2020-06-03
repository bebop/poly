package main

import "testing"

func BenchmarkParseGbk(b *testing.B) {
	for i := 0; i < b.N; i++ {
		ParseGbk("data/test.gbk")
	}
}

func BenchmarkParseGbk1(b *testing.B)     { BenchmarkParseGbk(b) }
func BenchmarkParseGbk10(b *testing.B)    { BenchmarkParseGbk(b) }
func BenchmarkParseGbk100(b *testing.B)   { BenchmarkParseGbk(b) }
func BenchmarkParseGbk1000(b *testing.B)  { BenchmarkParseGbk(b) }
func BenchmarkParseGbk10000(b *testing.B) { BenchmarkParseGbk(b) }
