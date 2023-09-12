package alphabet_test

import (
	"reflect"
	"testing"

	"github.com/TimothyStiles/poly/alphabet"
)

func TestAlphabet(t *testing.T) {
	symbols := []string{"A", "C", "G", "T"}
	a := alphabet.NewAlphabet(symbols)
	// Test encoding
	for i, symbol := range symbols {
		code, err := a.Encode(symbol)
		if err != nil {
			t.Errorf("Unexpected error encoding symbol %s: %v", symbol, err)
		}
		if code != i {
			t.Errorf("Incorrect encoding of symbol %s: expected %d, got %d", symbol, i, code)
		}
	}
	_, err := a.Encode("X")
	if err == nil {
		t.Error("Expected error for encoding symbol not in alphabet, but got nil")
	}

	// Test decoding
	for i, symbol := range symbols {
		decoded, err := a.Decode(i)
		if err != nil {
			t.Errorf("Unexpected error decoding code %d: %v", i, err)
		}
		if decoded != symbol {
			t.Errorf("Incorrect decoding of code %d: expected %s, got %s", i, symbol, decoded)
		}
	}
	_, err = a.Decode(len(symbols))
	if err == nil {
		t.Error("Expected error for decoding code not in alphabet, but got nil")
	}

	// Test extension
	extendedSymbols := []string{"N", "-", "*"}
	extendedAlphabet := a.Extend(extendedSymbols)
	for i, symbol := range symbols {
		code, err := extendedAlphabet.Encode(symbol)
		if err != nil {
			t.Errorf("Unexpected error encoding symbol %s: %v", symbol, err)
		}
		if code != i {
			t.Errorf("Incorrect encoding of symbol %s: expected %d, got %d", symbol, i, code)
		}
	}
	for i, symbol := range extendedSymbols {
		code, err := extendedAlphabet.Encode(symbol)
		if err != nil {
			t.Errorf("Unexpected error encoding symbol %s: %v", symbol, err)
		}
		if code != i+len(symbols) {
			t.Errorf("Incorrect encoding of symbol %s: expected %d, got %d", symbol, i+len(symbols), code)
		}
	}
}

func TestAlphabet_Symbols(t *testing.T) {
	// Test Symbols
	symbols := []string{"A", "C", "G", "T"}
	a := alphabet.NewAlphabet(symbols)
	if !reflect.DeepEqual(a.Symbols(), symbols) {
		t.Errorf("Symbols() = %v, want %v", a.Symbols(), symbols)
	}
}
