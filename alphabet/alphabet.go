/*
Package alphabet provides structs for defining biological sequence alphabets.
*/
package alphabet

// TODO: add Alphabet for codons

import "fmt"

// Alphabet is a struct that holds a list of symbols and a map of symbols to their index in the list.
type Alphabet struct {
	symbols  []string
	encoding map[interface{}]uint8
}

// Error is an error type that is returned when a symbol is not in the alphabet.
type Error struct {
	message string
}

// Error returns the error message for AlphabetError.
func (e *Error) Error() string {
	return e.message
}

// NewAlphabet creates a new alphabet from a list of symbols.
func NewAlphabet(symbols []string) *Alphabet {
	encoding := make(map[interface{}]uint8)
	for index, symbol := range symbols {
		encoding[symbol] = uint8(index)
		encoding[index] = uint8(index)
	}
	return &Alphabet{symbols, encoding}
}

// Encode returns the index of a symbol in the alphabet.
func (alphabet *Alphabet) Encode(symbol interface{}) (uint8, error) {
	c, ok := alphabet.encoding[symbol]
	if !ok {
		return 0, fmt.Errorf("Symbol %v not in alphabet", symbol)
	}
	return c, nil
}

// TODO: compress more when len(symbols) << 2^8
// TODO: DecodeAll
func (alphabet *Alphabet) EncodeAll(seq string) ([]uint8, error) {
	encoded := make([]uint8, len(seq))
	for i, r := range seq {
		encoding, err := alphabet.Encode(string(r))
		if err != nil {
			return nil, fmt.Errorf("Symbol %c in position %d not in alphabet", r, i)
		}
		encoded[i] = uint8(encoding)
	}
	return encoded, nil
}

func (alphabet *Alphabet) Check(seq string) int {
	for i, r := range seq {
		encoding, err := alphabet.Encode(string(r))
		if err != nil {
			return i
		}
	}
	return -1
}

// Decode returns the symbol at a given index in the alphabet.
func (alphabet *Alphabet) Decode(code interface{}) (string, error) {
	c, ok := code.(int)
	if !ok || c < 0 || c >= len(alphabet.symbols) {
		return "", &Error{fmt.Sprintf("Code %v not in alphabet", code)}
	}
	return alphabet.symbols[c], nil
}

// Extend returns a new alphabet that is the original alphabet extended with a list of symbols.
func (alphabet *Alphabet) Extend(symbols []string) *Alphabet {
	extended := append(alphabet.symbols, symbols...)
	return NewAlphabet(extended)
}

// Symbols returns the list of symbols in the alphabet.
func (alphabet *Alphabet) Symbols() []string {
	return alphabet.symbols
}

var (
	DNA     = NewAlphabet([]string{"A", "C", "G", "T"})
	RNA     = NewAlphabet([]string{"A", "C", "G", "U"})
	Protein = NewAlphabet([]string{"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"})
)
