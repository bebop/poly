package alphabet

import "fmt"

type Alphabet struct {
	symbols  []string
	encoding map[interface{}]int
}

type AlphabetError struct {
	message string
}

func (e *AlphabetError) Error() string {
	return e.message
}

func NewAlphabet(symbols []string) *Alphabet {
	encoding := make(map[interface{}]int)
	for i, symbol := range symbols {
		encoding[symbol] = i
		encoding[i] = i
	}
	return &Alphabet{symbols, encoding}
}

func (a *Alphabet) Encode(symbol interface{}) (int, error) {
	code, ok := a.encoding[symbol]
	if !ok {
		return 0, &AlphabetError{fmt.Sprintf("Symbol %v not in alphabet", symbol)}
	}
	return code, nil
}

func (a *Alphabet) Decode(code interface{}) (string, error) {
	c, ok := code.(int)
	if !ok || c < 0 || c >= len(a.symbols) {
		return "", &AlphabetError{fmt.Sprintf("Code %v not in alphabet", code)}
	}
	return a.symbols[c], nil
}

func (a *Alphabet) Extend(symbols []string) *Alphabet {
	extended := append(a.symbols, symbols...)
	return NewAlphabet(extended)
}

func (a *Alphabet) Symbols() []string {
	return a.symbols
}
