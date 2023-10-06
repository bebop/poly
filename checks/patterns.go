package checks

import (
	"regexp"
	"strings"
	"bytes"
)

var iupacToBases = map[string][]string{
	"A": {"A"},
	"C": {"C"},
	"G": {"G"},
	"T": {"T"},
	"R": {"A", "G"},
	"Y": {"C", "T"},
	"S": {"G", "C"},
	"W": {"A", "T"},
	"K": {"G", "T"},
	"M": {"A", "C"},
	"B": {"C", "G", "T"},
	"D": {"A", "G", "T"},
	"H": {"A", "C", "T"},
	"V": {"A", "C", "G"},
	"N": {"A", "C", "G", "T"},
}

var iupacToComplements = map[string]string{
	"A": "T",
	"C": "G",
	"G": "C",
	"T": "A",
	"R": "Y",
	"Y": "R",
	"S": "S",
	"W": "W",
	"K": "M",
	"M": "K",
	"B": "V",
	"D": "H",
	"H": "D",
	"V": "B",
	"N": "N",
}

// TODO allow passthru of special characters
func patternToRegexp(buf bytes.Buffer, pattern string) {
	buf.WriteString("(")
	for _, r := range pattern {
		bases := iupacToBases[string(r)]
		if len(bases) == 1 {
			buf.WriteString(bases[0])
		} else {
			buf.WriteString("[")
			buf.WriteString(strings.Join(bases, ""))
			buf.WriteString("]")
		}
	}
	buf.WriteString(")")
}

func patternsToRegexp(patterns []string) (*regexp.Regexp, error)  {
	var buf bytes.Buffer
	for i, pat := range patterns {
		if i != 0 {
			buf.WriteString("|")
		}
		patternToRegexp(buf, pat)
	}
	return regexp.Compile(buf.String())
}
