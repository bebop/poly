package fix

import (
	"fmt"
	"log"
	"regexp"
	"strings"

	"github.com/TimothyStiles/poly/transform/variants"
)

var translator = buildPatternTranslator(variants.IUPAC2Bases())

func buildPatternTranslator(ambiguityCodes map[rune][]rune) map[rune]string {
	trans := make(map[rune]string)
	for r, options := range ambiguityCodes {
		if len(options) > 1 {
			trans[r] = "[" + string(options) + "]"
		}
	}
	trans['('] = "(?:"

	return trans
}

// TODO support (n/m) notation for point of cleavage for non-palindromic enzymes
func patternToRegexp(buf *strings.Builder, pattern string, translator map[rune]string) {
	buf.WriteRune('(')
	for _, r := range pattern {
		if match, found := translator[r]; found {
			buf.WriteString(match)
		} else {
			buf.WriteRune(r)
		}
	}
	buf.WriteRune(')')
}

func PatternsToRegexp(patterns []string, doubleStranded bool) (*regexp.Regexp, error) {
	var buf strings.Builder
	if len(patterns) == 0 {
		return nil, nil
	}

	// expand to reduce reallocations
	estSize := -1
	for _, pat := range patterns {
		estSize += len(pat) + 3
	}
	buf.Grow(estSize)

	for i, pat := range patterns {
		if i != 0 {
			buf.WriteRune('|')
		}
		patternToRegexp(&buf, pat, translator)
		if doubleStranded {
			reverseComplementPat := transform.ReverseComplement(pat)
			if reverseComplementPat != pat {
				buf.WriteRune('|')
				patternToRegexp(&buf, reverseComplementPat, translator)
			}
		}
	}
	// TODO: only print during testing
	log.Println(buf.String())
	return regexp.Compile(buf.String())
}

func HomopolymericRuns(length int) string {
	return fmt.Sprintf("N{%d}", length)
}
