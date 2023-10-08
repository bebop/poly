package checks

import (
	"encoding/json"
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

	b, err := json.Marshal(trans)
	if err != nil {
		log.Fatal(err)
	}
	log.Println(string(b))

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

func patternsToRegexp(patterns []string) (*regexp.Regexp, error) {
	var buf strings.Builder
	if len(patterns) == 0 {
		return nil, nil
	}
	estSize := -1
	for _, pat := range patterns {
		estSize += len(pat) + 3
	}
	buf.Grow(estSize)
	for i, pat := range patterns {
		if i != 0 {
			buf.WriteString("|")
		}
		patternToRegexp(&buf, pat, translator)
	}
	log.Println(buf.String())
	return regexp.Compile(buf.String())
}
