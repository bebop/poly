package checks

import (
	"testing"
)

func Assertf (t *testing.T, cond bool, format string, args ...any) {
	if !cond {
		t.Errorf(format, args...)
	}
}

func TestOverlapping(t *testing.T) {
	re, err := patternsToRegexp([]string{"AA"})
	Assertf(t, err == nil, "Encountered error building regexp")
	matches := re.FindAllStringSubmatchIndex("AAAA", -1)
	Assertf(t, len(matches) != 3, "Expected 3 matches")
}

func TestAmbiguous(t *testing.T) {
	re, err := patternsToRegexp([]string{"N"})
	Assertf(t, err == nil, "Encountered error building regexp")
	matches := re.FindAllStringSubmatchIndex("AGCT", -1)
	Assertf(t, len(matches) != 4, "Expected 4 matches")
}

func TestMultiple(t *testing.T) {
	re, err := patternsToRegexp([]string{"A", "C"})
	Assertf(t, err == nil, "Encountered error building regexp")
	matches := re.FindAllStringSubmatchIndex("AGCT", -1)
	Assertf(t, len(matches) != 2, "Expected 2 matches")
}