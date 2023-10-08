package fix

import (
	"testing"

	"github.com/TimothyStiles/poly/random"
	"github.com/TimothyStiles/poly/transform/variants"
)

func Assertf(t *testing.T, cond bool, format string, args ...any) {
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
	Assertf(t, len(matches) == 4, "Expected 4 matches, got %d", len(matches))
}

func TestMultiple(t *testing.T) {
	re, err := patternsToRegexp([]string{"A", "C"})
	Assertf(t, err == nil, "Encountered error building regexp")
	matches := re.FindAllStringSubmatchIndex("AGCT", -1)
	Assertf(t, len(matches) == 2, "Expected 2 matches, got %d", len(matches))
}

func keys[K comparable, V any](m map[K]V) []K {
	array := make([]K, len(m))
	var i = 0
	for k := range m {
		array[i] = k
		i += 1
	}
	return array
}

func TestAcceptsRandomGen(t *testing.T) {
	allBases := string(keys[rune, []rune](variants.IUPAC2Bases()))
	re, err := patternsToRegexp([]string{allBases})
	Assertf(t, err == nil, "Encountered error building regexp")
	sample := random.DNASequenceFromPattern(allBases, 0)
	Assertf(t, len(re.FindAllStringSubmatchIndex(sample, -1)) == 1, "Expected sample DNA sequence %s to match pattern %s", sample, allBases)
}

func TestCheckHomopolymericRuns(t *testing.T) {
	testSeq := "GTAAAACGACGGCCAGT" // M13 fwd
	if run := checkHomopolymericRuns(testSeq); run == false {
		t.Errorf("checkHomopolymericRuns has changed on test. Got false instead of true")
	}

	testSeq = "ACGATGGCAGTAGCATGC" //"GTAAAACGACGGCCAGT" // M13 fwd

	if run := checkHomopolymericRuns(testSeq); run == true {
		t.Errorf("checkHomopolymericRuns has changed on test. Got true instead of false")
	}
}
