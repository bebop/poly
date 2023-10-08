package alphabet

import (
	"strings"
	"testing"

	"github.com/stretchr/testify/assert"

	"github.com/TimothyStiles/poly/transform/variants"
)

func Assertf(t *testing.T, cond bool, format string, args ...any) {
	if !cond {
		t.Errorf(format, args...)
	}
}

func TestKmerCounterEmpty(t *testing.T) {
	kc := NewKmerCounter(DNA, 3)

	count, err := LookupCount(kc, "A")
	assert.EqualValues(t, 0, count)
	assert.NoError(t, err)
	count, err = LookupCount(kc, "X")
	assert.Error(t, err)
}

func TestKmerCounterRepeated(t *testing.T) {
	kc := NewKmerCounter(DNA, 3)
	seq := strings.Repeat("A", 12)
	assert.Equal(t, 12, len(seq))
	err := Observe(kc, seq)
	assert.NoError(t, err)

	assert.Equal(t, 12, int(kc.total))

	onemers, _ := variants.AllVariantsIUPAC("N")
	for _, kmer := range onemers {
		count, err := LookupCount(kc, kmer)
		assert.NoError(t, err)
		if kmer == "A" {
			assert.EqualValues(t, 12, count, "Wrong count")
		} else {
			assert.EqualValues(t, 0, count, "Reports nonzero count for 1mer not present in sequence %v", kmer)
		}
	}
	twomers, _ := variants.AllVariantsIUPAC("NN")
	for _, kmer := range twomers {
		count, err := LookupCount(kc, kmer)
		assert.NoError(t, err)
		if kmer == "AA" {
			assert.EqualValues(t, 11, count, "Wrong count")
		} else {
			assert.EqualValues(t, 0, count, "Reports nonzero count for 2mer not present in sequence %v", kmer)
		}
	}
	threemers, _ := variants.AllVariantsIUPAC("NNN")
	for _, kmer := range threemers {
		count, err := LookupCount(kc, kmer)
		assert.NoError(t, err)
		if kmer == "AAA" {
			assert.EqualValues(t, 10, count, "Wrong count")
		} else {
			assert.EqualValues(t, 0, count, "Reports nonzero count for 3mer not present in sequence %v", kmer)
		}
	}
}
