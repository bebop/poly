package mash_test

import (
	"testing"

	"github.com/TimothyStiles/poly/mash"
)

func TestMash(t *testing.T) {
	fingerprint1 := mash.NewMash(17, 10)
	fingerprint1.Sketch("ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA")

	fingerprint2 := mash.NewMash(17, 10)
	fingerprint2.Sketch("ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA")

	distance := fingerprint1.Distance(fingerprint2)
	if distance != 0 {
		t.Errorf("Expected distance to be 0, got %f", distance)
	}
}
