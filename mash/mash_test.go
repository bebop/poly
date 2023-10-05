package mash_test

import (
	"testing"

	"github.com/TimothyStiles/poly/mash"
)

func TestMash(t *testing.T) {
	fingerprint1 := mash.NewMash(17, 10)
	fingerprint1.Sketch("ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA")

	fingerprint2 := mash.NewMash(17, 9)
	fingerprint2.Sketch("ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA")

	distance := fingerprint1.Distance(fingerprint2)
	if distance != 0 {
		t.Errorf("Expected distance to be 0, got %f", distance)
	}

	distance = fingerprint2.Distance(fingerprint1)
	if distance != 0 {
		t.Errorf("Expected distance to be 0, got %f", distance)
	}

	spoofedFingerprint := mash.NewMash(17, 10)
	spoofedFingerprint.Sketches[0] = 0

	distance = fingerprint1.Distance(spoofedFingerprint)
	if distance != 1 {
		t.Errorf("Expected distance to be 1, got %f", distance)
	}

	spoofedFingerprint = mash.NewMash(17, 9)

	distance = fingerprint1.Distance(spoofedFingerprint)
	if distance != 1 {
		t.Errorf("Expected distance to be 1, got %f", distance)
	}
}
