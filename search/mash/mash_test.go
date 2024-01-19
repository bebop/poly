package mash_test

import (
	"testing"

	"github.com/bebop/poly/search/mash"
)

func TestMash(t *testing.T) {
	fingerprint1 := mash.New(17, 10)
	fingerprint1.Sketch("ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA")

	fingerprint2 := mash.New(17, 9)
	fingerprint2.Sketch("ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA")

	distance := fingerprint1.Distance(fingerprint2)
	if distance != 0 {
		t.Errorf("Expected distance to be 0, got %f", distance)
	}

	distance = fingerprint2.Distance(fingerprint1)
	if distance != 0 {
		t.Errorf("Expected distance to be 0, got %f", distance)
	}

	spoofedFingerprint := mash.New(17, 10)
	spoofedFingerprint.Sketches[0] = 0

	distance = fingerprint1.Distance(spoofedFingerprint)
	if distance != 1 {
		t.Errorf("Expected distance to be 1, got %f", distance)
	}

	spoofedFingerprint = mash.New(17, 9)

	distance = fingerprint1.Distance(spoofedFingerprint)
	if distance != 1 {
		t.Errorf("Expected distance to be 1, got %f", distance)
	}

	fingerprint1 = mash.New(17, 10)
	fingerprint1.Sketch("ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA")

	fingerprint2 = mash.New(17, 5)
	fingerprint2.Sketch("ATCGATCGATCGATCGATCGATCGATCGATCGATCGAATGCGATCGATCGATCGATCGATCG")

	distance = fingerprint1.Distance(fingerprint2)
	if !(distance > 0.19 && distance < 0.21) {
		t.Errorf("Expected distance to be 0.19999999999999996, got %f", distance)
	}

	fingerprint1 = mash.New(17, 10)
	fingerprint1.Sketch("ATCGATCGATCGATCGATCGATCGATCGATCGATCGAATGCGATCGATCGATCGATCGATCG")

	fingerprint2 = mash.New(17, 5)
	fingerprint2.Sketch("ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA")

	distance = fingerprint1.Distance(fingerprint2)
	if distance != 0 {
		t.Errorf("Expected distance to be 0, got %f", distance)
	}
}

func BenchmarkMashDistancee(b *testing.B) {
	fingerprint1 := mash.New(17, 10)
	fingerprint1.Sketch("ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA")

	fingerprint2 := mash.New(17, 9)
	fingerprint2.Sketch("ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA")

	for i := 0; i < b.N; i++ {
		fingerprint1.Distance(fingerprint2)
	}
}
