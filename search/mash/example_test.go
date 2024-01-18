package mash_test

import (
	"fmt"

	"github.com/bebop/poly/search/mash"
)

func ExampleMash() {
	fingerprint1 := mash.New(17, 10)
	fingerprint1.Sketch("ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA")

	fingerprint2 := mash.New(17, 9)
	fingerprint2.Sketch("ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA")

	distance := fingerprint1.Distance(fingerprint2)

	fmt.Println(distance)

	// Output:
	// 0
}
