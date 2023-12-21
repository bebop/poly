package bwt_test

import (
	"fmt"
	"log"

	"github.com/bebop/poly/bwt"
	"golang.org/x/exp/slices"
)

func ExampleBWT_Count() {
	inputSequence := "AACCTGCCGTCGGGGCTGCCCGTCGCGGGACGTCGAAACGTGGGGCGAAACGTG"

	bwt, err := bwt.New(inputSequence)
	if err != nil {
		log.Fatal(err)
	}

	count, err := bwt.Count("CG")
	if err != nil {
		log.Fatal(err)
	}
	fmt.Println(count)
	// Output: 10
}

func ExampleBWT_Locate() {
	inputSequence := "AACCTGCCGTCGGGGCTGCCCGTCGCGGGACGTCGAAACGTGGGGCGAAACGTG"

	bwt, err := bwt.New(inputSequence)
	if err != nil {
		log.Fatal(err)
	}

	offsets, err := bwt.Locate("CG")
	if err != nil {
		log.Fatal(err)
	}
	slices.Sort(offsets)
	fmt.Println(offsets)
	// Output: [7 10 20 23 25 30 33 38 45 50]
}

func ExampleBWT_Extract() {
	inputSequence := "AACCTGCCGTCGGGGCTGCCCGTCGCGGGACGTCGAAACGTGGGGCGAAACGTG"

	bwt, err := bwt.New(inputSequence)
	if err != nil {
		log.Fatal(err)
	}

	extracted, err := bwt.Extract(48, 54)
	if err != nil {
		log.Fatal(err)
	}
	fmt.Println(extracted)
	// Output: AACGTG
}
