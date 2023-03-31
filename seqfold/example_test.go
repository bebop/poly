package seqfold_test

import (
	"fmt"

	"github.com/TimothyStiles/poly/seqfold"
)

func ExampleFold() {
	structs, _ := seqfold.Fold("ACCCCCUCCUUCCUUGGAUCAAGGGGCUCAA", 37.0)
	brackets := seqfold.DotBracket(structs)
	fmt.Println(brackets)
	// Output: .((((.(((......)))....))))
}
