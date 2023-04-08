package fold_test

import (
	"fmt"

	"github.com/TimothyStiles/poly/fold"
)

func ExampleFold() {
	structs, _ := fold.Fold("ACCCCCUCCUUCCUUGGAUCAAGGGGCUCAA", 37.0)
	brackets := fold.DotBracket(structs)
	fmt.Println(brackets)
	// Output: .((((.(((......)))....))))
}
