package fold_test

import (
	"fmt"

	"github.com/TimothyStiles/poly/fold"
)

func ExampleFold() {
	result, _ := fold.Fold("ACCCCCUCCUUCCUUGGAUCAAGGGGCUCAA", 37.0)
	brackets := result.DotBracket()
	fmt.Println(brackets)
	// Output: .((((.(((......)))....))))
}
