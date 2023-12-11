package fold_test

import (
	"fmt"

	"github.com/bebop/poly/fold"
)

func ExampleZuker() {
	result, _ := fold.Zuker("ACCCCCUCCUUCCUUGGAUCAAGGGGCUCAA", 37.0)
	brackets := result.DotBracket()
	fmt.Println(brackets)
	// Output: .((((.(((......)))....))))
}
