package finder

import (
	"fmt"
)

func ExampleForbiddenSequencesFinder() {
	sequence := "AAAAAATCGGTCGTAAAAAATT"
	var functions []func(string) []Match
	functions = append(functions, ForbiddenSequence([]string{"AAAAAA"}))

	problems := Find(sequence, functions)
	fmt.Println(problems)

	// Output: [{0 6 Forbidden sequence: AAAAAA} {14 20 Forbidden sequence: AAAAAA}]
}
