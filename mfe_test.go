package poly

import (
	"fmt"
)

func ExampleCalculateMfe() {
	// 2 3 1 4 2 1
	// A C G A U C
	// 1 2 3 1 4 2 1 3 1 3 1 4 2 1 3 1 3 2 1 4 1 2 3 1 2 1 3 2 1 3 1
	// A C G A U C A G
	mfe, _ := CalculateMfe("ACGAUCAGAGAUCAGAGCAUACGACAGCAG", "..((((...))))...((........))..")
	fmt.Printf("mfe: %v", mfe)
	// Output: mfe: .......................
}
