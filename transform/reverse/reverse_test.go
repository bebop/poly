package reverse

import "fmt"

func ExampleComplement() {
	sequence := "GATTACA"
	reverseComplement := Complement(sequence)
	fmt.Println(reverseComplement)

	// Output: TGTAATC
}
