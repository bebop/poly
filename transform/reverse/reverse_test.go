package reverse

import "fmt"

func ExampleComplement() {
	s := "ACGT"
	c := Complement(s)
	fmt.Println(c)

	// Output:
	// TACG
}
