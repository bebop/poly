package slow5_test

import (
	"fmt"
	"os"

	"github.com/TimothyStiles/poly/io/slow5"
)

func ExampleNewParser() {
	file, _ := os.Open("data/example.slow5")
	parser, _, _ := slow5.NewParser(file)

	var outputReads []slow5.Read
	for {
		read, err := parser.ParseNext()
		if err != nil {
			// Break at EOF
			break
		}
		outputReads = append(outputReads, read)
	}

	fmt.Println(outputReads[0].RawSignal[0:10])
	// Output: [430 472 463 467 454 465 463 450 450 449]
}
