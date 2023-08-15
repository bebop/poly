package slow5_test

import (
	"fmt"
	"os"

	"github.com/TimothyStiles/poly/io/slow5"
)

func ExampleNewParser() {
	// example.slow5 is a file I generated using slow5tools from nanopore fast5
	// run where I was testing using nanopore for doing COVID testing. It
	// contains real nanopore data.
	file, _ := os.Open("data/example.slow5")
	defer file.Close()
	// Set maxLineSize to 64kb. If you expect longer reads,
	// make maxLineSize longer!
	const maxLineSize = 2 * 32 * 1024
	parser, _, _ := slow5.NewParser(file, maxLineSize)

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
