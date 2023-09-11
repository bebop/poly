package slow5_test

import (
	"fmt"
	"os"

	"github.com/TimothyStiles/poly/bio/slow5"
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
	parser, _ := slow5.NewParser(file, maxLineSize)

	var outputReads []slow5.Read
	for {
		read, err := parser.Next()
		if err != nil {
			// Break at EOF
			break
		}
		outputReads = append(outputReads, *read)
	}

	fmt.Println(outputReads[0].RawSignal[0:10])
	// Output: [430 472 463 467 454 465 463 450 450 449]
}

func ExampleSvbCompressRawSignal() {
	// example.slow5 is a file I generated using slow5tools from nanopore fast5
	// run where I was testing using nanopore for doing COVID testing. It
	// contains real nanopore data.
	file, _ := os.Open("data/example.slow5")
	defer file.Close()
	// Set maxLineSize to 64kb. If you expect longer reads,
	// make maxLineSize longer!
	const maxLineSize = 2 * 32 * 1024
	parser, _ := slow5.NewParser(file, maxLineSize)
	read, _ := parser.Next()

	// Get the raw signal from a read
	rawSignal := read.RawSignal

	// Compress that raw signal into a mask and data
	mask, data := slow5.SvbCompressRawSignal(rawSignal)

	// Decompress mask and data back into raw signal
	rawSignalDecompressed := slow5.SvbDecompressRawSignal(len(rawSignal), mask, data)

	for idx := range rawSignal {
		if rawSignal[idx] != rawSignalDecompressed[idx] {
			fmt.Println("Compression failed!")
		}
	}
	fmt.Println(data[:10])

	// Output: [174 1 216 1 207 1 211 1 198 1]
}
