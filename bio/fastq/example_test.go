package fastq_test

import (
	_ "embed"
	"fmt"
	"strings"

	"github.com/TimothyStiles/poly/bio/fastq"
)

//go:embed data/nanosavseq.fastq
var baseFastq string

func ExampleParser() {
	parser := fastq.NewParser(strings.NewReader(baseFastq), 2*32*1024)
	for {
		fastq, err := parser.Next()
		if err != nil {
			fmt.Println(err)
			break
		}
		fmt.Println(fastq.Identifier)
	}
	//Output:
	//e3cc70d5-90ef-49b6-bbe1-cfef99537d73
	//92728f25-b658-426c-8cd7-d82dc70dbf71
	//60907b6b-5e38-498e-9c07-f036ebd8c658
	//990e110e-5e50-41a2-8ad5-92044d4465b8
	//EOF
}
