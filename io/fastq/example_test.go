package fastq_test

import (
	_ "embed"
	"fmt"
	"os"
	"strings"

	"github.com/TimothyStiles/poly/io/fastq"
)

//go:embed data/nanosavseq.fastq
var baseFastq string

// ExampleRead shows basic usage for Read.
func ExampleRead() {
	fastqs, _ := fastq.Read("data/nanosavseq.fastq")
	fmt.Println(fastqs[0].Identifier)
	//Output:
	//e3cc70d5-90ef-49b6-bbe1-cfef99537d73
}

// ExampleReadGz shows basic usage for ReadGz.
func ExampleReadGz() {
	fastqs, _ := fastq.ReadGz("data/nanosavseq.fastq.gz")
	fmt.Println(fastqs[0].Identifier)
	//Output:
	//e3cc70d5-90ef-49b6-bbe1-cfef99537d73
}

// ExampleWrite shows basic usage of the  writer.
func ExampleWrite() {
	fastqs, _ := fastq.Read("data/nanosavseq.fastq") // get example data
	_ = fastq.Write(fastqs, "data/test.fastq")       // write it out again
	testSequence, _ := fastq.Read("data/test.fastq") // read it in again

	os.Remove("data/test.fastq") // getting rid of test file

	fmt.Println(testSequence[0].Identifier)
	fmt.Println(testSequence[0].Sequence)
	fmt.Println(testSequence[0].Quality)
	//Output:
	//e3cc70d5-90ef-49b6-bbe1-cfef99537d73
	//GATGTGCGCCGTTCCAGTTGCGACGTACTATAATCCCCGGCAACACGGTGCTGATTCTCTTCCTGTTCCAGAAAGCATAAACAGATGCAAGTCTGGTGTGATTAACTTCACCAAAGGGCTGGTTGTAATATTAGGAAATCTAACAATAGATTCTGTTGGTTGGACTCTAAAATTAGAAATTTGATAGATTCCTTTTCCCAAATGAAAGTTTAACGTACACTTTGTTTCTAAAGGAAGGTCAAATTACAGTCTACAGCATCGTAATGGTTCATTTTCATTTATATTTTAATACTAGAAAAGTCCTAGGTTGAAGATAACCACATAATAAGCTGCAACTTCAGCTGTCCCAACCTGAAGAAGAATCGCAGGAGTCGAAATAACTTCTGTAAAGCAAGTAGTTTGAACCTATTGATGTTTCAACATGAGCAATACGTAACT
	//$$&%&%#$)*59;/767C378411,***,('11<;:,0039/0&()&'2(/*((4.1.09751).601+'#&&&,-**/0-+3558,/)+&)'&&%&$$'%'%'&*/5978<9;**'3*'&&A?99:;:97:278?=9B?CLJHGG=9<@AC@@=>?=>D>=3<>=>3362$%/((+/%&+//.-,%-4:+..000,&$#%$$%+*)&*0%.//*?<<;>DE>.8942&&//074&$033)*&&&%**)%)962133-%'&*99><<=1144??6.027639.011/-)($#$(/422*4;:=122>?@6964:.5'8:52)*675=:4@;323&&##'.-57*4597)+0&:7<7-550REGB21/0+*79/&/6538())+)+23665+(''$$$'-2(&&*-.-#$&%%$$,-)&$$#$'&,);;<C<@454)#'
}

func ExampleParser() {
	parser := fastq.NewParser(strings.NewReader(baseFastq), 2*32*1024)
	for {
		fastq, _, err := parser.ParseNext()
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
