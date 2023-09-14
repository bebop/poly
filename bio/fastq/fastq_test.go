package fastq

import (
	"io"
	"os"
	"strings"
	"testing"
)

func testException(t *testing.T, filePath string, errorString string) {
	file, err := os.Open(filePath)
	if err != nil {
		t.Errorf("Failed to read %s. Got error: %s", filePath, err)
	}
	const maxLineSize = 2 * 32 * 1024
	parser := NewParser(file, maxLineSize)
	for err == nil {
		_, err = parser.Next()
	}
	if err == nil {
		t.Errorf("%s parser should have gotten error: %s", filePath, errorString)
	}
}

func TestParseExceptions(t *testing.T) {
	testException(t, "data/nanosavseq_noseq.fastq", "no seq")
	testException(t, "data/nanosavseq_noquality.fastq", "no quality")
	testException(t, "data/nanosavseq_noidentifier.fastq", "no identifier")
	testException(t, "data/nanosavseq_emptyseq.fastq", "empty seq")
	testException(t, "data/nanosavseq_noplus.fastq", "no plus EOF")
	testException(t, "data/nanosavseq_noquality2.fastq", "no quality EOF")
}

func TestParser(t *testing.T) {
	file := strings.NewReader(`@e3cc70d5-90ef-49b6-bbe1-cfef99537d73 runid=99790f25859e24307203c25273f3a8be8283e7eb read=13956 ch=53 start_time=2020-11-11T01:49:01Z flow_cell_id=AEI083 protocol_group_id=NanoSav2 sample_id=nanosavseq2
GATGTGCGCCGTTCCAGTTGCGACGTACTATAATCCCCGGCAACACGGTGCTGATTCTCTTCCTGTTCCAGAAAGCATAAACAGATGCAAGTCTGGTGTGATTAACTTCACCAAAGGGCTGGTTGTAATATTAGGAAATCTAACAATAGATTCTGTTGGTTGGACTCTAAAATTAGAAATTTGATAGATTCCTTTTCCCAAATGAAAGTTTAACGTACACTTTGTTTCTAAAGGAAGGTCAAATTACAGTCTACAGCATCGTAATGGTTCATTTTCATTTATATTTTAATACTAGAAAAGTCCTAGGTTGAAGATAACCACATAATAAGCTGCAACTTCAGCTGTCCCAACCTGAAGAAGAATCGCAGGAGTCGAAATAACTTCTGTAAAGCAAGTAGTTTGAACCTATTGATGTTTCAACATGAGCAATACGTAACT
+
$$&%&%#$)*59;/767C378411,***,('11<;:,0039/0&()&'2(/*((4.1.09751).601+'#&&&,-**/0-+3558,/)+&)'&&%&$$'%'%'&*/5978<9;**'3*'&&A?99:;:97:278?=9B?CLJHGG=9<@AC@@=>?=>D>=3<>=>3362$%/((+/%&+//.-,%-4:+..000,&$#%$$%+*)&*0%.//*?<<;>DE>.8942&&//074&$033)*&&&%**)%)962133-%'&*99><<=1144??6.027639.011/-)($#$(/422*4;:=122>?@6964:.5'8:52)*675=:4@;323&&##'.-57*4597)+0&:7<7-550REGB21/0+*79/&/6538())+)+23665+(''$$$'-2(&&*-.-#$&%%$$,-)&$$#$'&,);;<C<@454)#'`)
	const maxLineSize = 2 * 32 * 1024
	parser := NewParser(file, maxLineSize)
	read, _ := parser.Next()
	_, err := parser.Next()
	if err != io.EOF {
		t.Errorf("Parser got unknown error: %s", err)
	}
	if read.Optionals["ch"] != "53" {
		t.Errorf("Optionals not parsed properly")
	}
}
