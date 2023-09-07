package bio_test

import (
	"bytes"
	"compress/gzip"
	"fmt"
	"os"
	"strings"

	"github.com/TimothyStiles/poly/bio"
	"github.com/TimothyStiles/poly/bio/fasta"
	"github.com/TimothyStiles/poly/bio/genbank"
	"github.com/TimothyStiles/poly/bio/gff"
	"github.com/TimothyStiles/poly/bio/polyjson"
)

// This is where the integration tests that make effed up cyclic dependencies go.

func Example() {

	// Poly can take in basic gff, gbk, fasta, and JSON.
	// We call the json package "pson" (poly JSON) to prevent namespace collision with Go's standard json package.

	gffInput, _ := gff.Read("../data/ecoli-mg1655-short.gff")
	gbkInput, _ := genbank.Read("../data/puc19.gbk")
	//fastaInput, _ := fasta.Read("fasta/data/base.fasta")
	jsonInput, _ := polyjson.Read("../data/cat.json")

	// Poly can also output these file formats. Every file format has a corresponding Write function.
	_ = gff.Write(gffInput, "test.gff")
	_ = genbank.Write(gbkInput, "test.gbk")
	//_ = fasta.WriteFile(fastaInput, "test.fasta")
	_ = polyjson.Write(jsonInput, "test.json")

	// Extra tips:

	// 1. All of these file formats can be read and written in JSON format using their native schemas.
	// 2. If you want to convert from one format to another (e.g. genbank to polyjson), you can easily do so with a for-loop and some field mapping.
	// 3. Every file format is unique but they all share a common interface so you can use them with almost every native function in Poly.
}

// ExampleRead shows an example of reading a file from disk.
func ExampleRead() {
	// Read lets you read files from disk into a parser.
	file, _ := os.Open("fasta/data/base.fasta")
	parser, _ := bio.NewFastaParser(file)

	records, _ := parser.Parse()

	fmt.Println(records[1].Sequence)
	// Output: ADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNGTIDFPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLTDEEVDEMIREADIDGDGQVNYEEFVQMMTAK*
}

func ExampleReadGz() {
	// ReadGz lets you read gzipped files into a parser.
	fileGz, _ := os.Open("fasta/data/base.fasta.gz")
	file, _ := gzip.NewReader(fileGz)
	parser, _ := bio.NewFastaParser(file)
	records, _ := parser.Parse()

	fmt.Println(records[1].Sequence)
	// Output: ADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNGTIDFPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLTDEEVDEMIREADIDGDGQVNYEEFVQMMTAK*
}

func ExampleNewParserGz() {
	// First, lets make a file that is gzip'd, represented by this
	// buffer.
	var file bytes.Buffer
	zw := gzip.NewWriter(&file)
	_, _ = zw.Write([]byte(`>gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]
LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLV
EWIWGGFSVDKATLNRFFAFHFILPFTMVALAGVHLTFLHETGSNNPLGLTSDSDKIPFHPYYTIKDFLG
LLILILLLLLLALLSPDMLGDPDNHMPADPLNTPLHIKPEWYFLFAYAILRSVPNKLGGVLALFLSIVIL
GLMPFLHTSKHRSMMLRPLSQALFWTLTMDLLTLTWIGSQPVEYPYTIIGQMASILYFSIILAFLPIAGX
IENY

>MCHU - Calmodulin - Human, rabbit, bovine, rat, and chicken
ADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNGTID
FPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLTDEEVDEMIREA
DIDGDGQVNYEEFVQMMTAK*`))
	zw.Close()

	fileDecompressed, _ := gzip.NewReader(&file) // Decompress the file
	parser, _ := bio.NewFastaParser(fileDecompressed)
	records, _ := parser.Parse() // Parse all data records from file

	fmt.Println(records[1].Sequence)
	// Output: ADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNGTIDFPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLTDEEVDEMIREADIDGDGQVNYEEFVQMMTAK*
}

func ExampleParseWithHeader() {
	// The following can be replaced with a any io.Reader. For example,
	// `file, err := os.Open(path)` for file would also work.
	file := strings.NewReader(`#slow5_version	0.2.0
#num_read_groups	1
@asic_id	4175987214
#char*	uint32_t	double	double	double	double	uint64_t	int16_t*	uint64_t	int32_t	uint8_t	double	enum{unknown,partial,mux_change,unblock_mux_change,data_service_unblock_mux_change,signal_positive,signal_negative}	char*
#read_id	read_group	digitisation	offset	range	sampling_rate	len_raw_signal	raw_signal	start_time	read_number	start_mux	median_before	end_reason	channel_number
0026631e-33a3-49ab-aa22-3ab157d71f8b	0	8192	16	1489.52832	4000	5347	430,472,463	8318394	5383	1	219.133423	5	10
`)
	parser, _ := bio.NewSlow5Parser(file)
	reads, header, _ := parser.ParseWithHeader() // Parse all data records from file

	fmt.Printf("%s, %s\n", header.HeaderValues[0].Slow5Version, reads[0].ReadID)
	// Output: 0.2.0, 0026631e-33a3-49ab-aa22-3ab157d71f8b
}

func ExampleParseToChannel() {
	// The following can be replaced with a any io.Reader. For example,
	// `file, err := os.Open(path)` for file would also work.
	file := strings.NewReader(`>gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]
LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLV
EWIWGGFSVDKATLNRFFAFHFILPFTMVALAGVHLTFLHETGSNNPLGLTSDSDKIPFHPYYTIKDFLG
LLILILLLLLLALLSPDMLGDPDNHMPADPLNTPLHIKPEWYFLFAYAILRSVPNKLGGVLALFLSIVIL
GLMPFLHTSKHRSMMLRPLSQALFWTLTMDLLTLTWIGSQPVEYPYTIIGQMASILYFSIILAFLPIAGX
IENY

>MCHU - Calmodulin - Human, rabbit, bovine, rat, and chicken
ADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNGTID
FPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLTDEEVDEMIREA
DIDGDGQVNYEEFVQMMTAK*`)
	parser, _ := bio.NewFastaParser(file)

	channel := make(chan fasta.Record)
	go parser.ParseToChannel(channel)

	var records []fasta.Record
	for record := range channel {
		records = append(records, record)
	}

	fmt.Println(records[1].Sequence)
	// Output: ADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNGTIDFPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLTDEEVDEMIREADIDGDGQVNYEEFVQMMTAK*
}

func ExampleWriteAll() {
	// The following can be replaced with a any io.Reader. For example,
	// `file, err := os.Open(path)` for file would also work.
	file := strings.NewReader(`#slow5_version	0.2.0
#num_read_groups	1
@asic_id	4175987214
#char*	uint32_t	double	double	double	double	uint64_t	int16_t*	uint64_t	int32_t	uint8_t	double	enum{unknown,partial,mux_change,unblock_mux_change,data_service_unblock_mux_change,signal_positive,signal_negative}	char*
#read_id	read_group	digitisation	offset	range	sampling_rate	len_raw_signal	raw_signal	start_time	read_number	start_mux	median_before	end_reason	channel_number
0026631e-33a3-49ab-aa22-3ab157d71f8b	0	8192	16	1489.52832	4000	5347	430,472,463	8318394	5383	1	219.133423	5	10
`)
	parser, _ := bio.NewSlow5Parser(file)
	reads, header, _ := parser.ParseWithHeader() // Parse all data records from file

	// Write the files to an io.Writer.
	// All headers and all records implement io.WriterTo interfaces.
	var buffer bytes.Buffer
	_, _ = header.WriteTo(&buffer)
	for _, read := range reads {
		read.WriteTo(&buffer)
	}

	fmt.Println(string(buffer.Bytes()))
	// Output: #slow5_version	0.2.0
	//#num_read_groups	1
	//@asic_id	4175987214
	//#char*	uint32_t	double	double	double	double	uint64_t	int16_t*	uint64_t	int32_t	uint8_t	double	enum{unknown,partial,mux_change,unblock_mux_change,data_service_unblock_mux_change,signal_positive,signal_negative}	char*
	//#read_id	read_group	digitisation	offset	range	sampling_rate	len_raw_signal	raw_signal	start_time	read_number	start_mux	median_before	end_reason	channel_number
	//0026631e-33a3-49ab-aa22-3ab157d71f8b	0	8192	16	1489.52832	4000	5347	430,472,463	8318394	5383	1	219.133423	5	10
	//
	//
}

func ExampleFasta() {
	// The following can be replaced with a any io.Reader. For example,
	// `file, err := os.Open(path)` for file would also work.
	file := strings.NewReader(`>gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]
LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLV
EWIWGGFSVDKATLNRFFAFHFILPFTMVALAGVHLTFLHETGSNNPLGLTSDSDKIPFHPYYTIKDFLG
LLILILLLLLLALLSPDMLGDPDNHMPADPLNTPLHIKPEWYFLFAYAILRSVPNKLGGVLALFLSIVIL
GLMPFLHTSKHRSMMLRPLSQALFWTLTMDLLTLTWIGSQPVEYPYTIIGQMASILYFSIILAFLPIAGX
IENY

>MCHU - Calmodulin - Human, rabbit, bovine, rat, and chicken
ADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNGTID
FPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLTDEEVDEMIREA
DIDGDGQVNYEEFVQMMTAK*`)
	parser, _ := bio.NewFastaParser(file)
	records, _ := parser.Parse() // Parse all data records from file

	fmt.Println(records[1].Sequence)
	// Output: ADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNGTIDFPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLTDEEVDEMIREADIDGDGQVNYEEFVQMMTAK*
}

func ExampleSlow5() {
	// The following can be replaced with a any io.Reader. For example,
	// `file, err := os.Open(path)` for file would also work.
	file := strings.NewReader(`#slow5_version	0.2.0
#num_read_groups	1
@asic_id	4175987214
#char*	uint32_t	double	double	double	double	uint64_t	int16_t*	uint64_t	int32_t	uint8_t	double	enum{unknown,partial,mux_change,unblock_mux_change,data_service_unblock_mux_change,signal_positive,signal_negative}	char*
#read_id	read_group	digitisation	offset	range	sampling_rate	len_raw_signal	raw_signal	start_time	read_number	start_mux	median_before	end_reason	channel_number
0026631e-33a3-49ab-aa22-3ab157d71f8b	0	8192	16	1489.52832	4000	5347	430,472,463	8318394	5383	1	219.133423	5	10
`)
	parser, _ := bio.NewSlow5Parser(file)
	reads, _ := parser.Parse() // Parse all data records from file

	fmt.Println(reads[0].RawSignal)
	// Output: [430 472 463]
}
