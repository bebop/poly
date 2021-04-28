package main

import (
	"bytes"
	"io/ioutil"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"testing"

	"github.com/TimothyStiles/poly"
	"github.com/google/go-cmp/cmp"
	"github.com/google/go-cmp/cmp/cmpopts"
)

/******************************************************************************
Oct, 15, 2020

Testing command line utilities can be annoying.

The way poly does it is by using the cli.app  object to spoof input and output
via cli.app.reader and cli.app.writer. This is the only way to get true stack
traceable coverage.
******************************************************************************/

var testFilePaths = []string{"../data/puc19.gbk", "../data/ecoli-mg1655-short.gff", "../data/sample.json", "../data/base.fasta"}

func TestMain(t *testing.T) {
	rescueStdout := os.Stdout
	_, w, _ := os.Pipe()
	os.Stdout = w

	arg := os.Args[0:1]
	os.Args = append(arg, "-h")
	main()
	os.Args = os.Args[0:1]
	w.Close()
	os.Stdout = rescueStdout
}
func TestConvertPipe(t *testing.T) {

	for _, match := range testFilePaths {
		extension := filepath.Ext(match)[1:]
		var writeBuffer bytes.Buffer
		app := application()
		app.Writer = &writeBuffer

		args := os.Args[0:1]                                    // Name of the program.
		args = append(args, "c", "-i", extension, "-o", "json") // Append a flag

		file, _ := ioutil.ReadFile(match)
		app.Reader = bytes.NewReader(file)

		err := app.Run(args)

		if err != nil {
			t.Fatalf("Run error: %s", err)
		}

		// getting test sequence from non-pipe io to compare against io to stdout
		baseTestSequence := parseExt(match)
		pipeOutputTestSequence := poly.ParseJSON(writeBuffer.Bytes())

		if diff := cmp.Diff(baseTestSequence, pipeOutputTestSequence, cmpopts.IgnoreFields(poly.Feature{}, "ParentSequence")); diff != "" {
			t.Errorf(" mismatch converting from %q to json (-want +got):\n%s", extension, diff)
		}
	}

}

func TestConvertIO(t *testing.T) {
	for _, match := range testFilePaths {
		extension := filepath.Ext(match)[1:]
		var writeBuffer bytes.Buffer
		app := application()
		app.Writer = &writeBuffer

		args := os.Args[0:1] // Name of the program.
		args = append(args, "c", "-i", extension, "-o", extension)
		file, _ := ioutil.ReadFile(match)
		app.Reader = bytes.NewReader(file)

		err := app.Run(args)

		if err != nil {
			t.Fatalf("Run error: %s", err)
		}

		// getting test sequence from non-pipe io to compare against io to stdout
		baseTestSequence := parseExt(match)

		pipeOutputTestSequence := parseFlag(writeBuffer.Bytes(), extension)

		if diff := cmp.Diff(baseTestSequence, pipeOutputTestSequence, cmpopts.IgnoreFields(poly.Feature{}, "ParentSequence")); diff != "" {
			t.Errorf(" mismatch reading and writing %q (-want +got):\n%s", extension, diff)
		}
	}
}

func TestConvertWriteFile(t *testing.T) {

	for _, match := range testFilePaths {
		extension := filepath.Ext(match)
		testOutputPath := "../data/test" + extension

		loopApp := application()

		args := os.Args[0:1] // Name of the program.
		args = append(args, "c", "-o", testOutputPath, match)

		err := loopApp.Run(args)

		if err != nil {
			t.Fatalf("Run error: %s", err)
		}

		// getting test sequence from non-pipe io to compare against io to stdout
		baseTestSequence := parseExt(match)

		outputSequence := parseExt(testOutputPath)
		os.Remove(testOutputPath)

		if diff := cmp.Diff(baseTestSequence, outputSequence, cmpopts.IgnoreFields(poly.Feature{}, "ParentSequence")); diff != "" {
			t.Errorf(" mismatch reading and writing %q (-want +got):\n%s", extension, diff)
		}
	}
}

// The test below passes on macosx but not ubuntu. Could someone try debugging it on a device running ubuntu?
func TestConvertWarning(t *testing.T) {

	testFilePath := "../data/puc19.gbk"
	var writeBuffer bytes.Buffer

	app := application()
	app.Writer = &writeBuffer

	// testing file matching hash
	args := os.Args[0:1]                                       // Name of the program.
	args = append(args, "c", "-o", testFilePath, testFilePath) // Append a flag

	err := app.Run(args)

	if err != nil {
		t.Fatalf("Run error: %s", err)
	}

	if writeBuffer.Len() == 0 {
		t.Error("TestConvertWarning did not write output to desired writer.")
	}

	warningOutputString := writeBuffer.String()

	warningTestString := "WARNING: " + "input path and output path match. File: " + testFilePath + ". Skipping to prevent possible data loss. Try providing a full path with a different name using the -o flag.\n"
	if warningOutputString != warningTestString {
		t.Errorf("TestConvertWarning has failed. Returned %q, want %q", warningOutputString, warningTestString)
	}
}

func TestConvertPipeString(t *testing.T) {

	var writeBuffer bytes.Buffer
	app := application()
	app.Writer = &writeBuffer

	args := os.Args[0:1] // Name of the program.
	args = append(args, "c", "-i", "string", "-o", "string")
	file := []byte("ACTGATCGATCAGCTACGACTAGCATCAGCATCAGCATCAGCTACGACTAG")
	app.Reader = bytes.NewReader(file)

	err := app.Run(args)

	if err != nil {
		t.Fatalf("Run error: %s", err)
	}

	// getting test sequence from non-pipe io to compare against io to stdout
	pipeOutputTestSequence := parseFlag(writeBuffer.Bytes(), "string")

	if string(file) != pipeOutputTestSequence.Sequence {
		t.Errorf(" mismatch reading and writing strings from stdin to stdout.")
	}
}

func TestConvertFile(t *testing.T) {

	app := application()

	args := os.Args[0:1]                                                                // Name of the program.
	args = append(args, "c", "-o", "json", "../data/puc19.gbk", "../data/t4_intron.gb") // Append a flag

	err := app.Run(args)

	if err != nil {
		t.Fatalf("Run error: %s", err)
	}

	puc19InputTestSequence := poly.ReadGbk("../data/puc19.gbk")
	puc19OutputTestSequence := poly.ReadJSON("../data/puc19.json")

	//clearing test data.
	os.Remove("../data/puc19.json")

	// compared input gff from resulting output json. Fail test and print diff if error.
	if diff := cmp.Diff(puc19InputTestSequence, puc19OutputTestSequence, cmpopts.IgnoreFields(poly.Feature{}, "ParentSequence")); diff != "" {
		t.Errorf(" mismatch from concurrent gbk input test (-want +got):\n%s", diff)
	}

	t4InputTestSequence := poly.ReadGbk("../data/t4_intron.gb")
	t4OutputTestSequence := poly.ReadJSON("../data/t4_intron.json")

	// clearing test data.
	os.Remove("../data/t4_intron.json")

	// compared input gbk from resulting output json. Fail test and print diff if error.
	if diff := cmp.Diff(t4InputTestSequence, t4OutputTestSequence, cmpopts.IgnoreFields(poly.Feature{}, "ParentSequence")); diff != "" {
		t.Errorf(" mismatch from concurrent gbk input test (-want +got):\n%s", diff)
	}
}

func TestHashFile(t *testing.T) {

	testFilePath := "../data/puc19.gbk"
	puc19GbkBlake3Hash := "4b0616d1b3fc632e42d78521deb38b44fba95cca9fde159e01cd567fa996ceb9"
	var writeBuffer bytes.Buffer

	app := application()
	app.Writer = &writeBuffer

	// testing file matching hash
	args := os.Args[0:1]                      // Name of the program.
	args = append(args, "hash", testFilePath) // Append a flag

	err := app.Run(args)

	if err != nil {
		t.Fatalf("Run error: %s", err)
	}

	if writeBuffer.Len() == 0 {
		t.Error("TestHash did not write output to desired writer.")
	}

	hashOutputString := writeBuffer.String()

	// use same function that outputs hash results to format test case.
	testHashString := formatHashOutput(puc19GbkBlake3Hash, testFilePath)
	if hashOutputString != testHashString {
		t.Errorf("TestHashFile has failed. Returned %q, want %q", hashOutputString, puc19GbkBlake3Hash)
	}

}

func TestHashPipe(t *testing.T) {

	puc19GbkBlake3Hash := "4b0616d1b3fc632e42d78521deb38b44fba95cca9fde159e01cd567fa996ceb9"
	testFilePath := "../data/puc19.gbk"

	var writeBuffer bytes.Buffer

	// create a mock application
	app := application()
	app.Writer = &writeBuffer
	file, _ := ioutil.ReadFile(testFilePath)
	app.Reader = bytes.NewReader(file)

	args := os.Args[0:1]                     // Name of the program.
	args = append(args, "hash", "-i", "gbk") // Append a flag

	err := app.Run(args)

	if err != nil {
		t.Fatalf("Run error: %s", err)
	}

	hashOutputString := writeBuffer.String()

	// use same function that outputs hash results to format test case.
	testHashString := formatHashOutput(puc19GbkBlake3Hash, "-")

	if hashOutputString != testHashString {
		t.Errorf("TestHashPipe has failed. Returned %q, want %q", hashOutputString, puc19GbkBlake3Hash)
	}

}

func TestHashFunctions(t *testing.T) {

	testFilePath := "../data/puc19.gbk"

	hashfunctionsTestMap := map[string]string{
		"":            "4b0616d1b3fc632e42d78521deb38b44fba95cca9fde159e01cd567fa996ceb9", // default / defaults to blake3
		"BLAKE2b_256": "51f33a2aaa3ad3884e18822b02ffd37ee27888c7d6103c97bf9628843e1fc445",
		"BLAKE2b_384": "c56ff217ca0e09b1d20a910ec51cb17fd828ae47cf6208bfca2381021cdff86ffd6d0975f6955593400232b1cafd96fa",
		"BLAKE2b_512": "9bc9066ea9f7bb58cc432175a9a99ebb2332247fc97b4f73318882e83bcff647d88503c6bfb02f6d6468af9b88dc303f3270b50922585e669e0c990fb384f4f6",
		"BLAKE2s_256": "3ca14269a42d9ea4c619ced109f879f1bba76469be519b9644ed977753e05353",
		"BLAKE3":      "4b0616d1b3fc632e42d78521deb38b44fba95cca9fde159e01cd567fa996ceb9",
		"MD5":         "2c6d6bdc39c1d59c9aa3d332ddc28190",
		"NO":          "aaaaaaaccaccgctaccagcggtggtttgtttgccggatcaagagctaccaactctttttccgaaggtaactggcttcagcagagcgcagataccaaatactgttcttctagtgtagccgtagttaggccaccacttcaagaactctgtagcaccgcctacatacctcgctctgctaatcctgttaccagtggctgctgccagtggcgataagtcgtgtcttaccgggttggactcaagacgatagttaccggataaggcgcagcggtcgggctgaacggggggttcgtgcacacagcccagcttggagcgaacgacctacaccgaactgagatacctacagcgtgagctatgagaaagcgccacgcttcccgaagggagaaaggcggacaggtatccggtaagcggcagggtcggaacaggagagcgcacgagggagcttccagggggaaacgcctggtatctttatagtcctgtcgggtttcgccacctctgacttgagcgtcgatttttgtgatgctcgtcaggggggcggagcctatggaaaaacgccagcaacgcggcctttttacggttcctggccttttgctggccttttgctcacatgttctttcctgcgttatcccctgattctgtggataaccgtattaccgcctttgagtgagctgataccgctcgccgcagccgaacgaccgagcgcagcgagtcagtgagcgaggaagcggaagagcgcccaatacgcaaaccgcctctccccgcgcgttggccgattcattaatgcagctggcacgacaggtttcccgactggaaagcgggcagtgagcgcaacgcaattaatgtgagttagctcactcattaggcaccccaggctttacactttatgcttccggctcgtatgttgtgtggaattgtgagcggataacaatttcacacaggaaacagctatgaccatgattacgccaagcttgcatgcctgcaggtcgactctagaggatccccgggtaccgagctcgaattcactggccgtcgttttacaacgtcgtgactgggaaaaccctggcgttacccaacttaatcgccttgcagcacatccccctttcgccagctggcgtaatagcgaagaggcccgcaccgatcgcccttcccaacagttgcgcagcctgaatggcgaatggcgcctgatgcggtattttctccttacgcatctgtgcggtatttcacaccgcatatggtgcactctcagtacaatctgctctgatgccgcatagttaagccagccccgacacccgccaacacccgctgacgcgccctgacgggcttgtctgctcccggcatccgcttacagacaagctgtgaccgtctccgggagctgcatgtgtcagaggttttcaccgtcatcaccgaaacgcgcgagacgaaagggcctcgtgatacgcctatttttataggttaatgtcatgataataatggtttcttagacgtcaggtggcacttttcggggaaatgtgcgcggaacccctatttgtttatttttctaaatacattcaaatatgtatccgctcatgagacaataaccctgataaatgcttcaataatattgaaaaaggaagagtatgagtattcaacatttccgtgtcgcccttattcccttttttgcggcattttgccttcctgtttttgctcacccagaaacgctggtgaaagtaaaagatgctgaagatcagttgggtgcacgagtgggttacatcgaactggatctcaacagcggtaagatccttgagagttttcgccccgaagaacgttttccaatgatgagcacttttaaagttctgctatgtggcgcggtattatcccgtattgacgccgggcaagagcaactcggtcgccgcatacactattctcagaatgacttggttgagtactcaccagtcacagaaaagcatcttacggatggcatgacagtaagagaattatgcagtgctgccataaccatgagtgataacactgcggccaacttacttctgacaacgatcggaggaccgaaggagctaaccgcttttttgcacaacatgggggatcatgtaactcgccttgatcgttgggaaccggagctgaatgaagccataccaaacgacgagcgtgacaccacgatgcctgtagcaatggcaacaacgttgcgcaaactattaactggcgaactacttactctagcttcccggcaacaattaatagactggatggaggcggataaagttgcaggaccacttctgcgctcggcccttccggctggctggtttattgctgataaatctggagccggtgagcgtgggtctcgcggtatcattgcagcactggggccagatggtaagccctcccgtatcgtagttatctacacgacggggagtcaggcaactatggatgaacgaaatagacagatcgctgagataggtgcctcactgattaagcattggtaactgtcagaccaagtttactcatatatactttagattgatttaaaacttcatttttaatttaaaaggatctaggtgaagatcctttttgataatctcatgaccaaaatcccttaacgtgagttttcgttccactgagcgtcagaccccgtagaaaagatcaaaggatcttcttgagatcctttttttctgcgcgtaatctgctgcttgcaaac",
		"SHA1":        "9fe3597d67e31e003bc0aa04054eb37e85e26ffd",
		"SHA244":      "9b674798798e8cfe1b341fe77f575357bb5e4038721ae9fbe41bdd3b",
		"SHA256":      "4a9da025ff81f65a25a73137e5e528ee8923dc47b29a6d12a64989d67aec5f2a",
		"SHA384":      "06f633efe87daf8f903bd1421670f8b4367394ee496408c31a45ae4866a82ec77237ad8f972717e5c0df8bc8bc2db181",
		"SHA3_224":    "2624ef4b8a88e9ad0504a33b6bb236ae9409e83a34e07c76fefbf20b",
		"SHA3_256":    "7ff3f5c0e95705e82701c62ab68e1a6af4591f26b413990eed37d661f4b3ffcb",
		"SHA3_384":    "3f631cbf030e8cbd5c72ec4ba2e775bae4f61c2f5c30c873178f62871d4da413266f908c4aa130ece5412b0243e15101",
		"SHA3_512":    "2f831d14e072044715c4ebfba80920aeefbf86728cf1132afe2688e7fe51a8d7e14203e95da88fd0bcefe064a588cda2497e5be69f304d35088708a0e5f7e843",
		"SHA512":      "66ee5400a01c238d6be1a04e1133f5f2a9fd178ff9c57c5e695b47191aadaddbff5eb885170b3e1715db9c06b1fd4ee88cd599ba98857978b2648e76c7c80ad9",
		"SHA512_224":  "1e543b46774dae76f8a8e3014a214b70b85cb49e7b3dba806b36263a",
		"SHA512_256":  "7691fa37e37ee5a56d0c844d8fde0682713ab4732eeedb973d1a485f29348ac6",
	}

	// sorting keys for easier to read debugging output.
	keys := make([]string, 0, len(hashfunctionsTestMap))
	for k := range hashfunctionsTestMap {
		keys = append(keys, k)
	}

	sort.Strings(keys)

	for _, key := range keys {

		var writeBuffer bytes.Buffer

		app := application()
		app.Writer = &writeBuffer
		args := os.Args[0:1]
		args = append(args, "hash", "-f", key, testFilePath)
		err := app.Run(args)

		if err != nil {
			t.Fatalf("Run error: %s", err)
		}

		if err != nil {
			t.Fatalf("Run error: %s", err)
		}

		if writeBuffer.Len() == 0 {
			t.Error("TestHash did not write output to desired writer.")
		}

		hashOutputString := writeBuffer.String()
		hashTestString := formatHashOutput(hashfunctionsTestMap[key], testFilePath)
		// fmt.Println("\""+key+"\"", ": ", "\""+hashOutputString+"\",") // <- prints out every
		if hashOutputString != hashTestString {
			t.Errorf("TestHashFunctions for function %q has failed. Returned %q, want %q", key, hashOutputString, hashTestString)
		}
	}
}

func TestOptimizeString(t *testing.T) {

	gfpTranslation := "MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK*"
	var writeBuffer bytes.Buffer

	app := application()
	app.Writer = &writeBuffer
	app.Reader = bytes.NewBufferString(gfpTranslation)

	args := os.Args[0:1]                                      // Name of the program.
	args = append(args, "op", "-aa", "-wt", "data/puc19.gbk") // Append a flag
	err := app.Run(args)

	if err != nil {
		t.Fatalf("Run error: %s", err)
	}
	app.Reader = os.Stdin

	// should return codon optimized sequence
	optimizeOutputString := strings.TrimSpace(writeBuffer.String())

	translation := poly.Translate(optimizeOutputString, poly.GetCodonTable(1))

	if translation != gfpTranslation {
		t.Errorf("TestOptimizeCommand for string output has failed. Returned %q, want %q", translation, gfpTranslation)
	}

}

func TestTranslationString(t *testing.T) {
	gfpDnaSequence := "ATGGCTAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCTACATACGGAAAGCTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCATATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAACGCACTATATCTTTCAAAGATGACGGGAACTACAAGACGCGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATCGTATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTCGGACACAAACTCGAGTACAACTATAACTCACACAATGTATACATCACGGCAGACAAACAAAAGAATGGAATCAAAGCTAACTTCAAAATTCGCCACAACATTGAAGATGGATCCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCGACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGCGTGACCACATGGTCCTTCTTGAGTTTGTAACTGCTGCTGGGATTACACATGGCATGGATGAGCTCTACAAATAA"
	gfpTranslation := "MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK*"

	var writeBuffer bytes.Buffer

	app := application()
	app.Writer = &writeBuffer
	app.Reader = bytes.NewBufferString(gfpDnaSequence)

	args := os.Args[0:1]                   // Name of the program.
	args = append(args, "tr", "-ct", "11") // Append a flag
	err := app.Run(args)

	if err != nil {
		t.Fatalf("Run error: %s", err)
	}
	app.Reader = os.Stdin

	translation := strings.TrimSpace(writeBuffer.String())

	if translation != gfpTranslation {
		t.Errorf("TestTranslationString for string output has failed. Returned %q, want %q", translation, gfpTranslation)
	}

}
