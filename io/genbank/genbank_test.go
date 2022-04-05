package genbank

import (
	"errors"
	"fmt"
	"io/ioutil"
	"os"
	"path/filepath"
	"regexp"
	"strings"
	"testing"

	"github.com/google/go-cmp/cmp"
	"github.com/google/go-cmp/cmp/cmpopts"
	"github.com/stretchr/testify/assert"
)

/******************************************************************************

Gbk/gb/genbank related benchmarks begin here.

******************************************************************************/

func ExampleRead() {
	sequence, _ := Read("../../data/puc19.gbk")
	fmt.Println(sequence.Meta.Locus.ModificationDate)
	// Output: 22-OCT-2019
}

func ExampleParse() {
	file, _ := ioutil.ReadFile("../../data/puc19.gbk")
	sequence, _ := Parse(file)

	fmt.Println(sequence.Meta.Locus.ModificationDate)
	// Output: 22-OCT-2019
}

func ExampleBuild() {
	sequence, _ := Read("../../data/puc19.gbk")
	gbkBytes := Build(sequence)
	testSequence, _ := Parse(gbkBytes)

	fmt.Println(testSequence.Meta.Locus.ModificationDate)
	// Output: 22-OCT-2019
}

func ExampleWrite() {
	tmpDataDir, err := ioutil.TempDir("", "data-*")
	if err != nil {
		fmt.Println(err.Error())
	}
	defer os.RemoveAll(tmpDataDir)

	sequence, _ := Read("../../data/puc19.gbk")

	tmpGbkFilePath := filepath.Join(tmpDataDir, "puc19.gbk")
	_ = Write(sequence, tmpGbkFilePath)

	testSequence, _ := Read(tmpGbkFilePath)

	fmt.Println(testSequence.Meta.Locus.ModificationDate)
	// Output: 22-OCT-2019
}

func TestGbkIO(t *testing.T) {
	tmpDataDir, err := ioutil.TempDir("", "data-*")
	if err != nil {
		t.Error(err)
	}
	defer os.RemoveAll(tmpDataDir)

	gbk, _ := Read("../../data/puc19.gbk")

	tmpGbkFilePath := filepath.Join(tmpDataDir, "puc19.gbk")
	_ = Write(gbk, tmpGbkFilePath)

	writeTestGbk, _ := Read(tmpGbkFilePath)
	if diff := cmp.Diff(gbk, writeTestGbk, []cmp.Option{cmpopts.IgnoreFields(Feature{}, "ParentSequence"), cmpopts.IgnoreFields(Meta{}, "CheckSum")}...); diff != "" {
		t.Errorf("Parsing the output of Build() does not produce the same output as parsing the original file read with Read(). Got this diff:\n%s", diff)
	}

	// Test multiline Genbank features
	pichia, _ := Read("../../data/pichia_chr1_head.gb")
	var multilineOutput string
	for _, feature := range pichia.Features {
		multilineOutput = feature.Location.GbkLocationString
	}
	if multilineOutput != "join(<459260..459456,459556..459637,459685..459739,459810..>460126)" {
		t.Errorf("Failed to parse multiline genbank feature string")
	}
}

func TestGbkLocationStringBuilder(t *testing.T) {
	tmpDataDir, err := ioutil.TempDir("", "data-*")
	if err != nil {
		t.Error(err)
	}
	defer os.RemoveAll(tmpDataDir)

	scrubbedGbk, _ := Read("../../data/sample.gbk")

	// removing gbkLocationString from features to allow testing for gbkLocationBuilder
	for featureIndex := range scrubbedGbk.Features {
		scrubbedGbk.Features[featureIndex].Location.GbkLocationString = ""
	}

	tmpGbkFilePath := filepath.Join(tmpDataDir, "sample.gbk")
	_ = Write(scrubbedGbk, tmpGbkFilePath)

	testInputGbk, _ := Read("../../data/sample.gbk")
	testOutputGbk, _ := Read(tmpGbkFilePath)

	if diff := cmp.Diff(testInputGbk, testOutputGbk, []cmp.Option{cmpopts.IgnoreFields(Feature{}, "ParentSequence"), cmpopts.IgnoreFields(Meta{}, "CheckSum")}...); diff != "" {
		t.Errorf("Issue with partial location building. Parsing the output of Build() does not produce the same output as parsing the original file read with Read(). Got this diff:\n%s", diff)
	}
}

func TestGbLocationStringBuilder(t *testing.T) {
	tmpDataDir, err := ioutil.TempDir("", "data-*")
	if err != nil {
		t.Error(err)
	}
	defer os.RemoveAll(tmpDataDir)

	scrubbedGb, _ := Read("../../data/t4_intron.gb")

	// removing gbkLocationString from features to allow testing for gbkLocationBuilder
	for featureIndex := range scrubbedGb.Features {
		scrubbedGb.Features[featureIndex].Location.GbkLocationString = ""
	}

	tmpGbFilePath := filepath.Join(tmpDataDir, "t4_intron_test.gb")
	_ = Write(scrubbedGb, tmpGbFilePath)

	testInputGb, _ := Read("../../data/t4_intron.gb")
	testOutputGb, _ := Read(tmpGbFilePath)

	if diff := cmp.Diff(testInputGb, testOutputGb, []cmp.Option{cmpopts.IgnoreFields(Feature{}, "ParentSequence"), cmpopts.IgnoreFields(Meta{}, "CheckSum")}...); diff != "" {
		t.Errorf("Issue with either Join or complement location building. Parsing the output of Build() does not produce the same output as parsing the original file read with Read(). Got this diff:\n%s", diff)
	}
}

func TestPartialLocationParseRegression(t *testing.T) {
	gbk, _ := Read("../../data/sample.gbk")

	for _, feature := range gbk.Features {
		if feature.Location.GbkLocationString == "687..3158>" && (feature.Location.Start != 686 || feature.Location.End != 3158) {
			t.Errorf("Partial location for three prime location parsing has failed. Parsing the output of Build() does not produce the same output as parsing the original file read with Read()")
		} else if feature.Location.GbkLocationString == "<1..206" && (feature.Location.Start != 0 || feature.Location.End != 206) {
			t.Errorf("Partial location for five prime location parsing has failed. Parsing the output of Build() does not produce the same output as parsing the original file read with Read().")
		}
	}
}

func TestSnapgeneGenbankRegression(t *testing.T) {
	snapgene, _ := Read("../../data/puc19_snapgene.gb")

	if snapgene.Sequence == "" {
		t.Errorf("Parsing snapgene returned an empty string")
	}
}

func TestGetSequenceMethod(t *testing.T) {

	gbk, _ := Read("../../data/t4_intron.gb")

	// Check to see if GetSequence method works on Features struct
	feature := gbk.Features[1].GetSequence()
	seq := "atgagattacaacgccagagcatcaaagattcagaagttagaggtaaatggtattttaatatcatcggtaaagattctgaacttgttgaaaaagctgaacatcttttacgtgatatgggatgggaagatgaatgcgatggatgtcctctttatgaagacggagaaagcgcaggattttggatttaccattctgacgtcgagcagtttaaagctgattggaaaattgtgaaaaagtctgtttga"
	if feature != seq {
		t.Errorf("Feature GetSequence method has failed. Got this:\n%s instead of \n%s", feature, seq)
	}

}

func TestLocationParser(t *testing.T) {
	gbk, _ := Read("../../data/t4_intron.gb")

	// Read 1..243
	feature := gbk.Features[1].GetSequence()
	seq := "atgagattacaacgccagagcatcaaagattcagaagttagaggtaaatggtattttaatatcatcggtaaagattctgaacttgttgaaaaagctgaacatcttttacgtgatatgggatgggaagatgaatgcgatggatgtcctctttatgaagacggagaaagcgcaggattttggatttaccattctgacgtcgagcagtttaaagctgattggaaaattgtgaaaaagtctgtttga"
	if feature != seq {
		t.Errorf("Feature sequence parser has changed on test '1..243'. Got this:\n%s instead of \n%s", feature, seq)
	}

	// Read join(893..1441,2459..2770)
	featureJoin := gbk.Features[6].GetSequence()
	seqJoin := "atgaaacaataccaagatttaattaaagacatttttgaaaatggttatgaaaccgatgatcgtacaggcacaggaacaattgctctgttcggatctaaattacgctgggatttaactaaaggttttcctgcggtaacaactaagaagctcgcctggaaagcttgcattgctgagctaatatggtttttatcaggaagcacaaatgtcaatgatttacgattaattcaacacgattcgttaatccaaggcaaaacagtctgggatgaaaattacgaaaatcaagcaaaagatttaggataccatagcggtgaacttggtccaatttatggaaaacagtggcgtgattttggtggtgtagaccaaattatagaagttattgatcgtattaaaaaactgccaaatgataggcgtcaaattgtttctgcatggaatccagctgaacttaaatatatggcattaccgccttgtcatatgttctatcagtttaatgtgcgtaatggctatttggatttgcagtggtatcaacgctcagtagatgttttcttgggtctaccgtttaatattgcgtcatatgctacgttagttcatattgtagctaagatgtgtaatcttattccaggggatttgatattttctggtggtaatactcatatctatatgaatcacgtagaacaatgtaaagaaattttgaggcgtgaacctaaagagctttgtgagctggtaataagtggtctaccttataaattccgatatctttctactaaagaacaattaaaatatgttcttaaacttaggcctaaagatttcgttcttaacaactatgtatcacaccctcctattaaaggaaagatggcggtgtaa"
	if featureJoin != seqJoin {
		t.Errorf("Feature sequence parser has changed on test 'join(893..1441,2459..2770)'. Got this:\n%s instead of \n%s", featureJoin, seqJoin)
	}

	// Read complement(2791..3054)
	featureComplement := gbk.Features[10].GetSequence()
	seqComplement := "ttattcactacccggcatagacggcccacgctggaataattcgtcatattgtttttccgttaaaacagtaatatcgtagtaacagtcagaagaagttttaactgtggaaattttattatcaaaatactcacgagtcattttatgagtatagtattttttaccataaatggtaataggctgttctggtcctggaacttctaactcgcttgggttaggaagtgtaaaaagaactacaccagaagtatctttaaatcgtaaaatcat"
	if featureComplement != seqComplement {
		t.Errorf("Feature sequence parser has changed on test 'complement(2791..3054)'. Got this:\n%s instead of \n%s", featureComplement, seqComplement)
	}

	// Read join(complement(315..330),complement(339..896))
	// Note: it is known that some software, like Snapgene, assumes that since both strands are in the reverse direction
	// that the first sequence should be appended to the reverse sequence, instead of the second sequence
	// getting appended to the first. Biopython appends the second sequence to the first, and that is logically
	// the most obvious thing to do, so we are implementing it that way.
	featureJoinComplement := gbk.Features[3].GetSequence()
	seqJoinComplement := "ataccaatttaatcattcatttatatactgattccgtaagggttgttacttcatctattttataccaatgcgtttcaaccatttcacgcttgcttatatcatcaagaaaacttgcgtctaattgaactgttgaattaacacgatgccttttaacgatgcgagaaacaactacttcatctgcataaggtaatgcagcatataacagagcaggcccgccaattacacttactttagaattctgatcaagcatagtttcgaatggtgcattagggcttgacacttgaatttcgccgccagaaatgtaagttatatattgctcccaagtaatatagaaatgtgctaaatcgccgtctttagttacaggataatcacgcgcaaggtcacacaccacaatatggctacgaccaggaagtaatgtaggcaatgactggaacgttttagcacccataatcataattgtgccttcagtacgagctttaaaattctggaggtcctttttaactcgtccccatggtaaaccatcacctaaaccgaatgctaattcattaaagccgtcgaccgttttagttggaga"
	if featureJoinComplement != seqJoinComplement {
		t.Errorf("Feature sequence parser has changed on test 'join(complement(315..330),complement(339..896))'. Got this:\n%s instead of \n%s", featureJoinComplement, seqJoinComplement)
	}

	// Read complement(join(893..1098,1101..2770))
	featureComplementJoin := gbk.Features[5].GetSequence()
	seqComplementJoin := "ttacaccgccatctttcctttaataggagggtgtgatacatagttgttaagaacgaaatctttaggcctaagtttaagaacatattttaattgttctttagtagaaagatatcggaatttataaggtagaccacttattaccagctcacaaagctctttaggttcacgcctcaaaatttctttacattgttctacgtgattcatatagatatgagtattaccaccagaaaatatcaaatcccctggaataagattacacatcttagctacaatatgaactaacgtagcatatgacgcaatattaaacggtagcattatgttcagataaggtcgttaatcttaccccggaattatatccagctgcatgtcaccatgcagagcagactatatctccaacttgttaaagcaagttgtctatcgtttcgagtcacttgaccctactccccaaagggatagtcgttaggcatttatgtagaaccaattccatttatcagattttacacgataagtaactaatccagacgaaattttaaaatgtctagctgcatctgctgcacaatcaaaaataaccccatcacatgaaatctttttaatattactaggctttttacctttcatcttttctgatattttagatttagttatgtctgaatgcttatgattaaagaatgaattattttcacctgaacgatttctgcatttactacaagtataagcagaagtttgtatgcgaacaccgcacttacaaaacttatgggtttctggattccaacgcccgtttttacttccgggtttactgtaaagagctttccgaccatcaggtccaagtttaagcatcttagctttaacagtttcagaacgtttcttaataatttcttcttttaatggatgcgtagaacatgtatcaccaaacgttgcatcagcaatattgtatccattaattttagaattaagctctttaatccaaaaattttctcgttcaataatcaaatctttctcatatggaatttcttccaaaatagaacattcaaacacattaccatgtttgttaaaagacctctgaagttttatagaagaatggcatcctttttctaaatctttaaaatgcctcttccatctcttttcaaaatctttagcacttcctacatatactttattgtttaaagtatttttaatctgataaattccgcttttcataaatacctctttaaatatagaagtatttattaaagggcaagtcctacaatttagcacgggattgtctactagagaggttccccgtttagatagattacaagtataagtcaccttatactcaggcctcaattaacccaagaaaacatctactgagcgttgataccactgcaaatccaaatagccattacgcacattaaactgatagaacatatgacaaggcggtaatgccatatatttaagttcagctggattccatgcagaaacaatttgacgcctatcatttggcagttttttaatacgatcaataacttctataatttggtctacaccaccaaaatcacgccactgttttccataaattggaccaagttcaccgctatggtatcctaaatcttttgcttgattttcgtaattttcatcccagactgttttgccttggattaacgaatcgtgttgaattaatcgtaaatcatacatttgtgcttcctgataaaaaccatattagctcagcaatgcaagctttccaggcgagcttcttagttgttaccgcaggaaaacctttagttaaatcccagcgtaatttagatccgaacagagcaattgttcctgtgcctgtacgatcatcggtttcataaccattttcaaaaatgtctttaattaaatcttggtattgtttcat"
	if featureComplementJoin != seqComplementJoin {
		t.Errorf("Feature sequence parser has changed on test 'complement(join(893..1098,1101..2770))'. Got this:\n%s instead of \n%s", featureComplementJoin, seqComplementJoin)
	}
}

func TestGenbankNewlineParsingRegression(t *testing.T) {
	gbk, _ := Read("../../data/puc19.gbk")

	for _, feature := range gbk.Features {
		if feature.Location.Start == 410 && feature.Location.End == 1750 && feature.Type == "CDS" {
			if feature.Attributes["product"] != "chromosomal replication initiator informational ATPase" {
				t.Errorf("Newline parsing has failed.")
			}
			break
		}
	}
}

func TestLocusParseRegression(t *testing.T) {
	gbk, _ := Read("../../data/puc19.gbk")
	staticGbk, _ := Read("../../data/puc19static.gbk")

	if diff := cmp.Diff(gbk, staticGbk, cmpopts.IgnoreFields(Feature{}, "ParentSequence")); diff != "" {
		t.Errorf("The meta parser has changed behaviour. Got this diff:\n%s", diff)
	}
}

func BenchmarkRead(b *testing.B) {
	for i := 0; i < b.N; i++ {
		_, _ = Read("../../data/bsub.gbk")
	}
}

func BenchmarkRead1(b *testing.B)     { BenchmarkRead(b) }
func BenchmarkRead10(b *testing.B)    { BenchmarkRead(b) }
func BenchmarkRead100(b *testing.B)   { BenchmarkRead(b) }
func BenchmarkRead1000(b *testing.B)  { BenchmarkRead(b) }
func BenchmarkRead10000(b *testing.B) { BenchmarkRead(b) }

/******************************************************************************

Gbk/gb/genbank related benchmarks end here.

******************************************************************************/

/******************************************************************************

GbkMulti/GbkFlat related tests begin here.

******************************************************************************/

func ExampleReadMulti() {
	sequences := ReadMulti("../../data/multiGbk_test.seq")
	var locus []string
	for _, sequence := range sequences {
		locus = append(locus, sequence.Meta.Locus.Name)
	}

	fmt.Println(strings.Join(locus, ", "))
	// Output: AB000100, AB000106
}

func ExampleReadFlat() {
	sequences := ReadFlat("../../data/long_comment.seq")
	var locus []string
	for _, sequence := range sequences {
		locus = append(locus, sequence.Meta.Locus.Name)
	}

	fmt.Println(strings.Join(locus, ", "))
	// Output: AB000100, AB000106
}

func ExampleReadFlatGz() {
	sequences := ReadFlatGz("../../data/flatGbk_test.seq.gz")
	//sequences := ReadFlatGz("../../data/gbbct358.seq.gz")
	var locus []string
	for _, sequence := range sequences {
		locus = append(locus, sequence.Meta.Locus.Name)
	}
	fmt.Println(strings.Join(locus, ", "))
	// Output: AB000100, AB000106
}

func ExampleParseMulti() {
	file, _ := ioutil.ReadFile("../../data/multiGbk_test.seq")
	sequences := ParseMulti(file)
	var locus []string
	for _, sequence := range sequences {
		locus = append(locus, sequence.Meta.Locus.Name)
	}

	fmt.Println(strings.Join(locus, ", "))
	// Output: AB000100, AB000106
}

func ExampleParseFlat() {
	file, _ := ioutil.ReadFile("../../data/flatGbk_test.seq")
	sequences := ParseFlat(file)
	var locus []string
	for _, sequence := range sequences {
		locus = append(locus, sequence.Meta.Locus.Name)
	}

	fmt.Println(strings.Join(locus, ", "))
	// Output: AB000100, AB000106
}

/******************************************************************************

GbkMulti/GbkFlat related tests end here.

******************************************************************************/

func TestBenchlingGenbank(t *testing.T) {
	sequence, _ := Read("../../data/benchling.gb")

	if len(sequence.Features) != 17 {
		t.Errorf("Parsing benchling genbank file not returned the correct quantity of features")
	}
}

func TestParse(t *testing.T) {
	genBank, err := Parse([]byte(`ORIGIN
ORIGIN`))
	assert.NotNil(t, genBank)
	assert.Nil(t, err)
}

func TestRead_error(t *testing.T) {
	readErr := errors.New("read error")
	oldReadFileFn := readFileFn
	readFileFn = func(filename string) ([]byte, error) {
		return nil, readErr
	}
	defer func() {
		readFileFn = oldReadFileFn
	}()
	_, err := Read("/tmp/file")
	assert.EqualError(t, err, readErr.Error())
}

func Test_getReference(t *testing.T) {
	ref := getReference([]string{
		"test a b c",
		"test g g g",
	}, []string{
		"REMARK something",
		"REMARK something else",
	})
	assert.NotNil(t, ref)
}

func Test_getSequence(t *testing.T) {
	var err error
	regexErr := errors.New("regex error")
	oldLogFatalFn := logFatalFn
	oldRegexpCompileFn := regexpCompileFn
	regexpCompileFn = func(expr string) (*regexp.Regexp, error) {
		r, _ := regexp.Compile("[^a-zA-Z]+")
		return r, regexErr
	}
	logFatalFn = func(v ...interface{}) {
		err = v[0].(error)
	}
	defer func() {
		logFatalFn = oldLogFatalFn
		regexpCompileFn = oldRegexpCompileFn
	}()
	getSequence([]string{"test"})
	assert.EqualError(t, err, regexErr.Error())
}

func TestBuildFeatureString(t *testing.T) {
	feature := Feature{
		Type:        "test type",
		Description: "a description",
		Location: Location{
			GbkLocationString: "gbk location",
		},
	}
	str := BuildFeatureString(feature)
	assert.Equal(t, str, "     test type       gbk location\n")
}

func Test_getSourceOrganism(t *testing.T) {
	for _, tc := range []struct {
		splitLine    []string
		subLines     []string
		wantSource   string
		wantOrganism string
	}{
		{
			splitLine:    []string{"test"},
			subLines:     []string{"subline 1", "subline 2"},
			wantSource:   "",
			wantOrganism: "1",
		},
		{
			splitLine:    []string{"test"},
			subLines:     []string{" ", "subline 2"},
			wantSource:   "",
			wantOrganism: "2",
		},
		{
			splitLine:    []string{"test"},
			subLines:     []string{" ", "subline 2"},
			wantSource:   "",
			wantOrganism: "2",
		},
		{
			splitLine:    []string{" something else "},
			subLines:     []string{" ORGANISM", "subline 2", "subline 3"},
			wantSource:   "",
			wantOrganism: "",
		},
	} {
		source, organism := getSourceOrganism(tc.splitLine, tc.subLines)
		assert.Equal(t, tc.wantSource, source)
		assert.Equal(t, tc.wantOrganism, organism)
	}
}

func Test_getAttributeValue(t *testing.T) {
	for _, tc := range []struct {
		attributeSplit     []string
		wantAttributeValue string
	}{
		{
			attributeSplit:     []string{"test"},
			wantAttributeValue: "",
		},
		{
			attributeSplit:     []string{"test", "something", "abcdef"},
			wantAttributeValue: "something",
		},
	} {
		attributeValue := getAttributeValue(tc.attributeSplit)
		assert.Equal(t, tc.wantAttributeValue, attributeValue)
	}
}

func Test_parseComplicatedJoin(t *testing.T) {
	for _, tc := range []struct {
		expression   string
		location     *Location
		wantLocation *Location
	}{
		{
			expression: "complement(12..34),complement(56..78)",
			location:   &Location{},
			wantLocation: &Location{
				SubLocations: []Location{
					{Start: 11, End: 34, Complement: true, Join: false, FivePrimePartial: false, ThreePrimePartial: false, GbkLocationString: "", SubLocations: []Location(nil)},
					{Start: 55, End: 78, Complement: true, Join: false, FivePrimePartial: false, ThreePrimePartial: false, GbkLocationString: "", SubLocations: []Location(nil)},
				},
			},
		},
		{
			expression: "complement(complement(12..34)),complement(56..78)",
			location:   &Location{},
			wantLocation: &Location{
				SubLocations: []Location{
					{Start: 11, End: 34, Complement: true, Join: false, FivePrimePartial: false, ThreePrimePartial: false, GbkLocationString: "", SubLocations: []Location(nil)},
					{Start: 55, End: 78, Complement: true, Join: false, FivePrimePartial: false, ThreePrimePartial: false, GbkLocationString: "", SubLocations: []Location(nil)},
				},
			},
		},
	} {
		parseComplicatedJoin(tc.expression, tc.location)
		assert.Equal(t, tc.wantLocation, tc.location)
	}
}
