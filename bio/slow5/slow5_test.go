package slow5

import (
	"errors"
	"io"
	"io/ioutil"
	"os"
	"strings"
	"testing"

	"github.com/google/go-cmp/cmp"
)

const maxLineSize = 2 * 32 * 1024

func TestParse(t *testing.T) {
	file, err := os.Open("data/example.slow5")
	if err != nil {
		t.Errorf("Failed to open example.slow5: %s", err)
	}
	parser, err := NewParser(file, maxLineSize)
	if err != nil {
		t.Errorf("Failed to parse headers of file: %s", err)
	}
	// Test headers
	headers, _ := parser.Header()
	if len(headers.HeaderValues) != 1 {
		t.Errorf("There should only be 1 read group. Got: %d", len(headers.HeaderValues))
	}
	if headers.HeaderValues[0].Attributes["@asic_id"] != "4175987214" {
		t.Errorf("Expected AsicId 4175987214. Got: %s", headers.HeaderValues[0].Attributes["asic_id"])
	}

	// Test reads
	var outputReads []Read
	for {
		read, err := parser.Next()
		if err != nil {
			if !errors.Is(err, io.EOF) {
				t.Errorf("Got unknown error: %s", err)
			}
			break
		}
		outputReads = append(outputReads, *read)
	}
	if outputReads[0].RawSignal[0] != 430 {
		t.Errorf("Expected first outputRead to have a raw_signal of 430. Got: %d", outputReads[0].RawSignal[0])
	}
}

func TestParseImproperHeaders(t *testing.T) {
	// Improper files
	file, err := os.Open("data/header_tests/test_header_without_tabs.slow5")
	if err != nil {
		t.Errorf("Failed to open file with error: %s", err)
	}
	_, err = NewParser(file, maxLineSize)
	if err == nil {
		t.Errorf("Test should have failed if header line doesn't have any tabs")
	}

	file, err = os.Open("data/header_tests/test_header_numReadGroups_bad.slow5")
	if err != nil {
		t.Errorf("Failed to open file with error: %s", err)
	}
	_, err = NewParser(file, maxLineSize)
	if err == nil {
		t.Errorf("Test should have failed if numReadGroup can't be converted to an int")
	}

	file, err = os.Open("data/header_tests/test_header_not_enough_attributes.slow5")
	if err != nil {
		t.Errorf("Failed to open file with error: %s", err)
	}
	_, err = NewParser(file, maxLineSize)
	if err == nil {
		t.Errorf("Test should have failed if the header doesn't have enough attributes for numReadGroup")
	}

	file, err = os.Open("data/header_tests/test_header_empty.slow5")
	if err != nil {
		t.Errorf("Failed to open file with error: %s", err)
	}
	_, err = NewParser(file, maxLineSize)
	if err == nil {
		t.Errorf("Test should have failed if the file is empty")
	}
}

func testParseReadsHelper(t *testing.T, fileTarget string, errorMessage string) {
	file, err := os.Open(fileTarget)
	if err != nil {
		t.Errorf("Failed to open file with error: %s", err)
	}
	parser, _ := NewParser(file, maxLineSize)
	var targetErr []error
	for {
		_, err := parser.Next()
		if err != nil {
			if !errors.Is(err, io.EOF) {
				targetErr = append(targetErr, err)
			}
			break
		}
	}
	if len(targetErr) == 0 {
		t.Errorf(errorMessage)
	}
}

func TestParseReads(t *testing.T) {
	file, err := os.Open("data/example.slow5")
	if err != nil {
		t.Errorf("Failed to open file with error: %s", err)
	}
	parser, _ := NewParser(file, maxLineSize)

	var outputReads []Read
	for {
		read, err := parser.Next()
		if err != nil {
			if !errors.Is(err, io.EOF) {
				t.Errorf("Got unknown error: %s", err)
			}
			break
		}
		outputReads = append(outputReads, *read)
	}
	if outputReads[0].ReadID != "0026631e-33a3-49ab-aa22-3ab157d71f8b" {
		t.Errorf("First read id should be 0026631e-33a3-49ab-aa22-3ab157d71f8b. Got: %s", outputReads[0].ReadID)
	}

	// Test improper files
	testParseReadsHelper(t, "data/read_tests/endReason.slow5", "Test should have failed if there are unknown end reasons")
	testParseReadsHelper(t, "data/read_tests/continue.slow5", "Test should have failed at terminate, but should have gone through a continue")
	testParseReadsHelper(t, "data/read_tests/read_group.slow5", "Test should have failed with bad read_group")
	testParseReadsHelper(t, "data/read_tests/digitisation.slow5", "Test should have failed with bad digitisation")
	testParseReadsHelper(t, "data/read_tests/offset.slow5", "Test should have failed with bad offset")
	testParseReadsHelper(t, "data/read_tests/range.slow5", "Test should have failed with bad range")
	testParseReadsHelper(t, "data/read_tests/sampling_rate.slow5", "Test should have failed with bad samping_rate")
	testParseReadsHelper(t, "data/read_tests/len_raw_signal.slow5", "Test should have failed with bad len_raw_signal")
	testParseReadsHelper(t, "data/read_tests/raw_signal.slow5", "Test should have failed with bad raw_signal")
	testParseReadsHelper(t, "data/read_tests/start_time.slow5", "Test should have failed with bad start_time")
	testParseReadsHelper(t, "data/read_tests/read_number.slow5", "Test should have failed with bad read_number")
	testParseReadsHelper(t, "data/read_tests/start_mux.slow5", "Test should have failed with bad start_mux")
	testParseReadsHelper(t, "data/read_tests/median_before.slow5", "Test should have failed with bad median_before")
	testParseReadsHelper(t, "data/read_tests/end_reason.slow5", "Test should have failed with if end reason can't be converted to int")
	testParseReadsHelper(t, "data/read_tests/end_reason_unknown.slow5", "Test should have failed with end reason out of range")
	testParseReadsHelper(t, "data/read_tests/unknown.slow5", "Test should have failed with unknown header")
}

func TestWrite(t *testing.T) {
	file, err := os.Open("data/example.slow5")
	if err != nil {
		t.Errorf("Failed to open file with error: %s", err)
	}
	// Parse headers
	parser, err := NewParser(file, maxLineSize)
	if err != nil {
		t.Errorf("Failed to parse headers with error: %s", err)
	}
	// Create empty file
	testFile, err := os.Create("data/test_write.slow5")
	if err != nil {
		t.Errorf("Failed to create temporary file")
	}

	// Send parsed reads into this channel
	reads := make(chan Read)
	go func() {
		for {
			read, err := parser.Next()
			if err != nil {
				// Break at EOF
				break
			}
			reads <- *read
		}
		close(reads)
	}()

	// Write header
	headers, _ := parser.Header()
	_, err = headers.WriteTo(testFile)
	if err != nil {
		t.Errorf("Failed to write slow5 header. Got error: %s", err)
	}
	for read := range reads {
		_, err = read.WriteTo(testFile)
		if err != nil {
			t.Errorf("Failed to write slow5 read. Got error: %s", err)
		}
	}

	// Compare both files
	example, err := ioutil.ReadFile("data/example.slow5")
	if err != nil {
		t.Errorf("Failed to read example file: %s", err)
	}
	testWrite, err := ioutil.ReadFile("data/test_write.slow5")
	if err != nil {
		t.Errorf("Failed to read test file: %s", err)
	}
	os.Remove("data/test_write.slow5")

	if string(example) != string(testWrite) {
		t.Errorf("Example and test write are different")
	}
}

func TestSimpleExample(t *testing.T) {
	file := strings.NewReader(`#slow5_version	0.2.0
#num_read_groups	1
@asic_id	4175987214
#char*	uint32_t	double	double	double	double	uint64_t	int16_t*	uint64_t	int32_t	uint8_t	double	enum{unknown,partial,mux_change,unblock_mux_change,data_service_unblock_mux_change,signal_positive,signal_negative}	char*
#read_id	read_group	digitisation	offset	range	sampling_rate	len_raw_signal	raw_signal	start_time	read_number	start_mux	median_before	end_reason	channel_number
0026631e-33a3-49ab-aa22-3ab157d71f8b	0	8192	16	1489.52832	4000	5347	430,472,463	8318394	5383	1	219.133423	5	10`)
	parser, err := NewParser(file, maxLineSize)
	if err != nil {
		t.Errorf("failed: %s", err)
	}
	read, _ := parser.Next()
	if read.RawSignal[0] != 430 {
		t.Errorf("Should have gotten 430. Got: %d", read.RawSignal[0])
	}
}

func TestSvb(t *testing.T) {
	file, _ := os.Open("data/example.slow5")
	defer file.Close()
	const maxLineSize = 2 * 32 * 1024
	parser, _ := NewParser(file, maxLineSize)
	var outputReads []Read
	for {
		read, err := parser.Next()
		if err != nil {
			if !errors.Is(err, io.EOF) {
				t.Errorf("Got unknown error: %s", err)
			}
			break
		}
		outputReads = append(outputReads, *read)
	}

	for readNum, read := range outputReads {
		rawSignal := read.RawSignal
		mask, data := SvbCompressRawSignal(rawSignal)
		rawSignalDecompressed := SvbDecompressRawSignal(len(rawSignal), mask, data)
		if !cmp.Equal(rawSignal, rawSignalDecompressed) {
			t.Errorf("Read signal at readNum %d (%v) didn't match decompressed signal %v", readNum, rawSignal, rawSignalDecompressed)
		}
	}
}
