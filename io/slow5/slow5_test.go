package slow5

import (
	"errors"
	"io"
	"os"
	"testing"
)

const maxLineSize = 2 * 32 * 1024

func TestParse(t *testing.T) {
	file, err := os.Open("data/example.slow5")
	if err != nil {
		t.Errorf("Failed to open example.slow5: %s", err)
	}
	parser, headers, err := NewParser(file, maxLineSize)
	if err != nil {
		t.Errorf("Failed to parse headers of file: %s", err)
	}
	// Test headers
	if len(headers) != 1 {
		t.Errorf("There should only be 1 read group. Got: %d", len(headers))
	}
	if headers[0].Attributes["@asic_id"] != "4175987214" {
		t.Errorf("Expected AsicId 4175987214. Got: %s", headers[0].Attributes["asic_id"])
	}

	// Test reads
	var outputReads []Read
	for {
		read, err := parser.ParseNext()
		if err != nil {
			if !errors.Is(err, io.EOF) {
				t.Errorf("Got unknown error: %s", err)
			}
			break
		}
		outputReads = append(outputReads, read)
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
	_, _, err = NewParser(file, maxLineSize)
	if err == nil {
		t.Errorf("Test should have failed if header line doesn't have any tabs")
	}

	file, err = os.Open("data/header_tests/test_header_numReadGroups_bad.slow5")
	if err != nil {
		t.Errorf("Failed to open file with error: %s", err)
	}
	_, _, err = NewParser(file, maxLineSize)
	if err == nil {
		t.Errorf("Test should have failed if numReadGroup can't be converted to an int")
	}

	file, err = os.Open("data/header_tests/test_header_not_enough_attributes.slow5")
	if err != nil {
		t.Errorf("Failed to open file with error: %s", err)
	}
	_, _, err = NewParser(file, maxLineSize)
	if err == nil {
		t.Errorf("Test should have failed if the header doesn't have enough attributes for numReadGroup")
	}

	file, err = os.Open("data/header_tests/test_header_empty.slow5")
	if err != nil {
		t.Errorf("Failed to open file with error: %s", err)
	}
	_, _, err = NewParser(file, maxLineSize)
	if err == nil {
		t.Errorf("Test should have failed if the file is empty")
	}
}

func testParseReadsHelper(t *testing.T, fileTarget string, errorMessage string) {
	file, err := os.Open(fileTarget)
	if err != nil {
		t.Errorf("Failed to open file with error: %s", err)
	}
	parser, _, _ := NewParser(file, maxLineSize)
	var targetErr []error
	for {
		read, err := parser.ParseNext()
		if err != nil {
			if !errors.Is(err, io.EOF) {
				t.Errorf("Got unknown error: %s", err)
			}
			break
		}
		err = read.Error
		if err != nil {
			targetErr = append(targetErr, err)
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
	parser, _, _ := NewParser(file, maxLineSize)

	var outputReads []Read
	for {
		read, err := parser.ParseNext()
		if err != nil {
			if !errors.Is(err, io.EOF) {
				t.Errorf("Got unknown error: %s", err)
			}
			break
		}
		err = read.Error
		if err != nil {
			t.Errorf("Failed ParseReads with error: %s", err)
		}
		outputReads = append(outputReads, read)
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
	parser, headers, err := NewParser(file, maxLineSize)
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
			read, err := parser.ParseNext()
			if err != nil {
				// Break at EOF
				break
			}
			reads <- read
		}
		close(reads)
	}()

	// Write
	err = Write(headers, reads, testFile)
	if err != nil {
		t.Errorf("Failed to write slow5 file. Got error: %s", err)
	}

	// Compare both files
	example, err := os.ReadFile("data/example.slow5")
	if err != nil {
		t.Errorf("Failed to read example file: %s", err)
	}
	testWrite, err := os.ReadFile("data/test_write.slow5")
	if err != nil {
		t.Errorf("Failed to read test file: %s", err)
	}
	os.Remove("data/test_write.slow5")

	if string(example) != string(testWrite) {
		t.Errorf("Example and test write are different")
	}
}
