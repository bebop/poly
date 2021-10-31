package slow5

import (
	"os"
	"testing"
)

func testParseReadsHelper(t *testing.T, fileTarget string, errorMessage string) {
	file, err := os.Open(fileTarget)
	if err != nil {
		t.Errorf("Failed to open file with error: %s", err)
	}
	reads := make(chan Read, 1000)
	errs := make(chan error, 1000)
	go ParseReads(file, reads, errs)
	var targetErr []error
	for err := range errs {
		targetErr = append(targetErr, err)
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
	reads := make(chan Read, 1000)
	errs := make(chan error, 1000)

	go ParseReads(file, reads, errs)

	var outputReads []Read
	for read := range reads {
		outputReads = append(outputReads, read)
	}

	for err := range errs {
		if err != nil {
			t.Errorf("Failed ParseReads with error: %s", err)
		}
	}

	if outputReads[0].ReadId != "0026631e-33a3-49ab-aa22-3ab157d71f8b" {
		t.Errorf("First read id should be 0026631e-33a3-49ab-aa22-3ab157d71f8b. Got: %s", outputReads[0].ReadId)
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

func TestParseHeader(t *testing.T) {
	// Proper file
	file, err := os.Open("data/example.slow5")
	if err != nil {
		t.Errorf("Failed to open file with error: %s", err)
	}
	readGroups, err := ParseHeader(file)
	if err != nil {
		t.Errorf("Failed ParseReadGroups with error: %s", err)
	}

	if len(readGroups) != 1 {
		t.Errorf("There should only be 1 read group. Got: %d", len(readGroups))
	}

	if readGroups[0].Attributes["@asic_id"] != "4175987214" {
		t.Errorf("Expected AsicId 4175987214. Got: %s", readGroups[0].Attributes["asic_id"])
	}

	// Improper files
	file, err = os.Open("data/header_tests/test_header_without_tabs.slow5")
	if err != nil {
		t.Errorf("Failed to open file with error: %s", err)
	}
	_, err = ParseHeader(file)
	if err == nil {
		t.Errorf("Test should have failed if header line doesn't have any tabs")
	}

	file, err = os.Open("data/header_tests/test_header_numReadGroups_bad.slow5")
	if err != nil {
		t.Errorf("Failed to open file with error: %s", err)
	}
	_, err = ParseHeader(file)
	if err == nil {
		t.Errorf("Test should have failed if numReadGroup can't be converted to an int")
	}

	file, err = os.Open("data/header_tests/test_header_not_enough_attributes.slow5")
	if err != nil {
		t.Errorf("Failed to open file with error: %s", err)
	}
	_, err = ParseHeader(file)
	if err == nil {
		t.Errorf("Test should have failed if the header doesn't have enough attributes for numReadGroup")
	}

	file, err = os.Open("data/header_tests/test_header_no_termination.slow5")
	if err != nil {
		t.Errorf("Failed to open file with error: %s", err)
	}
	_, err = ParseHeader(file)
	if err == nil {
		t.Errorf("Test should have failed if the the header doesn't terminate")
	}
}
