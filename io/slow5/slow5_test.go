package slow5

import (
	"os"
	"testing"
)

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
}

func TestParseHeader(t *testing.T) {
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
}
