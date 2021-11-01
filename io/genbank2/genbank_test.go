package genbank

import (
	"os"
	"testing"
)

func TestParse(t *testing.T) {
	file, err := os.Open("data/puc19.gbk")
	if err != nil {
		t.Errorf("Failed to open puc19.gbk test file. Got error: %s", err)
	}

	_, err = Parse(file)
	if err != nil {
		t.Errorf("%s", err)
	}
}
