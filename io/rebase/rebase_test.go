package rebase

import (
	"errors"
	"io"
	"os"
	"strings"
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestRead(t *testing.T) {
	_, err := Read("data/FAKE.txt")
	if err == nil {
		t.Errorf("Failed to error on fake file")
	}
}

func TestParse(t *testing.T) {
	file, err := os.Open("data/rebase_minimal.txt")
	if err != nil {
		t.Errorf("Failed to read rebase test file")
	}
	emzymeMap, err := parseFn(file)
	if err != nil {
		t.Errorf("Error parsing the file data")
	}
	exportedJSON, err := Export(emzymeMap)
	if err != nil {
		t.Errorf("Error exporing the enzyme map to JSON")
	}
	file2, err := os.Create("rebase_minimal.json")
	file2.WriteString(string(exportedJSON))
	file2.Close()
	
}

func TestParse_error(t *testing.T) {
	parseErr := errors.New("fake error")
	oldReadAllFn := readAllFn
	readAllFn = func(r io.Reader) ([]byte, error) {
		return nil, parseErr
	}
	defer func() {
		readAllFn = oldReadAllFn
	}()
	_, err := Parse(strings.NewReader(""))
	assert.EqualError(t, err, parseErr.Error())
}

func TestRead_error(t *testing.T) {
	readErr := errors.New("fake error")
	oldParseFn := parseFn
	parseFn = func(file io.Reader) (map[string]Enzyme, error) {
		return nil, readErr
	}
	defer func() {
		parseFn = oldParseFn
	}()
	_, err := Read("data/rebase_test.txt")
	assert.EqualError(t, err, readErr.Error())

}

func TestExport_error(t *testing.T) {
	exportErr := errors.New("fake error")
	oldMarshallFn := marshallFn
	marshallFn = func(v any) ([]byte, error) {
		return []byte{}, exportErr
	}
	defer func() {
		marshallFn = oldMarshallFn
	}()
	_, err := Export(map[string]Enzyme{})
	assert.EqualError(t, err, exportErr.Error())
}
