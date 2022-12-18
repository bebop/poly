package rebase

import (
	"encoding/json"
	"errors"
	"io"
	"io/ioutil"
	"os"
	"reflect"
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

	// Open the json version of the rebase file to compare against
	rebaseJSONFile, err := os.Open("data/rebase_minimal.json")
	if err != nil {
		t.Errorf("Failed to read rebase test file")
	}

	// Parse JSON file into the enzyme maps
	// read our opened jsonFile as a byte array.
	byteData, _ := ioutil.ReadAll(rebaseJSONFile)

	// we initialize our enzymes map
	var benchmarkEnzymeMap map[string]Enzyme

	// we unmarshal our byteArray which contains our
	// jsonFile's content into 'enzymes' which we defined above
	err = json.Unmarshal(byteData, &benchmarkEnzymeMap)
	if err != nil {
		t.Errorf("Failed to parse rebase json file into enzyme map")
	}

	// Open the rebase file
	rebaseFile, err := os.Open("data/rebase_minimal.txt")
	if err != nil {
		t.Errorf("Failed to read rebase test file")
	}

	
	// Parse the rebase file
	emzymeMap, err := parseFn(rebaseFile)
	if err != nil {
		t.Errorf("Error parsing the file data")
	}

	// Compare the enzyme maps
	eq := reflect.DeepEqual(emzymeMap, benchmarkEnzymeMap)
	if !eq {
		t.Errorf("Parsed enzyme map does not match benchmark")
	}

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
