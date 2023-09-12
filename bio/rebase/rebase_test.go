package rebase

import (
	"errors"
	"io"
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

func TestRead_multipleRefs(t *testing.T) {
	enzymeMap, err := Read("data/rebase_test.txt")
	if err != nil {
		t.Error("Failed to read test file")
	}

	assert.Equal(t,
		"Calleja, F., de Waard, A., Unpublished observations.\nHughes, S.G., Bruce, T., Murray, K., Unpublished observations.",
		enzymeMap["AcaI"].References)
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
