package uniprot

import (
	"compress/gzip"
	"encoding/xml"
	"errors"
	"fmt"
	"os"
	"testing"

	"github.com/stretchr/testify/assert"
)

func ExampleRead() {
	entries, _, _ := Read("data/uniprot_sprot_mini.xml.gz")

	var entry Entry
	for singleEntry := range entries {
		entry = singleEntry
	}
	fmt.Println(entry.Accession[0])
	// Output: O55723
}

func ExampleParse() {
	xmlFile, _ := os.Open("data/uniprot_sprot_mini.xml.gz")
	unzippedBytes, _ := gzip.NewReader(xmlFile)

	entries := make(chan Entry, 100) // if you don't have a buffered channel, nothing will be read in loops on the channel.
	decoderErrors := make(chan error, 100)
	decoder := xml.NewDecoder(unzippedBytes)
	go Parse(decoder, entries, decoderErrors)

	var entry Entry
	for singleEntry := range entries {
		entry = singleEntry
	}
	fmt.Println(entry.Accession[0])
	// Output: O55723
}

func TestRead(t *testing.T) {
	_, _, err := Read("data/test")
	if err == nil {
		t.Errorf("Failed to fail on non-gzipped file")
	}

	_, _, err = Read("data/FAKE")
	if err == nil {
		t.Errorf("Failed to fail on empty file")
	}

	_, errors, err := Read("data/uniprot_sprot_mini.xml.gz")
	if err != nil {
		t.Errorf("Failed on real file with error: %v", err)
	}

	for err := range errors {
		if err != nil {
			t.Errorf("Failed during parsing with error: %v", err)
		}
	}
}

func TestParse(t *testing.T) {
	t.Run("error getting a token", func(t *testing.T) {
		entries := make(chan Entry, 100)
		decoderErrors := make(chan error, 100)
		tokenErr := errors.New("token error")
		firstRun := true
		decoder := &mockDecoder{
			TokenFn: func() (xml.Token, error) {
				if firstRun {
					firstRun = false
					return nil, tokenErr
				}
				return nil, errors.New("EOF")
			},
		}
		Parse(decoder, entries, decoderErrors)
		assert.EqualError(t, <-decoderErrors, tokenErr.Error())
	})

	t.Run("error decoding after getting a token", func(t *testing.T) {
		entries := make(chan Entry, 100)
		decoderErrors := make(chan error, 100)
		decodeErr := errors.New("decode error")
		startElement := xml.StartElement{
			Name: xml.Name{
				Local: "entry",
			},
			Attr: nil,
		}
		firstRun := true
		decoder := &mockDecoder{
			DecodeElementFn: func(v interface{}, start *xml.StartElement) error {
				return decodeErr
			},
			TokenFn: func() (xml.Token, error) {
				if firstRun {
					firstRun = false
					return startElement, nil
				}
				return nil, errors.New("EOF")
			},
		}
		Parse(decoder, entries, decoderErrors)
		assert.EqualError(t, <-decoderErrors, decodeErr.Error())
	})
}

type mockDecoder struct {
	DecodeElementFn func(v interface{}, start *xml.StartElement) error
	TokenFn         func() (xml.Token, error)
}

func (d *mockDecoder) DecodeElement(v interface{}, start *xml.StartElement) error {
	return d.DecodeElementFn(v, start)
}

func (d *mockDecoder) Token() (xml.Token, error) {
	return d.TokenFn()
}
