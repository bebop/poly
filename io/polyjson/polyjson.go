/*
Package polyjson provides utilities to read and write poly.Sequence structs as JSON.

Poly's JSON schema is still in flux so be on the lookout for breaking changes as we
approach the 1.0 release.
*/
package polyjson

import (
	"encoding/json"
	"io/ioutil"
	"net/url"
	"time"

	"github.com/TimothyStiles/poly/io/poly"
)

/******************************************************************************

JSON specific IO related things begin here.

******************************************************************************/

type Poly struct {
	Meta     Meta      `json:"meta"`
	Features []Feature `json:"features"`
	Sequence string    `json:"sequence"`
}

type Meta struct {
	Name        string    `json:"name"`
	Hash        string    `json:"hash"`
	Description string    `json:"description"`
	URL         url.URL   `json:"url"`
	CreatedBy   string    `json:"created_by"`
	CreatedWith string    `json:"created_with"`
	CreatedOn   time.Time `json:"created_on"`
	Schema      string    `json:"schema"`
}
type Feature struct {
	Name           string            `json:"name"`
	Hash           string            `json:"hash"`
	Type           string            `json:"type"`
	Description    string            `json:"description"`
	Location       poly.Location     `json:"location"`
	Tags           map[string]string `json:"tags"`
	Sequence       string            `json:"sequence"`
	ParentSequence *Poly             `json:"-"`
}

func (sequence *Poly) AddFeature(feature *Feature) {
	feature.ParentSequence = sequence
	sequence.Features = append(sequence.Features, *feature)
}

func (sequence *Poly) GetFeatures() ([]Feature, error) {
	return sequence.Features, nil
}

func (sequence *Poly) GetSequence() (string, error) {
	return sequence.Sequence, nil
}

func (sequence *Poly) GetMeta() (Meta, error) {
	return sequence.Meta, nil
}

func (sequence *Poly) Write(path string) {
	Write(*sequence, path)
}

// Parse parses a Poly JSON file and adds appropriate pointers to struct.
func Parse(file []byte) Poly {
	var sequence Poly
	_ = json.Unmarshal([]byte(file), &sequence)
	legacyFeatures := sequence.Features
	sequence.Features = []Feature{}

	for _, feature := range legacyFeatures {
		sequence.AddFeature(&feature)
	}
	return sequence
}

// Read reads a Poly JSON file.
func Read(path string) Poly {
	file, _ := ioutil.ReadFile(path)
	sequence := Parse(file)
	return sequence
}

// Write writes a Poly struct out to json.
func Write(sequence Poly, path string) {
	file, _ := json.MarshalIndent(sequence, "", " ")
	_ = ioutil.WriteFile(path, file, 0644)
}

/******************************************************************************

JSON specific IO related things end here.

******************************************************************************/
