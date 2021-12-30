/*
Package polyjson provides utilities to read and write poly.Sequence structs as JSON.

Poly's JSON schema is still in flux so be on the lookout for breaking changes as we
approach the 1.0 release.
*/
package polyjson

import (
	"encoding/json"
	"io/ioutil"
	"time"
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
	URL         string    `json:"url"`
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
	Location       Location          `json:"location"`
	Tags           map[string]string `json:"tags"`
	Sequence       string            `json:"sequence"`
	ParentSequence *Poly             `json:"-"`
}

type Location struct {
	Start             int        `json:"start"`
	End               int        `json:"end"`
	Complement        bool       `json:"complement"`
	Join              bool       `json:"join"`
	FivePrimePartial  bool       `json:"five_prime_partial"`
	ThreePrimePartial bool       `json:"three_prime_partial"`
	SubLocations      []Location `json:"sub_locations"`
}

func (sequence *Poly) AddFeature(feature *Feature) error {
	feature.ParentSequence = sequence
	sequence.Features = append(sequence.Features, *feature)
	return nil
}

// Parse parses a Poly JSON file and adds appropriate pointers to struct.
func Parse(file []byte) (Poly, error) {
	var sequence Poly
	err := json.Unmarshal([]byte(file), &sequence)
	if err != nil {
		return sequence, err
	}
	legacyFeatures := sequence.Features
	sequence.Features = []Feature{}

	for _, feature := range legacyFeatures {
		sequence.AddFeature(&feature)
	}
	return sequence, nil
}

// Read reads a Poly JSON file.
func Read(path string) (Poly, error) {
	file, err := ioutil.ReadFile(path)
	if err != nil {
		return Poly{}, err
	}
	sequence, err := Parse(file)
	if err != nil {
		return Poly{}, err
	}
	return sequence, nil
}

// Write writes a Poly struct out to json.
func Write(sequence Poly, path string) error {
	file, err := json.MarshalIndent(sequence, "", " ")
	if err != nil {
		return err
	}
	err = ioutil.WriteFile(path, file, 0644)
	if err != nil {
		return err
	}
	return nil
}

/******************************************************************************

JSON specific IO related things end here.

******************************************************************************/
