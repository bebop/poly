/*
Package polyjson provides utilities to read and write poly.Sequence structs as JSON.

Poly's JSON schema is still in flux so be on the lookout for breaking changes as we
approach the 1.0 release.
*/

package polyjson

import (
	"encoding/json"
	"io/ioutil"

	"github.com/TimothyStiles/poly/io/poly"
)

/******************************************************************************

JSON specific IO related things begin here.

******************************************************************************/

// Parse parses a poly.Sequence JSON file and adds appropriate pointers to struct.
func Parse(file []byte) poly.Sequence {
	var sequence poly.Sequence
	_ = json.Unmarshal([]byte(file), &sequence)
	legacyFeatures := sequence.Features
	sequence.Features = []poly.Feature{}

	for _, feature := range legacyFeatures {
		sequence.AddFeature(&feature)
	}
	return sequence
}

// Read reads a poly.Sequence JSON file.
func Read(path string) poly.Sequence {
	file, _ := ioutil.ReadFile(path)
	sequence := Parse(file)
	return sequence
}

// Write writes a poly.Sequence struct out to json.
func Write(sequence poly.Sequence, path string) {
	file, _ := json.MarshalIndent(sequence, "", " ")
	_ = ioutil.WriteFile(path, file, 0644)
}

/******************************************************************************

JSON specific IO related things end here.

******************************************************************************/
