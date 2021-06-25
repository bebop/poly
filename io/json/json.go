package json

import (
	"encoding/json"
	"io/ioutil"

	"github.com/TimothyStiles/poly"
)

/******************************************************************************

JSON specific IO related things begin here.

******************************************************************************/

// ParseJSON parses an poly.Sequence JSON file and adds appropriate pointers to struct.
func ParseJSON(file []byte) poly.Sequence {
	var sequence poly.Sequence
	_ = json.Unmarshal([]byte(file), &sequence)
	legacyFeatures := sequence.Features
	sequence.Features = []poly.Feature{}

	for _, feature := range legacyFeatures {
		sequence.AddFeature(feature)
	}
	return sequence
}

// ReadJSON reads an poly.Sequence JSON file.
func ReadJSON(path string) poly.Sequence {
	file, _ := ioutil.ReadFile(path)
	sequence := ParseJSON(file)
	return sequence
}

// WriteJSON writes an poly.Sequence struct out to json.
func WriteJSON(sequence poly.Sequence, path string) {
	file, _ := json.MarshalIndent(sequence, "", " ")
	_ = ioutil.WriteFile(path, file, 0644)
}

/******************************************************************************

JSON specific IO related things end here.

******************************************************************************/
