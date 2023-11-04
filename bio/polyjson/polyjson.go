/*
Package polyjson provides utilities to read and write poly.Sequence structs as JSON.

Poly's JSON schema is still in flux so be on the lookout for breaking changes as we
approach the 1.0 release.
*/
package polyjson

import (
	"bytes"
	"encoding/json"
	"io"
	"os"
	"time"

	"github.com/TimothyStiles/poly/bio/genbank"
	"github.com/TimothyStiles/poly/transform"
)

/******************************************************************************

JSON specific IO related things begin here.

******************************************************************************/

var (
	marshalIndentFn = json.MarshalIndent
	readFileFn      = os.Open
	unmarshalFn     = json.Unmarshal
)

// Poly is poly's native JSON representation of a sequence.
type Poly struct {
	Meta     Meta      `json:"meta"`
	Features []Feature `json:"features"`
	Sequence string    `json:"sequence"`
}

// Meta contains all the metadata for a poly sequence struct.
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

// Feature contains all the feature data for a poly feature struct.
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

// Location contains all the location data for a poly feature's location.
type Location struct {
	Start             int        `json:"start"`
	End               int        `json:"end"`
	Complement        bool       `json:"complement"`
	Join              bool       `json:"join"`
	FivePrimePartial  bool       `json:"five_prime_partial"`
	ThreePrimePartial bool       `json:"three_prime_partial"`
	SubLocations      []Location `json:"sub_locations"`
}

// AddFeature adds a feature to a Poly struct. Does not add the feature's sequence
func (sequence *Poly) AddFeature(feature *Feature) error {
	feature.ParentSequence = sequence
	sequence.Features = append(sequence.Features, *feature)
	return nil
}

// GetSequence takes a feature and returns a sequence string for that feature.
func (feature Feature) GetSequence() (string, error) {
	return getFeatureSequence(feature, feature.Location)
}

// getFeatureSequence takes a feature and location object and returns a sequence string.
func getFeatureSequence(feature Feature, location Location) (string, error) {
	var sequenceBuffer bytes.Buffer
	var sequenceString string
	parentSequence := feature.ParentSequence.Sequence

	if len(location.SubLocations) == 0 {
		sequenceBuffer.WriteString(parentSequence[location.Start:location.End])
	} else {
		for _, subLocation := range location.SubLocations {
			sequence, err := getFeatureSequence(feature, subLocation)
			if err != nil { // todo: test error
				return "", err
			}
			sequenceBuffer.WriteString(sequence)
		}
	}

	// reverse complements resulting string if needed.
	if location.Complement {
		sequenceString = transform.ReverseComplement(sequenceBuffer.String())
	} else {
		sequenceString = sequenceBuffer.String()
	}

	return sequenceString, nil
}

// Parse parses a Poly JSON file and adds appropriate pointers to struct.
func Parse(file io.Reader) (Poly, error) {
	var sequence Poly
	buf := new(bytes.Buffer)
	_, err := buf.ReadFrom(file) // todo: test error
	if err != nil {
		return sequence, err
	}

	if err := unmarshalFn(buf.Bytes(), &sequence); err != nil {
		return sequence, err
	}

	legacyFeatures := sequence.Features
	sequence.Features = []Feature{}

	for _, feature := range legacyFeatures {
		err = sequence.AddFeature(&feature)
		if err != nil {
			return sequence, err
		}
	}
	return sequence, nil
}

// Read reads a Poly JSON file.
func Read(path string) (Poly, error) {
	file, err := readFileFn(path)
	if err != nil {
		return Poly{}, err
	}
	return Parse(file)
}

// Write writes a Poly struct out to json.
func Write(sequence Poly, path string) error {
	file, err := marshalIndentFn(sequence, "", " ")
	if err != nil {
		return err
	}
	return os.WriteFile(path, file, 0644)
}

// Utilities to convert polyjson objects -> their genbank equivalents
// TODO add convert <- genbank methods, which is currently difficult as most
// genbank Meta values are discarded due to lack of support for wildcard metadata in polyjson.

func (sequence *Poly) ToGenbank() genbank.Genbank {
	gb := genbank.Genbank{
		Meta:     sequence.Meta.ToGenbank(),
		Features: make([]genbank.Feature, len(sequence.Features)),
		Sequence: sequence.Sequence,
	}
	for i, f := range sequence.Features {
		gb.Features[i] = f.ToGenbank()
		gb.Features[i].ParentSequence = &gb
	}
	return gb
}

func (meta *Meta) ToGenbank() genbank.Meta {
	other := make(map[string]string)
	if meta.URL != "" {
		other["URL"] = meta.URL
	}
	if meta.CreatedBy != "" {
		other["CreatedBy"] = meta.CreatedBy
	}
	if meta.CreatedWith != "" {
		other["CreatedWith"] = meta.CreatedWith
	}
	other["CreatedOn"] = meta.CreatedOn.String()
	if meta.Schema != "" {
		other["Schema"] = meta.Schema
	}
	return genbank.Meta{
		Definition:   meta.Description,
		Source:       meta.CreatedBy,
		Origin:       meta.CreatedWith,
		Name:         meta.Name,
		SequenceHash: meta.Hash,
		Other:        other,
	}
}

func (feature *Feature) ToGenbank() genbank.Feature {
	attributes := genbank.NewMultiMap[string, string]()
	for key, value := range feature.Tags {
		genbank.Put(attributes, key, value)
	}
	genbank.Put(attributes, "Name", feature.Name)

	return genbank.Feature{
		Type:         feature.Type,
		Description:  feature.Description,
		Attributes:   attributes,
		SequenceHash: feature.Hash,
		Sequence:     feature.Sequence,
		Location:     feature.Location.ToGenbank(),
	}
}

func (location *Location) ToGenbank() genbank.Location {
	loc := genbank.Location{
		Start:             location.Start,
		End:               location.End,
		Complement:        location.Complement,
		Join:              location.Join,
		FivePrimePartial:  location.FivePrimePartial,
		ThreePrimePartial: location.ThreePrimePartial,
		SubLocations: genbank.MapSlice(
			location.SubLocations,
			func(s Location) genbank.Location { return s.ToGenbank() },
		),
	}
	loc.GbkLocationString = genbank.BuildLocationString(loc)
	return loc
}

/******************************************************************************

JSON specific IO related things end here.

******************************************************************************/
