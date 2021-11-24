/*
Package poly provides the core Sequence struct and methods for working with DNA, RNA, and Amino Acid sequences in Poly.

Each Sequence struct can be broken down into the following:

Meta: Author, Date, Description, Name, Version, etc.
Features:	A list of Feature structs containing feature locations and other information.
Sequece:	The actual sequence string.

As well as other info like source file checksums, etc.

This package will likely be overhauled before 1.0 release to be more flexible and robust.
So be on the lookup for breaking changes in future releases.
*/
package poly

/******************************************************************************

File is structured as so:

Structs:
	Sequence - main struct for sequence handling plus sub structs.

******************************************************************************/

/******************************************************************************

Sequence related structs begin here.

******************************************************************************/

type Meta interface {
	GetMeta() (Meta, error)
}
type Sequence interface {
	GetSequence() (string, error)
}

type Feature interface {
	GetParent() (*AnnotatedSequence, error)
	GetSequence() (Sequence, error)
	GetLocation() (Location, error)
	SetLocation(Location) error
	GetAttributes() (map[string]string, error)
}

type AnnotatedSequence interface {
	GetMeta() (Meta, error)
	GetFeatures() ([]Feature, error)
	GetSequence() (string, error)
	// AddFeature(Feature) error
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

// Meta Holds all the meta information of an Sequence struct.
// type Meta struct {
// 	Name                 string   `json:"name"`
// 	SequenceHash         string   `json:"sequence_hash"`
// 	SequenceHashFunction string   `json:"hash_function"`
// 	CheckSum             [32]byte `json:"checkSum"` // blake3 checksum of the parsed file itself. Useful for if you want to check if incoming genbank/gff files are different.
// }

// Location holds nested location info for sequence region.

// Feature holds a single annotation in a struct. from https://github.com/blachlylab/gff3/blob/master/gff3.go
// type Feature struct {
// 	Name string //Seqid in gff, name in gbk
// 	//gff specific
// 	Source               string            `json:"source"`
// 	Type                 string            `json:"type"`
// 	Score                string            `json:"score"`
// 	Strand               string            `json:"strand"`
// 	Phase                string            `json:"phase"`
// 	Attributes           map[string]string `json:"attributes"` //<- used in both genbank and gff
// 	GbkLocationString    string            `json:"gbk_location_string"`
// 	Sequence             string            `json:"sequence"`
// 	SequenceLocation     Location          `json:"sequence_location"`
// 	SequenceHash         string            `json:"sequence_hash"`
// 	Description          string            `json:"description"`
// 	SequenceHashFunction string            `json:"hash_function"`
// 	ParentSequence       *Sequence         `json:"-"`
// }

// Sequence holds all sequence information in a single struct.
// type Sequence struct {
// 	Meta        Meta      `json:"meta"`
// 	Description string    `json:"description"`
// 	Sequence    string    `json:"sequence"`
// 	Features    []Feature `json:"features"`
// }

// AddFeature is the canonical way to add a Feature into a Sequence struct. Appending a Feature struct directly to Sequence.Feature's will break .GetSequence() method.
// func (sequence *Sequence) AddFeature(feature *Feature) []Feature {
// 	feature.ParentSequence = sequence
// 	var featureCopy Feature = *feature
// 	sequence.Features = append(sequence.Features, featureCopy)
// 	return sequence.Features
// }

// GetSequence is a method wrapper to get a Feature's sequence. Mutates with Sequence.
// func (feature Feature) GetSequence() string {
// 	return getFeatureSequence(feature, feature.SequenceLocation)
// }

// // getFeatureSequence takes a feature and location object and returns a sequence string.
// func getFeatureSequence(feature Feature, location Location) string {
// 	var sequenceBuffer bytes.Buffer
// 	var sequenceString string
// 	parentSequence := feature.ParentSequence.Sequence

// 	if len(location.SubLocations) == 0 {
// 		sequenceBuffer.WriteString(parentSequence[location.Start:location.End])
// 	} else {

// 		for _, subLocation := range location.SubLocations {
// 			sequenceBuffer.WriteString(getFeatureSequence(feature, subLocation))
// 		}
// 	}

// 	// reverse complements resulting string if needed.
// 	if location.Complement {
// 		sequenceString = transform.ReverseComplement(sequenceBuffer.String())
// 	} else {
// 		sequenceString = sequenceBuffer.String()
// 	}

// 	return sequenceString
// }

/******************************************************************************

Sequence related structs end here.

******************************************************************************/
