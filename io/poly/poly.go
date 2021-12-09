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

Sequence related interfaces begin here.

******************************************************************************/

type Meta interface {
	GetMeta() (Meta, error)
}
type Sequence interface {
	GetSequence() (string, error)
}

type Feature interface {
	// GetParent() (*AnnotatedSequence, error)
	GetSequence() (string, error)
	// GetLocation() (Location, error)
	// SetLocation(Location) error
	GetType() (string, error)
	// GetAttributes() (map[string]string, error)
}

type AnnotatedSequence interface {
	GetMeta() (Meta, error)
	GetFeatures() ([]Feature, error)
	GetSequence() (string, error)
	// AddFeature(*Feature) error
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

/******************************************************************************

Sequence related structs end here.

******************************************************************************/
