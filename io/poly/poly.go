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

type AnnotatedSequence interface {
	AddFeature
}
