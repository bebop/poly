/*
Package bio provides utilities for reading and writing sequence data.
*/
package bio

import (
	"io"
)

/*
This package is supposed to be empty and only exists to provide a doc string.
Otherwise its namespace would collide with Go's native IO package.
*/

type Parser[T any, TH any] interface {
	Header() (TH, error)
	Next() (T, error)
	MaxLineCount() int64
}

type Writer interface {
	Write(io.Writer) error
}
