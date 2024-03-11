package genbank

import "fmt"

// A GenbankSyntaxError denotes a sytntax error in
// a Genbank flatfile.
type GenbankSyntaxError struct {
	Line     uint
	Context  string
	Msg      string
	InnerErr error
}

// Error returns a human-readable error message.
func (gse GenbankSyntaxError) Error() string {
	msg := gse.Msg
	if gse.InnerErr != nil {
		msg = fmt.Errorf("%v: %w", msg, gse.InnerErr).Error()
	}

	return fmt.Sprintf("syntax error at line %v: %v\n%v\t%v", gse.Line, msg, gse.Line, gse.Context)
}

// Unwrap returns any errors underlying the syntax error, if applicable.
func (gse GenbankSyntaxError) Unwrap() error {
	return gse.InnerErr
}
