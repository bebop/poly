package cif

import "fmt"

// A CIFSyntaxError is a syntax error produced while parsing a CIF file.
type CIFSyntaxError struct {
	Line int
	Msg  string
}

// Wrap creates a new CIFSyntaxError wrapped with msg.
func (s CIFSyntaxError) Wrap(format string, a ...any) error {
	return CIFSyntaxError{
		Line: s.Line,
		Msg:  fmt.Sprintf("%s: %s", fmt.Sprintf(format, a...), s.Msg),
	}
}

// Error returns the formatted error message.
func (s CIFSyntaxError) Error() string {
	return fmt.Sprintf("CIF syntax error at line %v: %s", s.Line, s.Msg)
}
