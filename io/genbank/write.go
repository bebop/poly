package genbank

import (
	"fmt"
	"strings"
)

const (
	dbNameIdx = 22
	dateIdx   = 27
	titleIdx  = 25
)

// String converts a Header to its Genbank flatfile representation.
func (h Header) String() string {
	builder := strings.Builder{}

	// Line 1: File name and database name
	builder.WriteString(h.FileName)
	builder.WriteString(strings.Repeat(" ", dbNameIdx-len(h.FileName)))
	builder.WriteString("Genetic Sequence Data Bank")
	builder.WriteRune('\n')

	// Line 2: Date
	builder.WriteString(strings.Repeat(" ", dateIdx))
	builder.WriteString(h.Date.Format(headerDateLayout))
	builder.WriteRune('\n')

	// Line 3: Blank
	builder.WriteRune('\n')

	// Line 4: Genbank release number
	builder.WriteString(fmt.Sprintf("                NCBI-GenBank Flat File Release %v.%v", h.MajorRelease, h.MinorRelease))
	builder.WriteRune('\n')

	// Line 5: Blank
	builder.WriteRune('\n')

	// Line 6: File title
	builder.WriteString(strings.Repeat(" ", titleIdx))
	builder.WriteString(h.Title)
	builder.WriteRune('\n')

	// Line 7: Blank
	builder.WriteRune('\n')

	// Line 8: File statistics
	builder.WriteString(fmt.Sprintf(
		"%8v loci, %11v bases, from %8v reported sequences",
		h.NumEntries,
		h.NumBases,
		h.NumSequences))
	builder.WriteRune('\n')

	return builder.String()
}
