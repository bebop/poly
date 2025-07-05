package genbank

import (
	"fmt"
	"strings"
)

const (
	dbNameIdx = 21
	dateIdx   = 26
	titleIdx  = 24
)

// String converts a Genbank to its Genbank flatfile representation.
func (g Genbank) String() string {
	builder := strings.Builder{}

	builder.WriteString(g.Header.String())

	for _, entry := range g.Entries {
		builder.WriteString(entry.String())
	}

	return builder.String()
}

func (e Entry) String() string {
	builder := strings.Builder{}

	// Write LOCUS entry.
	builder.WriteString(fmt.Sprintf("%-12s%-17s", "LOCUS", e.Name))
	if len(e.Name) >= 17 { // Ensure space after long names.
		builder.WriteRune(' ')
	}
	builder.WriteString(fmt.Sprintf(
		"%11d bp %3s%-4s %8s %3s %s",
		e.Length,
		e.Strandedness,
		e.MoleculeType,
		e.MoleculeToplogy,
		e.DivisionCode,
		strings.ToUpper(e.UpdateDate.Format("02-Jan-2006")),
	))

	builder.WriteString(fmt.Sprintf("%-12s%s", "DEFINITION", e.Definition))
	builder.WriteString(fmt.Sprintf("%-12s%s", "ACCESSION", e.Accession))
	builder.WriteString(fmt.Sprintf("%-12s%s.%d", "VERSION", e.Accession, e.AccessionVersion))

	wroteDBLinkKeyword := false
	for refType, refs := range e.DatabaseLinks {
		if !wroteDBLinkKeyword {
			builder.WriteString(fmt.Sprintf("%-12s%s:%s", "DBLINK", refType, strings.Join(refs, ",")))
			wroteDBLinkKeyword = true
		} else {
			builder.WriteString(fmt.Sprintf("%-12s%s:%s", "", refType, strings.Join(refs, ",")))
		}
	}

	builder.WriteString(fmt.Sprintf("%-12s", "KEYWORDS"))
	width := 0
	for i, keyword := range e.Keywords {
		width += len(keyword)
		if width > 80 { // If we exceed line width, newline before the keyword.
			builder.WriteString(fmt.Sprintf(";\n%-12s%s", "", keyword))
			width = 0
		} else if i > 0 {
			builder.WriteString(fmt.Sprintf("; %s", keyword))
		} else {
			builder.WriteString(keyword)
		}
	}
	builder.WriteRune('.')

	return builder.String()
}

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
