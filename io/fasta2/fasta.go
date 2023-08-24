package fasta2

import (
	"bytes"
	"fmt"
	"io"
	"os"
)

// Record is a struct representing a single Record element with a Name and its corresponding Sequence.
type Record struct {
	Header   string `json:"header"`
	Sequence string `json:"sequence"`
}

// buffer is a utility method to serialize the Record in a buffer.
// The buffer is reset before any writing happens
func (r Record) buffer(b *bytes.Buffer) {
	b.Reset()
	// grow the buffer to allocate just once, the numbers are in order:
	// the header + > + \n, the sequence + one \n for each 80 char, the last \n
	b.Grow(len(r.Header) + 2 + len(r.Sequence) + (len(r.Sequence) % 80) + 1)
	b.WriteByte('>')
	b.WriteString(r.Header)
	for i, c := range r.Sequence {
		// write the fasta sequence 80 characters at a time
		if i%80 == 0 {
			b.WriteByte('\n')
		}
		b.WriteRune(c)
	}
	b.WriteByte('\n')

	return
}

// returns the string representation of a Record.
func (r Record) String() string {
	var b bytes.Buffer
	r.buffer(&b)
	return b.String()
}

// returns the representation of a Record as []byte.
func (r Record) Bytes() []byte {
	var b bytes.Buffer
	r.buffer(&b)
	return b.Bytes()
}

// Writes the Record []byte representation to the passed io.Writer.
func (r Record) Write(w io.Writer) (int, error) {
	recBytes := r.Bytes()
	n, err := w.Write(recBytes)
	if err != nil {
		return n, fmt.Errorf("error writing record to io.Writer: %w", err)
	}

	return n, nil
}

// Write writes a fasta array to an io.Writer
func Write(recs []Record, w io.Writer) (int, error) {
	var (
		b       bytes.Buffer
		written int
	)
	for _, r := range recs {
		r.buffer(&b)
		n, err := w.Write(b.Bytes())
		written += n
		if err != nil {
			return written, err
		}
	}
	return written, nil
}

// WriteFile writes all the passed records to the file at path.
func WriteFile(recs []Record, path string) error {
	f, err := os.Create(path)
	if err != nil {
		return fmt.Errorf("error opening file %q: %w", path, err)
	}
	defer f.Close()
	for _, r := range recs {
		_, err := r.Write(f)
		if err != nil {
			return fmt.Errorf("error writing to file %q: %w", path, err)
		}
	}

	return nil
}
