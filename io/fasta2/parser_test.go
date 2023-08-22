package fasta2_test

import (
	"compress/gzip"
	_ "embed"
	"fmt"
	"os"
	"strings"
	"testing"

	"github.com/TimothyStiles/poly/io/fasta2"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

func TestParserNext(t *testing.T) {
	tests := []struct {
		name    string
		text    string
		err     error
		want    fasta2.Record
		lines   int
		hasNext bool
		wantErr bool
	}{
		{
			name: "OneSequenceNoNewlineTrmination",
			text: ">Header\nSEQLINE\nSEQLINE\n*",
			want: fasta2.Record{
				Header:   "Header",
				Sequence: "SEQLINESEQLINE*",
			},
			hasNext: false,
			lines:   4,
		},
		{
			name: "OneSequenceWithNewlineTrmination",
			text: ">Header\nSEQLINE\nSEQLINE\n*\n",
			want: fasta2.Record{
				Header:   "Header",
				Sequence: "SEQLINESEQLINE*",
			},
			hasNext: false,
			lines:   4,
		},
		{
			name: "OneSequencePrePostNewLines",
			text: "\n\n>Header\nSEQLINE\nSEQLINE\n*\n\n\n",
			err:  nil,
			want: fasta2.Record{
				Header:   "Header",
				Sequence: "SEQLINESEQLINE*",
			},
			hasNext: false,
			lines:   8,
		},
		{
			name: "TwoSequencesPrePostNewLines",
			text: "\n\n>Header\nSEQLINE\nSEQLINE\n*\n\n\n>Header\nSEQLINE\nSEQLINE\n*\n",
			err:  nil,
			want: fasta2.Record{
				Header:   "Header",
				Sequence: "SEQLINESEQLINE*",
			},
			hasNext: true,
			lines:   9,
		},
		{
			name:    "Empty",
			text:    "",
			want:    fasta2.Record{},
			hasNext: false,
			lines:   0,
		},
		{
			name:    "JustCommentsAndNewlines",
			text:    "\n\n;\n",
			want:    fasta2.Record{},
			hasNext: false,
			lines:   3,
		},
		{
			name: "InvalidHeader",
			text: ">Head",
			want: fasta2.Record{
				Header: "Head",
			},
			hasNext: false,
			lines:   1,
		},
		{
			name:    "InvalidNoHeadr",
			text:    "SEQLINE\nSEQLINE\n*\n",
			want:    fasta2.Record{},
			wantErr: true,
			err:     fmt.Errorf("invalid input: missing sequence header for sequence starting at line 1"),
			hasNext: true,
			lines:   1,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			buff := strings.NewReader(tt.text)
			p := fasta2.NewParser(buff)
			got, err := p.Next()
			if tt.wantErr {
				assert.EqualError(t, err, tt.err.Error())
			}
			assert.Equal(t, tt.want, got, "sequence")
			// this checks the state of the parser after a call to next
			assert.Equal(t, tt.hasNext, p.HasNext(), "HasNext")
			assert.Equal(t, tt.lines, p.Lines(), "Lines")
		})
	}
}

//go:embed testdata/sequences.fasta
var input string

func TestParser(t *testing.T) {
	expected := []fasta2.Record{
		{
			Header:   "Header 1",
			Sequence: "SEQLINESEQLINE*",
		},
		{
			Header:   "Header 2",
			Sequence: "SEQLINESEQLINE*",
		},
		{
			Header:   "Header 3",
			Sequence: "SEQLINESEQLINE*",
		},
	}
	t.Run("Success", func(t *testing.T) {
		var result []fasta2.Record
		par := fasta2.NewParser(strings.NewReader(input))
		for par.HasNext() {
			fas, err := par.Next()
			require.NoError(t, err)
			result = append(result, fas)
		}
		require.ElementsMatch(t, expected, result)
	})

	t.Run("SuccessGzip", func(t *testing.T) {
		f, err := os.Open("testdata/uniprot_1mb_test.fasta.gz")
		require.NoError(t, err)
		defer f.Close()

		gr, err := gzip.NewReader(f)
		require.NoError(t, err)
		defer gr.Close()

		recs, err := fasta2.ParseAll(gr)
		require.NoError(t, err)
		require.Len(t, recs, 12205)
	})
}
