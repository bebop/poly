package fasta2_test

import (
	"bytes"
	"io"
	"os"
	"path"
	"reflect"
	"testing"

	"github.com/TimothyStiles/poly/io/fasta2"
	"github.com/stretchr/testify/require"
)

func TestFastaString(t *testing.T) {
	type fields struct {
		Header   string
		Sequence string
	}
	tests := []struct {
		header string
		fields fields
		want   string
	}{
		{
			header: "success",
			fields: fields{
				Header:   "Cool Sequence",
				Sequence: "MDSKGSSQKGSRLLLLLVVSNLLLCQGVVSTPVCPNGPGNCQVSLRDLFDRAVMVSHYIHDLSSEMFNEFDKRYAQGKGFITMALNSCHTSSLPTPEDKEQAQQTHHEVLMSLILGLLRSWNDPLYHL",
			},
			want: ">Cool Sequence\nMDSKGSSQKGSRLLLLLVVSNLLLCQGVVSTPVCPNGPGNCQVSLRDLFDRAVMVSHYIHDLSSEMFNEFDKRYAQGKGF\nITMALNSCHTSSLPTPEDKEQAQQTHHEVLMSLILGLLRSWNDPLYHL\n",
		},
	}
	for _, tt := range tests {
		t.Run(tt.header, func(t *testing.T) {
			f := fasta2.Record{
				Header:   tt.fields.Header,
				Sequence: tt.fields.Sequence,
			}
			if got := f.String(); got != tt.want {
				t.Errorf("Record.String() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestRecord_Bytes(t *testing.T) {
	type fields struct {
		Header   string
		Sequence string
	}
	tests := []struct {
		name   string
		fields fields
		want   []byte
	}{
		{
			name: "success",
			fields: fields{
				Header:   "Cool Sequence",
				Sequence: "MDSKGSSQKGSRLLLLLVVSNLLLCQGVVSTPVCPNGPGNCQVSLRDLFDRAVMVSHYIHDLSSEMFNEFDKRYAQGKGFITMALNSCHTSSLPTPEDKEQAQQTHHEVLMSLILGLLRSWNDPLYHL",
			},
			want: []byte(">Cool Sequence\nMDSKGSSQKGSRLLLLLVVSNLLLCQGVVSTPVCPNGPGNCQVSLRDLFDRAVMVSHYIHDLSSEMFNEFDKRYAQGKGF\nITMALNSCHTSSLPTPEDKEQAQQTHHEVLMSLILGLLRSWNDPLYHL\n"),
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			r := fasta2.Record{
				Header:   tt.fields.Header,
				Sequence: tt.fields.Sequence,
			}
			if got := r.Bytes(); !reflect.DeepEqual(got, tt.want) {
				t.Errorf("Record.Bytes() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestRecordWrite(t *testing.T) {
	sequence := fasta2.Record{
		Header:   "Cool Sequence",
		Sequence: "MDSKGSSQKGSRLLLLLVVSNLLLCQGVVSTPVCPNGPGNCQVSLRDLFDRAVMVSHYIHDLSSEMFNEFDKRYAQGKGFITMALNSCHTSSLPTPEDKEQAQQTHHEVLMSLILGLLRSWNDPLYHL",
	}
	expected := ">Cool Sequence\nMDSKGSSQKGSRLLLLLVVSNLLLCQGVVSTPVCPNGPGNCQVSLRDLFDRAVMVSHYIHDLSSEMFNEFDKRYAQGKGF\nITMALNSCHTSSLPTPEDKEQAQQTHHEVLMSLILGLLRSWNDPLYHL\n"
	t.Run("success", func(t *testing.T) {
		r := sequence
		w := &bytes.Buffer{}
		err := r.Write(w)
		require.NoError(t, err)
		require.Equal(t, expected, w.String())
	})
	t.Run("fail truncated", func(t *testing.T) {
		r := sequence
		w := errorWriter{}
		err := r.Write(w)
		require.Error(t, err)
	})
}

func TestWrite(t *testing.T) {
	recs := []fasta2.Record{
		{
			Header:   "name1",
			Sequence: "seq1",
		},
		{
			Header:   "name2",
			Sequence: "seq2",
		},
	}
	t.Run("success", func(t *testing.T) {
		w := &bytes.Buffer{}
		err := fasta2.Write(recs, w)
		require.NoError(t, err)
		require.Equal(t, ">name1\nseq1\n>name2\nseq2\n", w.String())
	})
	t.Run("fail EOF", func(t *testing.T) {
		w := errorWriter{}
		err := fasta2.Write(recs, w)
		require.Error(t, err)
	})
}

func TestWriteFile(t *testing.T) {
	path := path.Join(os.TempDir(), "fasta_test")
	defer os.Remove(path) // clean up
	recs := []fasta2.Record{
		{
			Header:   "name1",
			Sequence: "seq1",
		},
		{
			Header:   "name2",
			Sequence: "seq2",
		},
	}
	err := fasta2.WriteFile(recs, path)
	require.NoError(t, err)
}

// errorWriter is a test utility to have errors
type errorWriter struct{}

func (ew errorWriter) Write(p []byte) (n int, err error) {
	return len(p), io.EOF
}
