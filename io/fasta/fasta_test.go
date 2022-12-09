package fasta

import (
	"compress/gzip"
	"errors"
	"io"
	"os"
	"reflect"
	"strings"
	"testing"

	"github.com/stretchr/testify/assert"
)

// Initialized at TestMain.
var uniprotFasta string

func TestMain(m *testing.M) {
	const uniprotFastaGzFilePath = "data/uniprot_1mb_test.fasta.gz"
	// unzip uniprot data and create uniprotFasta string for benchmarks and testing.
	uniprotFastaGzFile, err := os.Open(uniprotFastaGzFilePath)
	if err != nil {
		panic(uniprotFastaGzFilePath + " required for tests!")
	}
	defer uniprotFastaGzFile.Close()

	uniprotFastaGzReader, _ := gzip.NewReader(uniprotFastaGzFile)
	defer uniprotFastaGzReader.Close()

	uniprotFastaBytes, err := io.ReadAll(uniprotFastaGzReader)
	if err != nil {
		panic(err)
	}

	uniprotFasta = string(uniprotFastaBytes)
	m.Run()
}

func BenchmarkFastaLegacy(b *testing.B) {
	var fastas []Fasta
	var err error
	for i := 0; i < b.N; i++ {
		fastas, err = Parse(strings.NewReader(uniprotFasta))
		if err != nil {
			b.Fatal(err)
		}
	}
	_ = fastas
}

func BenchmarkParser(b *testing.B) {
	var fastaRecords []Fasta
	for i := 0; i < b.N; i++ {
		parser := NewParser(strings.NewReader(uniprotFasta), 256)
		for {
			fasta, _, err := parser.ParseNext()
			if err != nil {
				if !errors.Is(err, io.EOF) {
					b.Fatal(err)
				}
				break
			}
			fastaRecords = append(fastaRecords, fasta)
		}
		fastaRecords = nil // Reset memory
	}
	_ = fastaRecords
}

func TestParse(t *testing.T) {
	type args struct {
		r io.Reader
	}
	tests := []struct {
		name    string
		args    args
		want    []Fasta
		wantErr bool
	}{
		// TODO: Add test cases.
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := Parse(tt.args.r)
			if (err != nil) != tt.wantErr {
				t.Errorf("Parse() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if !reflect.DeepEqual(got, tt.want) {
				t.Errorf("Parse() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestParseConcurrent(t *testing.T) {
	type args struct {
		r         io.Reader
		sequences chan<- Fasta
	}
	tests := []struct {
		name string
		args args
	}{
		// TODO: Add test cases.
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			ParseConcurrent(tt.args.r, tt.args.sequences)
		})
	}
}

func TestReadGzConcurrent(t *testing.T) {
	type args struct {
		path      string
		sequences chan<- Fasta
	}
	tests := []struct {
		name string
		args args
	}{
		// TODO: Add test cases.
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			ReadGzConcurrent(tt.args.path, tt.args.sequences)
		})
	}
}

func TestReadConcurrent(t *testing.T) {
	type args struct {
		path      string
		sequences chan<- Fasta
	}
	tests := []struct {
		name string
		args args
	}{
		// TODO: Add test cases.
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			ReadConcurrent(tt.args.path, tt.args.sequences)
		})
	}
}

func TestReadGz(t *testing.T) {
	type args struct {
		path string
	}
	tests := []struct {
		name    string
		args    args
		want    []Fasta
		wantErr bool
	}{
		// TODO: Add test cases.
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := ReadGz(tt.args.path)
			if (err != nil) != tt.wantErr {
				t.Errorf("ReadGz() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if !reflect.DeepEqual(got, tt.want) {
				t.Errorf("ReadGz() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestRead(t *testing.T) {
	type args struct {
		path string
	}
	tests := []struct {
		name    string
		args    args
		want    []Fasta
		wantErr bool
	}{
		// TODO: Add test cases.
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := Read(tt.args.path)
			if (err != nil) != tt.wantErr {
				t.Errorf("Read() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if !reflect.DeepEqual(got, tt.want) {
				t.Errorf("Read() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestBuild(t *testing.T) {
	type args struct {
		fastas []Fasta
	}
	tests := []struct {
		name    string
		args    args
		want    []byte
		wantErr bool
	}{
		// TODO: Add test cases.
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := Build(tt.args.fastas)
			if (err != nil) != tt.wantErr {
				t.Errorf("Build() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if !reflect.DeepEqual(got, tt.want) {
				t.Errorf("Build() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestWrite(t *testing.T) {
	type args struct {
		fastas []Fasta
		path   string
	}
	tests := []struct {
		name    string
		args    args
		wantErr bool
	}{
		// TODO: Add test cases.
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if err := Write(tt.args.fastas, tt.args.path); (err != nil) != tt.wantErr {
				t.Errorf("Write() error = %v, wantErr %v", err, tt.wantErr)
			}
		})
	}
}

func TestRead_error(t *testing.T) {
	t.Run("Read errors opening the file", func(t *testing.T) {
		openErr := errors.New("open error")
		oldOpenFn := openFn
		openFn = func(name string) (*os.File, error) {
			return nil, openErr
		}
		defer func() {
			openFn = oldOpenFn
		}()
		_, err := Read("/tmp/file")
		assert.EqualError(t, err, openErr.Error())
	})

	t.Run("ReadGz errors opening the file", func(t *testing.T) {
		openErr := errors.New("open error")
		oldOpenFn := openFn
		openFn = func(name string) (*os.File, error) {
			return nil, openErr
		}
		defer func() {
			openFn = oldOpenFn
		}()
		_, err := ReadGz("/tmp/file")
		assert.EqualError(t, err, openErr.Error())
	})

	t.Run("ReadGz errors reading the file", func(t *testing.T) {
		readErr := errors.New("read error")
		oldOpenFn := openFn
		oldGzipReaderFn := gzipReaderFn
		openFn = func(name string) (*os.File, error) {
			return &os.File{}, nil
		}
		gzipReaderFn = func(r io.Reader) (*gzip.Reader, error) {
			return nil, readErr
		}
		defer func() {
			openFn = oldOpenFn
			gzipReaderFn = oldGzipReaderFn
		}()
		_, err := ReadGz("/tmp/file")
		assert.EqualError(t, err, readErr.Error())
	})
}

func TestWrite_error(t *testing.T) {
	buildErr := errors.New("build error")
	oldBuildFn := buildFn
	buildFn = func(fastas []Fasta) ([]byte, error) {
		return nil, buildErr
	}
	defer func() {
		buildFn = oldBuildFn
	}()
	err := Write([]Fasta{}, "/tmp/file")
	assert.EqualError(t, err, buildErr.Error())
}

func TestParser(t *testing.T) {
	parser := NewParser(nil, 256)
	for testIndex, test := range []struct {
		content  string
		expected []Fasta
	}{
		{
			content:  ">humen\nGATTACA\nCATGAT", // EOF-ended Fasta not valid
			expected: []Fasta{},
		},
		{
			content:  ">humen\nGATTACA\nCATGAT\n",
			expected: []Fasta{{Name: "humen", Sequence: "GATTACACATGAT"}},
		},
		{
			content: ">doggy or something\nGATTACA\n\nCATGAT\n" +
				">homunculus\nAAAA\n",
			expected: []Fasta{
				{Name: "doggy or something", Sequence: "GATTACACATGAT"},
				{Name: "homunculus", Sequence: "AAAA"},
			},
		},
	} {
		parser.Reset(strings.NewReader(test.content))
		fastas, err := parser.ParseAll()
		if err != nil {
			t.Fatal(err)
		}
		if len(fastas) != len(test.expected) {
			t.Errorf("case index %d: got %d fastas, expected %d", testIndex, len(fastas), len(test.expected))
			continue
		}
		for index, gotFasta := range fastas {
			expected := test.expected[index]
			if expected != gotFasta {
				t.Errorf("got!=expected: %+v != %+v", gotFasta, expected)
			}
		}
	}
}

func TestParseBytes(t *testing.T) {
	// Partial read test.
	const testFasta = ">0\nGAT\n>1\nCAC\n"
	p := NewParser(strings.NewReader(testFasta), 256)
	result1, bytesRead, err := p.ParseByteLimited(1)
	if err != nil {
		t.Fatal(err)
	}
	if len(result1) != 1 {
		t.Error("expected result of length 1 (partial read)")
	}
	expectBytesRead := 1 + strings.Index(testFasta[1:], ">")
	if int(bytesRead) != expectBytesRead {
		t.Errorf("expected %d bytes read, got %d bytes read", expectBytesRead, bytesRead)
	}

	// Full read test.
	p.Reset(strings.NewReader(testFasta))
	result1, bytesRead, err = p.ParseByteLimited(100)
	if err != nil {
		t.Fatal(err)
	}
	if len(result1) != 2 {
		t.Error("expected result of length 2 (full read)")
	}
	expectBytesRead = len(testFasta)
	if int(bytesRead) != expectBytesRead {
		t.Errorf("expected %d bytes read, got %d bytes read", expectBytesRead, bytesRead)
	}
}

// TestReadEmptyFasta tests that an empty fasta file is parsed correctly.
func TestReadEmptyFasta(t *testing.T) {
	fastas, err := Read("data/empty.fasta")
	if err != nil {
		t.Fatal(err)
	}
	if len(fastas) != 1 {
		t.Errorf("expected 1 fastas, got %d", len(fastas))
	}
}
