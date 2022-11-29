package fasta

import (
	"compress/gzip"
	"errors"
	"io"
	"os"
	"reflect"
	"testing"

	"github.com/stretchr/testify/assert"
)

func BenchmarkFastaLegacy(b *testing.B) {
	var fastas []Fasta
	for i := 0; i < b.N; i++ {
		fastas, _ = Read("data/uniprot_1mb_test.fasta")
	}
	_ = fastas
}

func BenchmarkParser(b *testing.B) {
	var fastas []Fasta
	for i := 0; i < b.N; i++ {
		fp, _ := os.Open("data/uniprot_1mb_test.fasta")
		p := NewParser(fp)
		for {
			fasta, _, err := p.ParseNext()
			if err != nil {
				break
			}
			fastas = append(fastas, fasta)
		}
		fastas = nil // Reset memory
	}
	_ = fastas
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
