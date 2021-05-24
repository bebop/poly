package uniprot

import (
	"compress/gzip"
	"encoding/xml"
	"io"
	"os"
)

/******************************************************************************
May 18, 2021

Uniprot is comprehensive, high-quality and freely accessible resource of protein
sequence and functional information. It the best(1) protein database out there.

Uniprot database dumps are available as gzipped FASTA files or gzipped XML files.
The XML files have significantly more information than the FASTA files, and this
parser specifically works on the gzipped XML files from Uniprot.

Uniprot provides an XML schema of their data dumps(3), which is useful for
autogeneration of Golang structs. I used xsdgen(4) to automatically generate
xml.go from uniprot.xsd

Each protein in Uniprot is known as an "Entry" (as defined in xml.go).

The function ParseUniprot stream-reads Uniprot into an Entry channel, from which
you can use the entries however you want. ReadUniprot simplifies reading gzipped
files from a disk into an Entry channel, essentially just preparing the reader for
ParseUniprot.

Cheers,
Keoni

(1) Opinion of Keoni Gandall as of May 18, 2021
(2) https://www.uniprot.org/downloads
(3) https://www.uniprot.org/docs/uniprot.xsd

******************************************************************************/

// ReadUniprot reads a gzipped Uniprot XML dump.
func ReadUniprot(path string) (chan Entry, chan error, error) {
	entries := make(chan Entry)
	errors := make(chan error)
	xmlFile, err := os.Open(path)
	if err != nil {
		return entries, errors, err
	}
	r, err := gzip.NewReader(xmlFile)
	if err != nil {
		return entries, errors, err
	}
	go ParseUniprot(r, entries, errors)
	return entries, errors, nil
}

// ParseUniprot parses Uniprot entries into a channel.
func ParseUniprot(r io.Reader, entries chan<- Entry, errors chan<- error) {
	decoder := xml.NewDecoder(r)
	for {
		decoderToken, err := decoder.Token()
		if err != nil {
			errors <- err
		}
		if decoderToken == nil {
			break
		}
		switch startElement := decoderToken.(type) {
		case xml.StartElement:
			if startElement.Name.Local == "entry" {
				var e Entry
				decoder.DecodeElement(&e, &startElement)
				entries <- e
			}
		}
	}
	close(entries)
	close(errors)
}
