package uniprot

import (
	"compress/gzip"
	"encoding/xml"
	"io"
	"os"
)

// ReadUniprot reads a gzipped Uniprot XML dump.
func ReadUniprot(path string, entries chan Entry) {
	xmlFile, _ := os.Open(path)
	r, _ := gzip.NewReader(xmlFile)
	go ParseUniprot(r, entries)
}

// ParseUniprot parses Uniprot entries into a channel.
func ParseUniprot(r io.Reader, entries chan<- Entry) {
	decoder := xml.NewDecoder(r)
	for {
		t, _ := decoder.Token()
		if t == nil {
			break
		}
		switch se := t.(type) {
		case xml.StartElement:
			switch se.Name.Local {
			case "entry":
				var e Entry
				decoder.DecodeElement(&e, &se)
				entries <- e
			}
		}
	}
	close(entries)
}
