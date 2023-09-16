/*
Package slow5 contains slow5 parsers and writers.

Right now, only parsing slow5 files is supported. Support for writing and blow5
coming soon.

slow5 is a file format alternative to fast5, which is the file format outputted
by Oxford Nanopore sequencing devices. fast5 uses hdf5, which is a complex file
format that can only be read and written with a single software library built
in 1998. On the other hand, slow5 uses a .tsv file format, which is easy to
both parse and write.

slow5 files contain both general metadata about the sequencing run and raw
signal reads from the sequencing run. This raw signal can be used directly or
basecalled and used for alignment.

More information on slow5 can be found here: https://github.com/hasindu2008/slow5tools
*/
package slow5

import (
	"bufio"
	"fmt"
	"io"
	"sort"
	"strconv"
	"strings"

	"github.com/koeng101/svb"
)

/******************************************************************************
Oct 10, 2021

slow5 parser begins here. Specification below:
https://hasindu2008.github.io/slow5specs/slow5-v1.0.0.pdf

slow5 is able to combine multiple sequencing runs into a single file format,
but we read each sequencing run separately. Each sequencing run contains a
header with metadata and a list of reads. However, unlike many other file
formats, a slow5 file should almost never be read into a common struct, since
most runs are large, and will take a ton of memory. Instead, the default way
to parse slow5 files is to produce a list of headers and a channel of raw
reads. In order to connect the two (if needed), use ReadGroupID.

Nanopore changes the attributes found in the header quite often, so we store
most of these attributes in a map for future proofing. Even the binary file
format, blow5, does not have types for these attributes, and just stores them
as a long string.

Reads have 8 required columns, and a few auxiliary. These are typed, since they
are what will probably be used in real software.

Cheers mate,

Keoni

******************************************************************************/

// Header contains metadata about the sequencing run in general.
type Header struct {
	HeaderValues []HeaderValue
}

type HeaderValue struct {
	ReadGroupID        uint32
	Slow5Version       string
	Attributes         map[string]string
	EndReasonHeaderMap map[string]int
}

// Read contains metadata and raw signal strengths for a single nanopore read.
type Read struct {
	ReadID       string
	ReadGroupID  uint32
	Digitisation float64
	Offset       float64
	Range        float64
	SamplingRate float64
	LenRawSignal uint64
	RawSignal    []int16

	// Auxiliary fields
	ChannelNumber string
	MedianBefore  float64
	ReadNumber    int32
	StartMux      uint8
	StartTime     uint64
	EndReason     string // enum{unknown,partial,mux_change,unblock_mux_change,data_service_unblock_mux_change,signal_positive,signal_negative}

	EndReasonMap map[string]int // Used for writing
}

var knownEndReasons = map[string]bool{"unknown": true,
	"partial":                         true,
	"mux_change":                      true,
	"unblock_mux_change":              true,
	"data_service_unblock_mux_change": true,
	"signal_positive":                 true,
	"signal_negative":                 true,
}

// Parser is a flexible parser that provides ample
// control over reading slow5 sequences.
// It is initialized with NewParser.
type Parser struct {
	// reader keeps state of current reader.
	reader             bufio.Reader
	line               uint
	headerMap          map[int]string
	endReasonMap       map[int]string
	endReasonHeaderMap map[string]int
	header             Header
	hitEOF             bool
}

// Header returns the header
func (parser *Parser) Header() (*Header, error) {
	return &parser.header, nil
}

// NewParser parsers a slow5 file.
func NewParser(r io.Reader, maxLineSize int) (*Parser, error) {
	parser := &Parser{
		reader: *bufio.NewReaderSize(r, maxLineSize),
		line:   0,
	}
	var headers []HeaderValue
	var slow5Version string
	var numReadGroups uint32
	headerMap := make(map[int]string)
	endReasonMap := make(map[int]string)
	endReasonHeaderMap := make(map[string]int)

	for {
		lineBytes, err := parser.reader.ReadSlice('\n')
		if err != nil {
			return parser, err
		}
		line := strings.TrimSpace(string(lineBytes))
		parser.line++
		values := strings.Split(line, "\t")
		if len(values) < 2 {
			return parser, fmt.Errorf("Got following line without tabs: %s", line)
		}

		// First, we need to identify the number of read groups. This number will be the length of our
		// ReadGroups output, and we will need it for iteration through the rest of the header.
		if numReadGroups == 0 {
			switch values[0] {
			case "#slow5_version":
				slow5Version = values[1]
			case "#num_read_groups":
				numReadGroupsUint, err := strconv.ParseUint(values[1], 10, 32)
				if err != nil {
					return parser, err
				}
				numReadGroups = uint32(numReadGroupsUint)
				for id := uint32(0); id < numReadGroups; id++ {
					headers = append(headers, HeaderValue{Slow5Version: slow5Version, ReadGroupID: id, Attributes: make(map[string]string)})
				}
			}
			continue
		}
		// Terminate if we hit the beginning of the raw read headers
		// Get endReasonEnums. This is simply a string between enum{} that is used for the reasons that a read could have ended.
		if values[0] == "#char*" {
			for _, typeInfo := range values {
				if strings.Contains(typeInfo, "enum") {
					endReasonEnumsMinusPrefix := strings.TrimPrefix(typeInfo, "enum{")
					endReasonEnumsMinusSuffix := strings.TrimSuffix(endReasonEnumsMinusPrefix, "}")
					endReasons := strings.Split(endReasonEnumsMinusSuffix, ",")

					for endReasonIndex, endReason := range endReasons {
						if _, ok := knownEndReasons[endReason]; !ok {
							return parser, fmt.Errorf("unknown end reason '%s' found in end_reason enum. Please report", endReason)
						}
						endReasonMap[endReasonIndex] = endReason
						endReasonHeaderMap[endReason] = endReasonIndex
					}
				}
			}
			// Add endReasonMap to each output header. This helps the writer
			for headerIndex := range headers {
				headers[headerIndex].EndReasonHeaderMap = endReasonHeaderMap
			}
			continue
		}

		// Get the read headers and their identifiers. Though the primary read headers are in a defined order, the auxiliary headers are not.
		if values[0] == "#read_id" {
			headerMap[0] = "read_id"
			for headerNum := 1; headerNum < len(values); headerNum++ {
				headerMap[headerNum] = values[headerNum]
			}
			break
		}

		// Check to make sure we have the right amount of information for the num_read_groups
		if len(values) != int(numReadGroups+1) {
			return parser, fmt.Errorf("Improper amount of information for read groups. Needed %d, got %d, in line: %s", numReadGroups+1, len(values), line)
		}
		for id := 0; id < int(numReadGroups); id++ {
			headers[id].Attributes[values[0]] = values[id+1]
		}
		continue
	}
	parser.headerMap = headerMap
	parser.endReasonMap = endReasonMap
	parser.endReasonHeaderMap = endReasonHeaderMap
	parser.header = Header{HeaderValues: headers}
	return parser, nil
}

// Next parses the next read from a parser.
func (parser *Parser) Next() (*Read, error) {
	if parser.hitEOF {
		return &Read{}, io.EOF
	}
	lineBytes, err := parser.reader.ReadSlice('\n')
	if err != nil {
		if err == io.EOF {
			parser.hitEOF = true
			if strings.TrimSpace(string(lineBytes)) == "" {
				// If EOF line is empty, return. If not, continue!
				return nil, io.EOF
			}
		} else {
			return nil, err
		}
	}
	parser.line++
	line := strings.TrimSpace(string(lineBytes))

	values := strings.Split(line, "\t")
	// Reads have started.
	// Once we have the read headers, start to parse the actual reads
	var newRead Read
	newRead.EndReasonMap = parser.endReasonHeaderMap
	for valueIndex := 0; valueIndex < len(values); valueIndex++ {
		fieldValue := parser.headerMap[valueIndex]
		if values[valueIndex] == "." {
			continue
		}
		switch fieldValue {
		case "read_id":
			newRead.ReadID = values[valueIndex]
		case "read_group":
			readGroupID, err := strconv.ParseUint(values[valueIndex], 10, 32)
			if err != nil {
				return nil, fmt.Errorf("Failed convert read_group '%s' to uint on line %d. Got Error: %w", values[valueIndex], parser.line, err)
			}
			newRead.ReadGroupID = uint32(readGroupID)
		case "digitisation":
			digitisation, err := strconv.ParseFloat(values[valueIndex], 64)
			if err != nil {
				return nil, fmt.Errorf("Failed to convert digitisation '%s' to float on line %d. Got Error: %w", values[valueIndex], parser.line, err)
			}
			newRead.Digitisation = digitisation
		case "offset":
			offset, err := strconv.ParseFloat(values[valueIndex], 64)
			if err != nil {
				return nil, fmt.Errorf("Failed to convert offset '%s' to float on line %d. Got Error: %w", values[valueIndex], parser.line, err)
			}
			newRead.Offset = offset
		case "range":
			nanoporeRange, err := strconv.ParseFloat(values[valueIndex], 64)
			if err != nil {
				return nil, fmt.Errorf("Failed to convert range '%s' to float on line %d. Got Error: %w", values[valueIndex], parser.line, err)
			}
			newRead.Range = nanoporeRange
		case "sampling_rate":
			samplingRate, err := strconv.ParseFloat(values[valueIndex], 64)
			if err != nil {
				return nil, fmt.Errorf("Failed to convert sampling_rate '%s' to float on line %d. Got Error: %w", values[valueIndex], parser.line, err)
			}
			newRead.SamplingRate = samplingRate
		case "len_raw_signal":
			lenRawSignal, err := strconv.ParseUint(values[valueIndex], 10, 64)
			if err != nil {
				return nil, fmt.Errorf("Failed to convert len_raw_signal '%s' to float on line %d. Got Error: %w", values[valueIndex], parser.line, err)
			}
			newRead.LenRawSignal = lenRawSignal
		case "raw_signal":
			var rawSignals []int16
			for rawSignalIndex, rawSignalString := range strings.Split(values[valueIndex], ",") {
				rawSignal, err := strconv.ParseInt(rawSignalString, 10, 16)
				if err != nil {
					return nil, fmt.Errorf("Failed to convert raw signal '%s' to int on line %d, signal index %d. Got Error: %w", rawSignalString, parser.line, rawSignalIndex, err)
				}
				rawSignals = append(rawSignals, int16(rawSignal))
			}
			newRead.RawSignal = rawSignals
		case "start_time":
			startTime, err := strconv.ParseUint(values[valueIndex], 10, 64)
			if err != nil {
				return nil, fmt.Errorf("Failed to convert start_time '%s' to uint on line %d. Got Error: %w", values[valueIndex], parser.line, err)
			}
			newRead.StartTime = startTime
		case "read_number":
			readNumber, err := strconv.ParseInt(values[valueIndex], 10, 32)
			if err != nil {
				return nil, fmt.Errorf("Failed to convert read_number '%s' to int on line %d. Got Error: %w", values[valueIndex], parser.line, err)
			}
			newRead.ReadNumber = int32(readNumber)
		case "start_mux":
			startMux, err := strconv.ParseUint(values[valueIndex], 10, 8)
			if err != nil {
				return nil, fmt.Errorf("Failed to convert start_mux '%s' to uint on line %d. Got Error: %w", values[valueIndex], parser.line, err)
			}
			newRead.StartMux = uint8(startMux)
		case "median_before":
			medianBefore, err := strconv.ParseFloat(values[valueIndex], 64)
			if err != nil {
				return nil, fmt.Errorf("Failed to convert median_before '%s' to float on line %d. Got Error: %w", values[valueIndex], parser.line, err)
			}
			newRead.MedianBefore = medianBefore
		case "end_reason":
			endReasonIndex, err := strconv.ParseInt(values[valueIndex], 10, 64)
			if err != nil {
				return nil, fmt.Errorf("Failed to convert end_reason '%s' to int on line %d. Got Error: %w", values[valueIndex], parser.line, err)
			}
			if _, ok := parser.endReasonMap[int(endReasonIndex)]; !ok {
				return nil, fmt.Errorf("End reason out of range. Got '%d' on line %d. Cannot find valid enum reason", int(endReasonIndex), parser.line)
			}
			newRead.EndReason = parser.endReasonMap[int(endReasonIndex)]
		case "channel_number":
			// For whatever reason, this is a string.
			newRead.ChannelNumber = values[valueIndex]
		default:
			return nil, fmt.Errorf("Unknown field to parser '%s' found on line %d. Please report to github.com/TimothyStiles/poly", fieldValue, parser.line)
		}
	}
	return &newRead, nil
}

/******************************************************************************
March 26, 2023

Start of Write functions

Slow5 write takes in a header, a channel of reads, and an io.Writer output.
The intended use case of slow5 write is reading from a location, such as a
database, and then directly writing the output to somewhere (stdout, a file,
etc).

A channel is used here so that reading and writing of slow5 can be done
concurrently. In almost all cases, you do not want to have slow5 files entirely
in memory, because they're freakin' huge.

Cheers,

Keoni

******************************************************************************/

func (header *Header) WriteTo(w io.Writer) (int64, error) {
	headers := header.HeaderValues
	var written int64
	// First, write the slow5 version number
	slow5Version := headers[0].Slow5Version
	endReasonHeaderMap := headers[0].EndReasonHeaderMap
	newWriteN, err := fmt.Fprintf(w, "#slow5_version\t%s\n", slow5Version)
	written += int64(newWriteN)
	if err != nil {
		return written, err
	}
	// Then, write the number of read groups (ie, the number of headers)
	newWriteN, err = fmt.Fprintf(w, "#num_read_groups\t%d\n", len(headers))
	written += int64(newWriteN)
	if err != nil {
		return written, err
	}
	// Next, we need a map of what attribute values are available
	possibleAttributeKeys := make(map[string]bool)
	for _, header := range headers {
		for key := range header.Attributes {
			possibleAttributeKeys[key] = true
		}
	}
	// Now that we know what attribute values are possible, lets build a map
	// with those values with "." placeholders (as defined in slow5 spec)
	headerAttributes := make(map[string][]string)
	for key := range possibleAttributeKeys {
		newBlankAttributes := make([]string, len(headers))
		for blankAttributeIndex := range newBlankAttributes {
			newBlankAttributes[blankAttributeIndex] = "."
		}
		headerAttributes[key] = newBlankAttributes
	}
	// Build a list with all header values
	var headerAttributeStrings []string
	for headerIndex, header := range headers {
		for key, value := range header.Attributes {
			headerAttributes[key][headerIndex] = value
		}
	}
	for key, value := range headerAttributes {
		var attributeRow strings.Builder
		attributeRow.Write([]byte(key))
		for _, attributeValue := range value {
			attributeRow.Write([]byte("\t"))
			attributeRow.Write([]byte(attributeValue))
		}
		headerAttributeStrings = append(headerAttributeStrings, attributeRow.String())
	}
	// Sort the header strings
	sort.Strings(headerAttributeStrings)

	// Write the header attribute strings to the output
	for _, headerAttributeString := range headerAttributeStrings {
		newWriteN, err = fmt.Fprintf(w, "%s\n", headerAttributeString)
		written += int64(newWriteN)
		if err != nil {
			return written, err
		}
	}

	// Now we handle endReasons / failureReasons. In the slow5 spec, these are
	// enums depending on what is present in the FAST5 files. This means the
	// labels may not be consistent and there is no exhaustive list of enum
	// labels. So, we have to create this from the endReasonMap each time we
	// write the slow5 file (not a const!)
	// Invert the endReasonMap
	endReasonStringList := make([]string, len(endReasonHeaderMap))
	for endReasonString, endReasonIndex := range endReasonHeaderMap {
		endReasonStringList[endReasonIndex] = endReasonString
	}
	// Build endReasonString
	var endReasonStringBuilder strings.Builder
	for _, endReason := range endReasonStringList {
		fmt.Fprintf(&endReasonStringBuilder, "%s,", endReason)
	}
	endReasonString := endReasonStringBuilder.String()
	endReasonString = endReasonString[:len(endReasonString)-1] // Remove trailing comma

	// Write the read headers
	// These are according to the slow5 specifications
	newWriteN, err = fmt.Fprintf(w, "#char*	uint32_t	double	double	double	double	uint64_t	int16_t*	uint64_t	int32_t	uint8_t	double	enum{%s}	char*\n", endReasonString)
	written += int64(newWriteN)
	if err != nil {
		return written, err
	}
	newWriteN, err = fmt.Fprintln(w, "#read_id	read_group	digitisation	offset	range	sampling_rate	len_raw_signal	raw_signal	start_time	read_number	start_mux	median_before	end_reason	channel_number")
	written += int64(newWriteN)
	if err != nil {
		return written, err
	}
	return written, nil
}

func (read *Read) WriteTo(w io.Writer) (int64, error) {
	var err error
	var written int64
	// converts []int16 to string
	var rawSignalStringBuilder strings.Builder
	for signalIndex, signal := range read.RawSignal {
		_, err = fmt.Fprint(&rawSignalStringBuilder, signal)
		if err != nil {
			return written, err
		}
		if signalIndex != len(read.RawSignal)-1 { // Don't add a comma to last number
			_, err = fmt.Fprint(&rawSignalStringBuilder, ",")
			if err != nil {
				return written, err
			}
		}
	}
	// Look at above output.Write("#read_id ... for the values here.
	newWriteN, err := fmt.Fprintf(w, "%s\t%d\t%g\t%g\t%g\t%g\t%d\t%s\t%d\t%d\t%d\t%g\t%d\t%s\n", read.ReadID, read.ReadGroupID, read.Digitisation, read.Offset, read.Range, read.SamplingRate, read.LenRawSignal, rawSignalStringBuilder.String(), read.StartTime, read.ReadNumber, read.StartMux, read.MedianBefore, read.EndReasonMap[read.EndReason], read.ChannelNumber)
	written += int64(newWriteN)
	if err != nil {
		return written, err
	}
	return written, nil
}

/******************************************************************************
Aug 15, 2023

StreamVByte (svb) compression of raw signal strength is used to decrease the
overall size of records in blow5 files. In my tests using a SQLite database
that contained both slow5 and fastq files, switching from raw signal TEXT to
svb compressed BLOBs brought the total database size from 12GB to 7.1GB, a
~40% reduction in size.

This is the primary method that can be used to decrease the size of storing
slow5 records in alternative datastores to blow5 (SQL databases, etc).

svb is an integer compression algorithm, which are specialized in compressing
integer arrays, which perfectly fits raw signal strength data. When converting
from slow5 to blow5, files are additionally compressed with zstd or zlib on top
of raw signal compression with svb. Still, svb compression is where most the
data size saving comes from.

Stream VByte paper: https://doi.org/10.48550/arXiv.1709.08990

Keoni

******************************************************************************/

// SvbCompressRawSignal takes a read and converts its raw signal field to two
// arrays: a mask array and a data array. Both are needed for decompression.
func SvbCompressRawSignal(rawSignal []int16) (mask, data []byte) {
	rawSignalUint32 := make([]uint32, len(rawSignal))
	for idx := range rawSignal {
		rawSignalUint32[idx] = uint32(rawSignal[idx])
	}
	return svb.Uint32Encode(rawSignalUint32)
}

// SvbDecompressRawSignal decompresses raw signal back to a []int16. It
// requires not only the mask array and data array returned by
// SvbCompressRawSignal, but also the length of the raw signals.
func SvbDecompressRawSignal(lenRawSignal int, mask, data []byte) []int16 {
	rawSignalUint32 := make([]uint32, lenRawSignal)
	rawSignal := make([]int16, lenRawSignal)
	svb.Uint32Decode32(mask, data, rawSignalUint32)
	for idx := 0; idx < lenRawSignal; idx++ {
		rawSignal[idx] = int16(rawSignalUint32[idx])
	}
	return rawSignal
}
