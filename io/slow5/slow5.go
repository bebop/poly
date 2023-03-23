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
	"errors"
	"fmt"
	"io"
	"strconv"
	"strings"
)

/******************************************************************************
Oct 10, 2021

slow5 parser begins here. Specification below:
https://hasindu2008.github.io/slow5specs/slow5-v0.1.0.pdf

slow5 is able to combine multiple sequencing runs into a single file format,
but we read each sequencing run separately. Each sequencing run contains a
header with metadata and a list of reads. However, unlike many other file
formats, a slow5 file should almost never be read into a common struct, since
most runs are large, and will take a ton of memory. Instead, the default way
to parse slow5 files is to produce a list of headers and a channel of raw
reads. In order to connect the two (if needed), use ReadGroupId.

The channel of raw reads can take advantage of Go's concurrency to do
multi-core processing of the raw reads.

Nanopore changes the attributes found in the header quite often, so we store
most of these attributes in a map for future proofing. Even the binary file
format, blow5, does not have types for these attributes, and just stores them
as a long string.

Reads have 8 required columns, and a few auxillary. These are typed, since they
are what will probably be used in real software.

Cheers mate,

Keoni

******************************************************************************/

// Header contains metadata about the sequencing run in general.
type Header struct {
	ReadGroupId  uint32
	Slow5Version string
	Attributes   map[string]string
}

// Read contains metadata and raw signal strengths for a single nanopore read.
type Read struct {
	ReadId       string
	ReadGroupId  uint32
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
}

// ParseHeader parses the header of a slow5 file.
func ParseHeader(r io.Reader) ([]Header, error) {
	var headers []Header
	scanner := bufio.NewScanner(r)
	var slow5Version string
	var numReadGroups uint32
	var lineNum int
	for scanner.Scan() {
		line := scanner.Text()
		values := strings.Split(line, "\t")
		if len(values) < 2 {
			return []Header{}, fmt.Errorf("Got following line without tabs: %s", line)
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
					return []Header{}, err
				}
				numReadGroups = uint32(numReadGroupsUint)
				for id := uint32(0); id < numReadGroups; id++ {
					emptyMap := make(map[string]string)
					headers = append(headers, Header{Slow5Version: slow5Version, ReadGroupId: id, Attributes: emptyMap})
				}
			}
			continue
		}
		// Terminate if we hit the beginning of the raw read headers
		if values[0] == "#char*" || values[0] == "#read_id" {
			return headers, nil
		}
		// Check to make sure we have the right amount of information for the num_read_groups
		if len(values) != int(numReadGroups+1) {
			return []Header{}, fmt.Errorf("Improper amount of information for read groups. Needed %d, got %d, in line: %s", numReadGroups+1, len(values), line)
		}
		for id := 0; id < int(numReadGroups); id++ {
			headers[id].Attributes[values[0]] = values[id+1]
		}
		lineNum++
	}
	return headers, errors.New("ReadGroup header never terminated. Improper file?")
}

var knownEndReasons = map[string]bool{"unknown": true,
	"partial":                         true,
	"mux_change":                      true,
	"unblock_mux_change":              true,
	"data_service_unblock_mux_change": true,
	"signal_positive":                 true,
	"signal_negative":                 true,
}

// ParseReads parses all reads in a slow5 file.
func ParseReads(r io.Reader, reads chan<- Read, errorsChan chan<- error) {
	headerMap := make(map[int]string)
	endReasonMap := make(map[int]string)
	scanner := bufio.NewScanner(r)
	start := false
	lineNum := 0
	for scanner.Scan() {
		line := scanner.Text()
		values := strings.Split(line, "\t")
		if !start {
			// Get endReasonEnums. This is simply a string between enum{} that is used for the reasons that a read could have ended.
			if values[0] == "#char*" {
				for _, typeInfo := range values {
					if strings.Contains(typeInfo, "enum") {
						endReasonEnumsMinusPrefix := strings.TrimPrefix(typeInfo, "enum{")
						endReasonEnumsMinusSuffix := strings.TrimSuffix(endReasonEnumsMinusPrefix, "}")
						endReasons := strings.Split(endReasonEnumsMinusSuffix, ",")

						for endReasonIndex, endReason := range endReasons {
							if _, ok := knownEndReasons[endReason]; !ok {
								errorsChan <- fmt.Errorf("Unknown end reason '%s' found in end_reason enum. Please report.", endReason)
							}
							endReasonMap[endReasonIndex] = endReason
						}
					}
				}
			}

			// Get the read headers and their identifiers. Though the primary read headers are in a defined order, the auxiliary headers are not.
			if values[0] == "#read_id" {
				headerMap[0] = "read_id"
				for headerNum := 1; headerNum < len(values); headerNum++ {
					headerMap[headerNum] = values[headerNum]
				}
				start = true
			}
			continue
		}

		// Once we have the read headers, start to parse the actual reads
		var newRead Read
		for valueIndex := 0; valueIndex < len(values); valueIndex++ {
			fieldValue := headerMap[valueIndex]
			if values[valueIndex] == "." {
				continue
			}
			switch fieldValue {
			case "read_id":
				newRead.ReadId = values[valueIndex]
			case "read_group":
				readGroupId, err := strconv.ParseUint(values[valueIndex], 10, 32)
				if err != nil {
					errorsChan <- fmt.Errorf("Failed convert read_group '%s' to uint on line %d. Got error: %s", values[valueIndex], lineNum, err)
				}
				newRead.ReadGroupId = uint32(readGroupId)
			case "digitisation":
				digitisation, err := strconv.ParseFloat(values[valueIndex], 64)
				if err != nil {
					errorsChan <- fmt.Errorf("Failed to convert digitisation '%s' to float on line %d. Got error: %s", values[valueIndex], lineNum, err)
				}
				newRead.Digitisation = digitisation
			case "offset":
				offset, err := strconv.ParseFloat(values[valueIndex], 64)
				if err != nil {
					errorsChan <- fmt.Errorf("Failed to convert offset '%s' to float on line %d. Got error: %s", values[valueIndex], lineNum, err)
				}
				newRead.Offset = offset
			case "range":
				nanoporeRange, err := strconv.ParseFloat(values[valueIndex], 64)
				if err != nil {
					errorsChan <- fmt.Errorf("Failed to convert range '%s' to float on line %d. Got error: %s", values[valueIndex], lineNum, err)
				}
				newRead.Range = nanoporeRange
			case "sampling_rate":
				samplingRate, err := strconv.ParseFloat(values[valueIndex], 64)
				if err != nil {
					errorsChan <- fmt.Errorf("Failed to convert sampling_rate '%s' to float on line %d. Got error: %s", values[valueIndex], lineNum, err)
				}
				newRead.SamplingRate = samplingRate
			case "len_raw_signal":
				lenRawSignal, err := strconv.ParseUint(values[valueIndex], 10, 64)
				if err != nil {
					errorsChan <- fmt.Errorf("Failed to convert len_raw_signal '%s' to float on line %d. Got error: %s", values[valueIndex], lineNum, err)
				}
				newRead.LenRawSignal = lenRawSignal
			case "raw_signal":
				var rawSignals []int16
				for rawSignalIndex, rawSignalString := range strings.Split(values[valueIndex], ",") {
					rawSignal, err := strconv.ParseInt(rawSignalString, 10, 16)
					if err != nil {
						errorsChan <- fmt.Errorf("Failed to convert raw signal '%s' to int on line %d, signal index %d. Got error: %s", rawSignalString, lineNum, rawSignalIndex, err)
					}
					rawSignals = append(rawSignals, int16(rawSignal))
				}
				newRead.RawSignal = rawSignals
			case "start_time":
				startTime, err := strconv.ParseUint(values[valueIndex], 10, 64)
				if err != nil {
					errorsChan <- fmt.Errorf("Failed to convert start_time '%s' to uint on line %d. Got error: %s", values[valueIndex], lineNum, err)
				}
				newRead.StartTime = startTime
			case "read_number":
				readNumber, err := strconv.ParseInt(values[valueIndex], 10, 32)
				if err != nil {
					errorsChan <- fmt.Errorf("Failed to convert read_number '%s' to int on line %d. Got error: %s", values[valueIndex], lineNum, err)
				}
				newRead.ReadNumber = int32(readNumber)
			case "start_mux":
				startMux, err := strconv.ParseUint(values[valueIndex], 10, 8)
				if err != nil {
					errorsChan <- fmt.Errorf("Failed to convert start_mux '%s' to uint on line %d. Got error: %s", values[valueIndex], lineNum, err)
				}
				newRead.StartMux = uint8(startMux)
			case "median_before":
				medianBefore, err := strconv.ParseFloat(values[valueIndex], 64)
				if err != nil {
					errorsChan <- fmt.Errorf("Failed to convert median_before '%s' to float on line %d. Got error: %s", values[valueIndex], lineNum, err)
				}
				newRead.MedianBefore = medianBefore
			case "end_reason":
				endReasonIndex, err := strconv.ParseInt(values[valueIndex], 10, 64)
				if err != nil {
					errorsChan <- fmt.Errorf("Failed to convert end_reason '%s' to int on line %d. Got error: %s", values[valueIndex], lineNum, err)
				}
				if _, ok := endReasonMap[int(endReasonIndex)]; !ok {
					errorsChan <- fmt.Errorf("End reason out of range. Got '%d' on line %d. Cannot find valid enum reason", int(endReasonIndex), lineNum)
				}
				newRead.EndReason = endReasonMap[int(endReasonIndex)]
			case "channel_number":
				// For whatever reason, this is a string.
				newRead.ChannelNumber = values[valueIndex]
			default:
				errorsChan <- fmt.Errorf("Unknown field to parser '%s' found on line %d. Please report to github.com/TimothyStiles/poly", fieldValue, lineNum)
			}
		}
		reads <- newRead

		lineNum++
	}
	close(reads)
	close(errorsChan)
}
