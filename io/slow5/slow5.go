package slow5

import (
	"bufio"
	"errors"
	"fmt"
	"io"
	"strconv"
	"strings"
	"time"
)

type ReadGroup struct {
	ReadGroupId uint32

	Slow5Version               string
	AsicId                     int `json:"asic_id"`
	asic_id_eeprom             int
	asic_temp                  float64
	asic_version               string
	auto_update                bool
	auto_update_source         string
	barcoding_enabled          bool
	bream_is_standard          bool
	configuration_version      string
	device_id                  string
	device_type                string
	distribution_status        string
	distribution_version       string
	exp_script_name            string
	exp_script_purpose         string
	exp_start_time             time.Time
	experiment_duration_set    int
	experiment_type            string
	file_type                  string
	file_version               string
	flongle_adapter_id         string
	flow_cell_id               string
	flow_cell_product_code     string
	guppy_version              string
	heatsink_temp              float64
	host_product_code          string
	host_product_serial_number string
	hostname                   string
	installation_type          string
	local_basecalling          bool
	local_firmware_file        int // no information on what this exactly is
	operating_system           string
	Package                    string
	package_version            string
	pore_type                  string
	protocol_group_id          string
	protocol_run_id            string
	protocol_start_time        time.Time
	protocols_version          string
	run_id                     string
	sample_frequency           int
	sample_id                  string
	sequencing_kit             string
	usb_config                 string
	version                    string
}

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

func ParseReadGroups(r io.Reader) ([]ReadGroup, error) {
	var readGroups []ReadGroup
	scanner := bufio.NewScanner(r)
	var slow5Version string
	var numReadGroups uint32
	var lineNum int
	for scanner.Scan() {
		line := scanner.Text()
		values := strings.Split(line, "\t")
		if len(values) < 2 {
			return []ReadGroup{}, fmt.Errorf("Got following line without tabs: %s", line)
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
					return []ReadGroup{}, err
				}
				numReadGroups = uint32(numReadGroupsUint)
				for id := uint32(0); id < numReadGroups; id++ {
					readGroups = append(readGroups, ReadGroup{Slow5Version: slow5Version, ReadGroupId: id})
				}
			}
		} else {
			// Terminate if we hit the beginning of the raw read headers
			if values[0] == "#char*" || values[0] == "#read_id" {
				return readGroups, nil
			}
			// Check to make sure we have the right amount of information for the num_read_groups
			if len(values) != int(numReadGroups+1) {
				return []ReadGroup{}, fmt.Errorf("Improper amount of information for read groups. Needed %d, got %d, in line: %s", numReadGroups+1, len(values), line)
			}
			switch values[0] {
			case "@asic_id":
				for id := 0; id < int(numReadGroups); id++ {
					asic_id, err := strconv.ParseInt(values[id+1], 10, 64)
					if err != nil {
						return []ReadGroup{}, fmt.Errorf("Failed to convert asic_id '%s' to int on line %d. Got error: %s", values[id+1], lineNum, err)
					}
					readGroups[id].AsicId = int(asic_id)
				}
			}
		}
		lineNum++
	}
	return readGroups, errors.New("ReadGroup header never terminated. Improper file?")
}

var knownEndReasons = map[string]bool{"unknown": true,
	"partial":                         true,
	"mux_change":                      true,
	"unblock_mux_change":              true,
	"data_service_unblock_mux_change": true,
	"signal_positive":                 true,
	"signal_negative":                 true,
}

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
				continue
			}
		}

		// Once we have the read headers, start to parse the actual reads
		if start {
			var newRead Read
			for valueIndex := 0; valueIndex < len(values); valueIndex++ {
				fieldValue := headerMap[valueIndex]
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
					newRead.EndReason = endReasonMap[int(endReasonIndex)]
				case "channel_number":
					// For whatever reason, this is a string.
					newRead.ChannelNumber = values[valueIndex]
				default:
					errorsChan <- fmt.Errorf("Unknown field to parser '%s' found on line %d. Please report to github.com/TimothyStiles/poly", fieldValue, lineNum)
				}
			}
			reads <- newRead
		}
		lineNum++
	}
	close(reads)
	close(errorsChan)
}
