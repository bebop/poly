package rbs_calculator

import (
	"encoding/csv"
	"io"
	"os"
	"strconv"
)

var rbsDataDirectory string = "./data"

func parseCSVIntoLookupTable(file string, lookupTable *map[string]map[string]float64) error {
	f, err := os.Open(file)
	if err != nil {
		return err
	}
	defer f.Close()

	csvr := csv.NewReader(f)

	// consume header line
	csvr.Read()

	for {
		row, err := csvr.Read()
		if err != nil {
			if err == io.EOF {
				err = nil
			}
			return err
		}

		var dG float64
		if dG, err = strconv.ParseFloat(row[2], 64); err != nil {
			return err
		}

		if (*lookupTable)[row[0]] == nil {
			(*lookupTable)[row[0]] = make(map[string]float64)
		}

		(*lookupTable)[row[0]][row[1]] = dG
	}
}

func LookupTable() (map[string]map[string]float64, error) {
	lookupTable := make(map[string]map[string]float64)
	csvFiles := csvFiles(rbsDataDirectory)
	for _, csv := range csvFiles {
		err := parseCSVIntoLookupTable(csv, &lookupTable)
		if err != nil {
			return nil, err
		}
	}
	return lookupTable, nil
}
