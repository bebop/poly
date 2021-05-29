package rbs_calculator

import (
	"embed"
	"encoding/csv"
	"io"
	"os"
	"strconv"

	"github.com/TimothyStiles/poly/rbs_calculator/csv_helper"
)

//go:embed lookup_data/*
var embededRBSLookupDataDirectory embed.FS
var rbsLookupDataDirectory = "lookup_data"

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
	csvFiles := csv_helper.CSVFilesFromEmbededFS(embededRBSLookupDataDirectory, rbsLookupDataDirectory)

	for _, csv := range csvFiles {
		err := parseCSVIntoLookupTable(csv, &lookupTable)
		if err != nil {
			return nil, err
		}
	}
	return lookupTable, nil
}
