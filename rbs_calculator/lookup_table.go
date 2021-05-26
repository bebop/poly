package rbs_calculator

import (
	"encoding/csv"
	"io"
	"io/fs"
	"os"
	"path/filepath"
	"regexp"
	"strconv"
)

var csv_files []string

func GetFileName(s string, d fs.DirEntry, e error) error {
	if e != nil {
		return e
	}
	if !d.IsDir() {
		csv_files = append(csv_files, s)
	}
	return nil
}

func PopulateCsvFiles() {
	// Walks DATA_DIR and adds name of files to global var csv_files
	filepath.WalkDir(DATA_DIR, GetFileName)

	// Filter to remove non .csv files
	re := regexp.MustCompile(`.*\.csv`)
	csv_files = Filter(csv_files, re.MatchString)
}

var DATA_DIR string = "./data"

func Filter(ss []string, filter func(string) bool) (ret []string) {
	for _, s := range ss {
		if filter(s) {
			ret = append(ret, s)
		}
	}
	return
}

func Map(list []string, f func(string) string) []string {
	result := make([]string, len(list))
	for i, item := range list {
		result[i] = f(item)
	}
	return result
}

func parseValues(file string, lookupTable *map[string]map[string]float64) error {
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

		// dG := &float64{}
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
	PopulateCsvFiles()
	for _, csv := range csv_files {
		err := parseValues(csv, &lookupTable)
		if err != nil {
			return nil, err
		}
	}
	return lookupTable, nil
}
