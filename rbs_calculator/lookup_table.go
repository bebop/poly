package rbs_calculator

import (
	"encoding/csv"
	"io"
	"log"
	"os"
	"regexp"
	"strconv"
)

var DATA_DIR string = "./data"
var lookup_table map[string](map[string]float64)

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

func getCsvsInDataDir() []string {
	// get file names from DATA_DIR
	file, err := os.Open(DATA_DIR)
	if err != nil {
		log.Fatalf("failed opening directory: %s", err)
	}
	defer file.Close()
	files, _ := file.Readdirnames(0)

	// Filter to remove non .csv files
	re := regexp.MustCompile(`.*\.csv`)
	csvs := Filter(files, re.MatchString)

	// Add relative path to filename
	csvs = Map(csvs, func(file string) string { return DATA_DIR + "/" + file })
	return csvs
}

func parseValues(file string) error {
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

		if lookup_table[row[0]] == nil {
			lookup_table[row[0]] = make(map[string]float64)
		}

		lookup_table[row[0]][row[1]] = dG
	}
}

func Initalize() (map[string]map[string]float64, error) {
	lookup_table = make(map[string]map[string]float64)
	csvs := getCsvsInDataDir()
	for _, csv := range csvs {
		err := parseValues(csv)
		if err != nil {
			return nil, err
		}
	}
	return lookup_table, nil
}
