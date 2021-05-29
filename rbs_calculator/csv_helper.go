package rbs_calculator

import (
	"io/fs"
	"path/filepath"
	"regexp"
)

var csvFilesInDirectory []string

func fileName(s string, d fs.DirEntry, e error) error {
	if e != nil {
		return e
	}
	if !d.IsDir() {
		csvFilesInDirectory = append(csvFilesInDirectory, s)
	}
	return nil
}

func csvFiles(directory string) []string {
	// Walks directory and adds name of files to global var csvFilesInDirectory
	filepath.WalkDir(directory, fileName)

	// Filter to remove non .csv files
	re := regexp.MustCompile(`.*\.csv`)
	return filter(csvFilesInDirectory, re.MatchString)
}

func filter(ss []string, filter func(string) bool) (ret []string) {
	for _, s := range ss {
		if filter(s) {
			ret = append(ret, s)
		}
	}
	return
}
