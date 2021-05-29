package csv_helper

import (
	"embed"
	"io/fs"
	"path/filepath"
	"regexp"
)

var csvFilesInDirectory []string

// appends `dirName` to `csvFilesInDirectory` if `dir` is a file
func fileName(dirName string, dir fs.DirEntry, e error) error {
	if e != nil {
		return e
	}
	if !dir.IsDir() {
		csvFilesInDirectory = append(csvFilesInDirectory, dirName)
	}
	return nil
}

// CSVFilesFromEmbededFS returns all csv files in a directory from an embeded FS
func CSVFilesFromEmbededFS(embededFS embed.FS, directory string) []string {
	csvFilesInDirectory = make([]string, 0)
	fs.WalkDir(embededFS, directory, fileName)

	re := regexp.MustCompile(`.*\.csv`)
	return filter(csvFilesInDirectory, re.MatchString)
}

// CSVFiles returns all csv files in a directory (and doesn't include csv files
// with 'processed' in its name)
func CSVFiles(directory string) []string {
	csvFilesInDirectory = make([]string, 0)
	// Walks directory and adds name of files to global var csvFilesInDirectory
	filepath.WalkDir(directory, fileName)

	// Filter to remove non .csv files
	re := regexp.MustCompile(`.*\.csv`)
	csvFiles := filter(csvFilesInDirectory, re.MatchString)

	// remove all csv files with 'processed' in its name
	re = regexp.MustCompile(`.*processed.*`)
	return filter(csvFiles, func(file string) bool { return !re.MatchString(file) })
}

func filter(ss []string, filter func(string) bool) (ret []string) {
	for _, s := range ss {
		if filter(s) {
			ret = append(ret, s)
		}
	}
	return
}
