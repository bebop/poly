package csv_helper

import (
	"embed"
	"encoding/csv"
	"fmt"
	"io"
	"io/fs"
	"os"
	"path/filepath"
	"regexp"
	"strings"
	"sync"
	"time"
)

/*****************************************************************************

This file contains helper functions to work with CSV files.
The functions defined here are used by the `model`, and
`shine_dalgarno_binding_site` subpackages.

*****************************************************************************/

// CSVFilesFromEmbeddedFS returns all csv files in a directory from an embedded FS
func CSVFilesFromEmbeddedFS(embeddedFS embed.FS, directory string) (csvFilesInDirectory []string) {
	filesInDirectory = make([]string, 0)
	fs.WalkDir(embeddedFS, directory, fileName)

	csvFileRegex := regexp.MustCompile(`.*\.csv`)
	csvFilesInDirectory = filter(filesInDirectory, csvFileRegex.MatchString)
	return
}

// appends `dirName` to global variable `filesInDirectory` if `dir` is a file
func fileName(dirName string, dir fs.DirEntry, e error) error {
	if e != nil {
		return e
	}

	if !dir.IsDir() {
		filesInDirectory = append(filesInDirectory, dirName)
	}
	return nil
}

// When walking a directory with `fs.WalkDir` (in the func
// `CSVFilesFromEmbeddedFS`), we use the func `fileName` as the function that
// gets called for every directory. Unfortunately, the walking func cannot
// return anything so we have to use this global un-exported variable to
// keep track of the files in the directory
var filesInDirectory []string

// returns all strings that return true for the `filter` func
func filter(ss []string, filter func(string) bool) (ret []string) {
	for _, s := range ss {
		if filter(s) {
			ret = append(ret, s)
		}
	}
	return
}

// CSVOutputFileName returns name of csv file suffixed with the current unix
// timestamp and with the directory `directory` prefixed to its path
// Example:
// inputFile: "./data/something.csv"
// output: "./data/`directory`/something_`current unix timestamp`.csv"
func CSVOutputFileName(directory, inputFile string) string {
	directorySeperatedFileName := strings.Split(inputFile, "/")
	fileName := directorySeperatedFileName[len(directorySeperatedFileName)-1]
	fileName = FileNameWithUNIXTimestamp(fileName)
	directorySeperatedFileName[len(directorySeperatedFileName)-1] = directory
	directorySeperatedFileName = append(directorySeperatedFileName, fileName)
	return strings.Join(directorySeperatedFileName, "/")
}

// returns the file name concatenated with the current unix timestamp before
// the file extension
// Example:
// file: "something.csv"
// output: "something_`current unix timestamp`.csv"
func FileNameWithUNIXTimestamp(file string) string {
	colonSeperatedFileName := strings.Split(file, ".")
	fileName := colonSeperatedFileName[len(colonSeperatedFileName)-2]
	colonSeperatedFileName[len(colonSeperatedFileName)-2] = fileName + fmt.Sprintf("_%v", time.Now().Unix())
	return strings.Join(colonSeperatedFileName, ".")
}

// Re-export the package `os`'s constants so that packages that depend on
// csv_helper's `OpenFile` func don't need to depend on package `os` as well
const (
	APPEND int = os.O_APPEND
	CREATE int = os.O_CREATE
)

// OpenFile creates the file `file` and all the directories specified in
// `file`'s path if they don't exist.
// The available options for `osFlag` are `APPEND` (which will open the
// specified file without overwriting it) and `CREATE` (which will overwrite
// the specified file if it exists)
func OpenFile(file string, osFlag int) (*os.File, error) {
	if err := os.MkdirAll(filepath.Dir(file), 0770); err != nil {
		return nil, err
	}

	flags := os.O_RDWR | osFlag
	return os.OpenFile(file, flags, 0770)
}

// WriteToCSV writes to the specified output file with the contents from
// the channel `csvOutputChannel`. Writing is stopped when `csvOutputChannel`
// is closed. Please read the doc of `OpenFile` to understand the available
// options for the `osFlags` argument.
func WriteToCSV(outputFile string, osFlags int, csvOutputChannel chan []string, wg *sync.WaitGroup) {
	wFile, err := OpenFile(outputFile, osFlags)
	if err != nil {
		panic(err)
	}
	defer wFile.Close()
	csvWriter := csv.NewWriter(wFile)
	defer csvWriter.Flush()
	defer wg.Done()

	for {
		select {
		case csvLine, ok := <-csvOutputChannel:
			if !ok {
				return
			}

			if err := csvWriter.Write(csvLine); err != nil {
				panic(err)
			}
		}
	}
}

// ReadCSV reads a csv file into the channel `csvRowsChannel`.
func ReadCSV(fs embed.FS, csvFile string, csvRowsChannel chan []string,
	wg *sync.WaitGroup, nbSkipHeader int) {
	defer close(csvRowsChannel)
	defer wg.Done()

	f, err := fs.Open(csvFile)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	csvReader := csv.NewReader(f)

	for i := 0; i < nbSkipHeader; i++ {
		csvReader.Read()
	}

	csvRow, err := csvReader.Read()
	for err == nil {
		csvRowsChannel <- csvRow
		csvRow, err = csvReader.Read()
	}

	if err != io.EOF {
		panic(err)
	}
}

// ReadHeader returns the specified number of rows from the start of a csv file
func ReadHeader(fs embed.FS, csvFile string, nbHeaderRows int) (header [][]string) {
	f, err := fs.Open(csvFile)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	csvReader := csv.NewReader(f)

	for i := 0; i < nbHeaderRows; i++ {
		csvRow, err := csvReader.Read()
		if err != nil {
			if err != io.EOF {
				panic(err)
			} else {
				panic("reached EOF before parsing all header rows")
			}
		}
		header = append(header, csvRow)
	}
	return
}
