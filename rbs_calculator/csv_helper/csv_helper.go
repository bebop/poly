package csv_helper

import (
	"embed"
	"encoding/csv"
	"fmt"
	"io"
	"io/fs"
	"log"
	"os"
	"path/filepath"
	"regexp"
	"strings"
	"time"
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
// with strings in `disclude` in its name)
func CSVFiles(directory string, disclude []string) []string {
	csvFilesInDirectory = make([]string, 0)
	// Walks directory and adds name of files to global var csvFilesInDirectory
	filepath.WalkDir(directory, fileName)

	// Filter to remove non .csv files
	re := regexp.MustCompile(`.*\.csv`)
	csvFiles := filter(csvFilesInDirectory, re.MatchString)

	// remove all csv files with strings in `disclude` in its name
	for name := range disclude {
		re = regexp.MustCompile(fmt.Sprintf(".*%v.*", name))
		csvFiles = filter(csvFiles, func(file string) bool { return !re.MatchString(file) })
	}
	return csvFiles
}

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
	fileName = fileNameWithUNIXTimestamp(fileName)
	directorySeperatedFileName[len(directorySeperatedFileName)-1] = directory
	directorySeperatedFileName = append(directorySeperatedFileName, fileName)
	return strings.Join(directorySeperatedFileName, "/")
}

// returns the file name concatenated with the current unix timestamp before
// the file extension
// Example:
// file: "something.csv"
// output: "something_`current unix timestamp`.csv"
func fileNameWithUNIXTimestamp(file string) string {
	colonSeperatedFileName := strings.Split(file, ".")
	fileName := colonSeperatedFileName[len(colonSeperatedFileName)-2]
	colonSeperatedFileName[len(colonSeperatedFileName)-2] = fileName + fmt.Sprintf("_%v", time.Now().Unix())
	return strings.Join(colonSeperatedFileName, ".")
}

// CreateFile creates the file `file` and all the directories sepcified in
// `file`'s path if they don't exist
func CreateFile(file string) (*os.File, error) {
	if err := os.MkdirAll(filepath.Dir(file), 0770); err != nil {
		return nil, err
	}
	return os.Create(file)
}

// keepColumns saves the csv file `inputFile` as `outputFile` with only the
// specified columns
func KeepColumns(inputFile, outputFile string, cols []int) error {
	f, err := os.Open(inputFile)
	if err != nil {
		return err
	}
	defer f.Close()
	csvr := csv.NewReader(f)

	wFile, err := CreateFile(outputFile)
	if err != nil {
		return err
	}
	defer wFile.Close()
	csvWriter := csv.NewWriter(wFile)
	defer csvWriter.Flush()

	for {
		row, err := csvr.Read()
		if err != nil {
			if err == io.EOF {
				err = nil
			}
			return err
		}

		sublist := make([]string, 0)
		for _, col := range cols {
			sublist = append(sublist, row[col])
		}
		// row = append(row, fmt.Sprint(bindingSite.MRNAStructure))
		if err := csvWriter.Write(sublist); err != nil {
			log.Fatalln("error writing record to file", err)
		}
	}
}
