package rbs_analysis

import (
	"encoding/csv"
	"fmt"
	"io"
	"log"
	"os"
	"path/filepath"
	"strings"
	"time"

	"github.com/TimothyStiles/poly/rbs_calculator"
)

// A file that contains helper functions needed to compare poly's RBS calculator
// with other implementations of RBS calculators and with data from experiments

var defaultTranslationInitiationRateHeader = "rate"
var defaultOutputCSVFile = "DEFAULT"
var rbsAnalysisDataDirectory = "./data"

// processCSV reads a CSV (skipping the first header line), and passes mRNA sequence
// at column number 0-indexed `mRNAColNumber` to `RibosomeBindingSite()` with rRNA
// of `EColiRNA`. The output is then written to the last available column in the
// CSV under the header `translationInitiationRateHeader` and saved in `csvOutputFile`.
func processCSV(file string, mRNAColNumber int, translationInitiationRateHeader string, csvOutputFile string) error {
	if csvOutputFile == defaultOutputCSVFile {
		csvOutputFile = csvOutputFileName(file)
	}

	f, err := os.Open(file)
	if err != nil {
		return err
	}
	defer f.Close()
	csvr := csv.NewReader(f)

	wFile, err := createFile(csvOutputFile)
	if err != nil {
		return err
	}
	defer f.Close()
	csvWriter := csv.NewWriter(wFile)
	defer csvWriter.Flush()

	// consume header line
	header, err := csvr.Read()
	if err != nil {
		if err == io.EOF {
			err = nil
		}
		return err
	}
	header = append(header, translationInitiationRateHeader)
	header = append(header, "dG_total")
	header = append(header, "dG_mRNA")
	header = append(header, "dG_mRNA_rRNA")
	header = append(header, "bindingSiteSequence")
	header = append(header, "foundBindingSite")
	// header = append(header, "mRNAStructure")
	if err := csvWriter.Write(header); err != nil {
		log.Fatalln("error writing record to file", err)
	}

	for {
		row, err := csvr.Read()
		if err != nil {
			if err == io.EOF {
				err = nil
			}
			return err
		}

		// parse values in each csv line
		rRNA := rbs_calculator.EColiRNA
		mRNA := row[mRNAColNumber]
		bindingSite, err := rbs_calculator.RibosomeBindingSite(rRNA, mRNA)
		if err != nil {
			return err
		}
		row = append(row, fmt.Sprint(bindingSite.TranslationInitiationRate))
		row = append(row, fmt.Sprint(bindingSite.MinimumFreeEnergy))
		row = append(row, fmt.Sprint(bindingSite.MRNAFreeEnergy))
		row = append(row, fmt.Sprint(bindingSite.BindingSiteFreeEnergy))

		bindingSiteSequence := bindingSite.MRNA[bindingSite.FivePrimeIdx:bindingSite.ThreePrimeIdx]
		bindingSiteSequence = fmt.Sprintf("%v (%v,%v)", bindingSiteSequence, bindingSite.FivePrimeIdx, bindingSite.ThreePrimeIdx)
		row = append(row, fmt.Sprint(bindingSiteSequence))

		row = append(row, fmt.Sprint(bindingSite.FoundBindingSite))
		// row = append(row, fmt.Sprint(bindingSite.MRNAStructure))
		if err := csvWriter.Write(row); err != nil {
			log.Fatalln("error writing record to file", err)
		}
	}
}

// returns the file name concatenated with the current unix timestamp before
// the file extension
func processedFileName(file string) string {
	colonSeperatedFileName := strings.Split(file, ".")
	fileName := colonSeperatedFileName[len(colonSeperatedFileName)-2]
	colonSeperatedFileName[len(colonSeperatedFileName)-2] = fileName + fmt.Sprintf("_%v", time.Now().Unix())
	return strings.Join(colonSeperatedFileName, ".")
}

func csvOutputFileName(inputFile string) string {
	directorySeperatedFileName := strings.Split(inputFile, "/")
	fileName := directorySeperatedFileName[len(directorySeperatedFileName)-1]
	fileName = processedFileName(fileName)
	directorySeperatedFileName[len(directorySeperatedFileName)-1] = "processed"
	directorySeperatedFileName = append(directorySeperatedFileName, fileName)
	return strings.Join(directorySeperatedFileName, "/")
}

func createFile(p string) (*os.File, error) {
	if err := os.MkdirAll(filepath.Dir(p), 0770); err != nil {
		return nil, err
	}
	return os.Create(p)
}
