package dataset

import (
	"encoding/csv"
	"io"
	"log"
	"os"

	"github.com/TimothyStiles/poly"
	"github.com/TimothyStiles/poly/rbs_calculator/csv_helper"
)

/**********************
Helper functions to add the folded structure of a sequence to a dataset
***********************/

func AddFoldedStructure(dataset, datasetOutputFile string, mRNAColNumber int) error {
	f, err := os.Open(dataset)
	if err != nil {
		return err
	}
	defer f.Close()
	csvr := csv.NewReader(f)

	wFile, err := csv_helper.CreateFile(datasetOutputFile)
	if err != nil {
		return err
	}
	defer wFile.Close()
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

	header = append(header, "STRUCTURE")
	if err := csvWriter.Write(header); err != nil {
		log.Fatalln("error writing record to file", err)
	}

	for {
		csvRow, err := csvr.Read()
		if err != nil {
			if err == io.EOF {
				err = nil
			}
			return err
		}

		mRNA := csvRow[mRNAColNumber]
		structure, _ := poly.LinearFold(mRNA, poly.DefaultBeamSize)
		csvRow = append(csvRow, structure)

		if err := csvWriter.Write(csvRow); err != nil {
			log.Fatalln("error writing record to file", err)
		}
	}
}
