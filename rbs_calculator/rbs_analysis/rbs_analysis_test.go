package rbs_analysis

import (
	"fmt"

	"github.com/TimothyStiles/poly/rbs_calculator/csv_helper"
)

func ExampleProcessCSV() {
	csvFiles := csv_helper.CSVFiles(rbsAnalysisDataDirectory)
	for _, csv := range csvFiles {
		err := processCSV(csv, 0, defaultTranslationInitiationRateHeader, defaultOutputCSVFile)
		if err != nil {
			fmt.Println(err)
		}
	}
	fmt.Println("hi")
	// Output:
	// 11
}
