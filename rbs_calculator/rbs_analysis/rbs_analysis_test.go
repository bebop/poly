package rbs_analysis

import (
	"fmt"

	"github.com/TimothyStiles/poly/rbs_calculator"
)

func ExampleProcessCSV1() {
	dG_rRNA_mRNA, err := rbs_calculator.LookupTable()
	if err != nil {
		fmt.Printf("Failed to initialize lookup table: %s", err)
		return
	}
	fmt.Println(dG_rRNA_mRNA["ACCTCCTTA"]["CAAGGAGGGTG"])
	// Output: -11.9
}

func ExampleProcessCSV() {
	// csvFiles := csv_helper.CSVFiles(rbsAnalysisDataDirectory)
	// for _, csv := range csvFiles {
	err := processCSV("./data/salis_lab.csv", 6, "rbs_calc_translation_rate", defaultOutputCSVFile)
	if err != nil {
		fmt.Println(err)
	}
	// }
	fmt.Println("hi")
	// Output:
	// 11
}
