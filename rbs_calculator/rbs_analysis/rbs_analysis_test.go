package rbs_analysis

import (
	"fmt"
)

func ExampleProcessCSV() {
	// csvFiles := csv_helper.CSVFiles(rbsAnalysisDataDirectory)
	// for _, csv := range csvFiles {
	// err := processCSV("./data/uASPIre_RBS_10k_r1.csv", 0, "rbs_calc_translation_rate", defaultOutputCSVFile)
	err := processCSV("./data/RBS_Calculator_v2.1.csv", 16, "rbs_calc_translation_rate", defaultOutputCSVFile)
	// err := processCSV("./data/salis_lab.csv", 4, "rbs_calc_translation_rate", defaultOutputCSVFile)
	if err != nil {
		fmt.Println(err)
	}
	// }
	fmt.Println("hi")
	// Output:
	// 11
}
