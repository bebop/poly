package dataset

import (
	"sync"

	"github.com/TimothyStiles/poly/rbs_calculator/csv_helper"
)

// func ExampleGenerateVariables_1014IC() {
// 	fmt.Println("done")
// 	generateVariables("./with_folded_structure/1014IC.csv", "./processed_dataset/1014IC_vars.csv", 5, 3, 4, 8, 6, DefaultProperies)
// Output: AB000100
// }

// run with `go test -timeout 0 -run ^ExampleAddFoldedStructure_1014IC`
func ExampleAddFoldedStructure_1014IC() {
	datasetCSV := "./1014IC.csv"
	// The (1-indexed) row at which to start parsing datasetCSV from. Min is 0
	// (which means parsing the mRNA sequences from line number 1). Min: 0 Max:
	// Number of mRNA sequences to parse in CSV (Example: 1015 for `1014IC.csv`)
	datasetCSVLineNumber := 1015
	datasetWithStructureCSV := "./with_folded_structure/1014IC.csv"
	sequenceColNum, fivePrimeUTRColNum, cdsColNum := 5, 3, 4
	AddFoldedStructure(datasetCSV, datasetWithStructureCSV, datasetCSVLineNumber, sequenceColNum, fivePrimeUTRColNum, cdsColNum)
	// Random output to cause test to run

	// Output: 1
	// 2
}

// run with `go test -timeout 0 -run ^ExampleComputeProperties_1014IC`
func ExampleComputeProperties_1014IC() {
	// datasetName := "1014IC_1"
	datasetName := "1014IC"
	datasetCSV := "./" + datasetName + ".csv"
	var wg sync.WaitGroup

	// run this function and it outputs mRNAs to a channel
	wg.Add(1)
	mrnaChannel := make(chan MRNA)
	sequenceColNum, fivePrimeUTRColNum, cdsColNum, proteinMeanColNum, proteinStdColNum, orgnaismColNum, proteinColNum := 5, 3, 4, 6, 7, 1, 2
	go createDataset(datasetCSV, sequenceColNum, fivePrimeUTRColNum, cdsColNum, proteinMeanColNum, proteinStdColNum, orgnaismColNum, proteinColNum, mrnaChannel, &wg)

	// take the mRNA channel and then output to a computed properties channel
	wg.Add(1)
	csvOutputChannel := make(chan []string)
	go computeProperties2(mrnaChannel, DefaultProperies, csvOutputChannel, &wg)

	// take the computed properties channel and print to CSV
	wg.Add(1)
	datasetOutputFile := "./with_variables/" + datasetName + ".csv"
	go writeToCSV(datasetOutputFile, csv_helper.CREATE, csvOutputChannel, &wg)

	wg.Wait()

	// datasetWithProperties := computeProperties(dataset, DefaultProperies)
	// writeDatasetWithPropertiesToCSV(datasetWithProperties, datasetOutputFile)

	// Output: AB000100
	// Hello
}

func ExampleABD() {
	_, err := csv_helper.OpenFile("./with_variables/1014IC.csv", csv_helper.CREATE)
	if err != nil {
		panic(err)
	}

	// Output: AB000100
	// Hello
}
