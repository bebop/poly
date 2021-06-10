package dataset

import "sync"

// func ExampleGenerateVariables_1014IC() {
// 	fmt.Println("done")
// 	generateVariables("./processed_dataset/1014IC.csv", "./processed_dataset/1014IC_vars.csv", 5, 3, 4, 8, 6, DefaultProperies)
// Output: AB000100
// }

func ExampleAddFoldedStructure_1014IC() {
	datasetCSV := "./1014IC.csv"
	// Column number of mRNA sequence in datasetCSV
	sequenceColNum := 5
	datasetWithStructureCSV := "./processed_dataset/1014IC.csv"
	AddFoldedStructure(datasetCSV, datasetWithStructureCSV, sequenceColNum)
}

func ExampleComputeProperties_1014IC() {
	datasetWithStructureCSV := "./processed_dataset/1014IC.csv"
	sequenceColNum, fivePrimeUTRColNum, cdsColNum, structureColNum, outputColNum := 5, 3, 4, 8, 6

	var wg sync.WaitGroup

	// run this function and it outputs mRNAs to a channel
	wg.Add(1)
	mrnaChannel := make(chan MRNA)
	go func() {
		createDataset(datasetWithStructureCSV, sequenceColNum, fivePrimeUTRColNum, cdsColNum, structureColNum, outputColNum, mrnaChannel)
		wg.Done()
	}()

	// take the mRNA channel and then output to a computed properties channel
	wg.Add(1)
	csvOutputChannel := make(chan []string)
	go func() {
		computeProperties2(mrnaChannel, DefaultProperies, csvOutputChannel)
		wg.Done()
	}()

	// take the computed properties channel and print to CSV
	wg.Add(1)
	datasetOutputFile := "./processed_dataset/1014IC_vars.csv"
	go func() {
		writeToCSV(datasetOutputFile, csvOutputChannel)
		wg.Done()
	}()

	wg.Wait()

	// datasetWithProperties := computeProperties(dataset, DefaultProperies)
	// writeDatasetWithPropertiesToCSV(datasetWithProperties, datasetOutputFile)

	// Output: AB000100
	// Hello
}
