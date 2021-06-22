package dataset

import (
	"encoding/csv"
	"fmt"
	"io"
	"os"
	"reflect"
	"runtime"
	"strconv"
	"strings"
	"sync"

	"github.com/TimothyStiles/poly"
	"github.com/TimothyStiles/poly/rbs_calculator/csv_helper"
)

/**********************
Step 1: Add folded structrue to dataset

Helper functions to add the folded structure of a sequence to a dataset
***********************/

func AddFoldedStructure(dataset, datasetOutputFile string, csvStartRowIdx, sequenceColNum, fivePrimeUTRColNum, cdsColNum int) {
	if csvStartRowIdx < 1 {
		panic("csvStartRowIdx should be >= 1")
	}

	var wg sync.WaitGroup

	wg.Add(1)
	csvOutputChannel := make(chan []string)
	if csvStartRowIdx == 1 {
		go writeToCSV(datasetOutputFile, csv_helper.CREATE, csvOutputChannel, &wg)
	} else {
		go writeToCSV(datasetOutputFile, csv_helper.APPEND, csvOutputChannel, &wg)
	}

	defer close(csvOutputChannel)
	defer wg.Wait()

	f, err := os.Open(dataset)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	csvr := csv.NewReader(f)

	// read header line only if starting at first row
	if csvStartRowIdx == 1 {

		// read header line
		header, err := csvr.Read()
		if err != nil {
			if err == io.EOF {
				err = nil
			}
			panic(err)
		}

		header = append(header, "STRUCTURE")
		// header := getMRNAFields()
		csvOutputChannel <- header
	}

	for i := 1; i < csvStartRowIdx; i++ {
		csvr.Read()
	}

	csvRow, err := csvr.Read()
	for err == nil {
		for _, rnaColNum := range []int{sequenceColNum, fivePrimeUTRColNum, cdsColNum} {
			csvRow[rnaColNum] = cleanRNA(csvRow[rnaColNum])
		}

		structure, _ := poly.LinearFold(csvRow[sequenceColNum])
		panic(structure)

		csvRow = append(csvRow, structure)

		csvOutputChannel <- csvRow

		csvRow, err = csvr.Read()
	}

	if err != io.EOF {
		panic(err)
	}
}

/**********************
Step 2: Add variables of model to visualize and explore

Functions to add variables to visualize and explore to the dataset
***********************/

// MRNA contains all the 'ground-truth' information of a mRNA sequence
type MRNA struct {
	Sequence string
	// Structure             string
	FivePrimeUTR          string
	ProteinCodingSequence string
	ProteinMean           float64
	ProteinStd            float64
	Organism              string
	Protein               string
	properties            map[string]interface{}
	// energyContributions   []mfe.EnergyContribution
}

// Dataset is a list of 'ground-truth' mRNA sequences
type Dataset = []MRNA

// cleanRNA converts DNA to RNA and makes sequence upper case
func cleanRNA(rna string) string {
	rna = strings.ToUpper(rna)
	rna = strings.ReplaceAll(rna, "T", "U")
	return rna
}

// createDataset parses a csv to a Dataset struct
func createDataset(csvFile string, sequenceColNum, fivePrimeUTRColNum, cdsColNum, proteinMeanColNum, proteinStdColNum, orgnaismColNum, proteinColNum int, mrnaChannel chan MRNA, wg *sync.WaitGroup) {
	defer close(mrnaChannel)
	defer wg.Done()

	f, err := os.Open(csvFile)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	csvr := csv.NewReader(f)

	// consume header line
	csvr.Read()

	// var mrnas []MRNA

	csvRow, err := csvr.Read()
	for err == nil {
		// fmt.Println(i)
		// log.Println(i)
		var proteinMean, proteinStd float64 = 0.0, 0.0
		if proteinMean, err = strconv.ParseFloat(csvRow[proteinMeanColNum], 64); err != nil {
			panic(err)
		}
		if proteinStd, err = strconv.ParseFloat(csvRow[proteinStdColNum], 64); err != nil {
			if csvRow[proteinStdColNum] != "" {
				// only panic if there is value in row, else keep value as 0
				panic(err)
			}
		}
		// sequence, structure := cleanRNA(csvRow[sequenceColNum]), csvRow[structureColNum]
		// _, energyContributions, err := mfe.MinimumFreeEnergy(sequence, structure, mfe.DefaultTemperature)
		// if err != nil {
		// 	panic(err)
		// }
		mRNA := MRNA{
			Sequence: cleanRNA(csvRow[sequenceColNum]),
			// Structure:             csvRow[structureColNum],
			ProteinMean:           proteinMean,
			FivePrimeUTR:          cleanRNA(csvRow[fivePrimeUTRColNum]),
			ProteinCodingSequence: cleanRNA(csvRow[cdsColNum]),
			ProteinStd:            proteinStd,
			Organism:              csvRow[orgnaismColNum],
			Protein:               csvRow[proteinColNum],
			properties:            make(map[string]interface{}),
			// energyContributions:   energyContributions,
		}
		// mrnas = append(mrnas, mRNA)
		mrnaChannel <- mRNA
		csvRow, err = csvr.Read()
	}

	if err != io.EOF {
		panic(err)
	}
	// return
}

// createDataset parses a csv to a Dataset struct
func streamDataset(csvFile string, sequenceColNum, fivePrimeUTRColNum, cdsColNum, outputColNum, batchSize, batchNum int) (Dataset, bool) {
	f, err := os.Open(csvFile)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	csvr := csv.NewReader(f)

	// consume header line only for the first batch
	if batchNum == 0 {
		csvr.Read()
	}

	var mrnas []MRNA

	csvRow, err := csvr.Read()
	for i := batchSize; err == nil && i < batchSize; i++ {
		// fmt.Println(i)
		// log.Println(i)
		var output float64
		if output, err = strconv.ParseFloat(csvRow[outputColNum], 64); err != nil {
			panic(err)
		}
		// sequence, structure := cleanRNA(csvRow[sequenceColNum]), csvRow[structureColNum]
		// _, energyContributions, err := mfe.MinimumFreeEnergy(sequence, structure, mfe.DefaultTemperature)
		if err != nil {
			panic(err)
		}
		mRNA := MRNA{
			Sequence: cleanRNA(csvRow[sequenceColNum]),
			// Structure:             csvRow[structureColNum],
			ProteinMean:           output,
			FivePrimeUTR:          cleanRNA(csvRow[fivePrimeUTRColNum]),
			ProteinCodingSequence: cleanRNA(csvRow[cdsColNum]),
			// energyContributions:   energyContributions,
		}
		mrnas = append(mrnas, mRNA)
		csvRow, err = csvr.Read()
	}

	if err == io.EOF {
		return mrnas, true
	}
	panic(err)
}

// MRNAWithProperties contains a 'ground-truth' MRNA with its computed property
// values
type MRNAWithProperties struct {
	mrna MRNA
	// list all properties that need to be analyzed here
	// propertyNames  []string
	propertyValues []string
}

// DatasetWithProperties is a Dataset with the computed property values
type DatasetWithProperties struct {
	mrnaWithProperties []MRNAWithProperties
	propertyNames      []string
}

// func computeProperties(dataset Dataset, properties [](func(MRNA) string)) DatasetWithProperties {
// 	var mrnasWithProperties []MRNAWithProperties

// 	for _, mrna := range dataset {
// 		var propertyValues []string
// 		for _, property := range properties {
// 			propertyValues = append(propertyValues, property(mrna))
// 		}

// 		mrnaWithProperties := MRNAWithProperties{
// 			mrna:           mrna,
// 			propertyValues: propertyValues,
// 		}
// 		mrnasWithProperties = append(mrnasWithProperties, mrnaWithProperties)
// 	}

// 	return DatasetWithProperties{
// 		mrnaWithProperties: mrnasWithProperties,
// 		propertyNames:      GetFunctionNames(properties),
// 	}
// }

func computeProperties2(mrnaChannel chan MRNA, properties []Property, csvOutputChannel chan []string, wg *sync.WaitGroup) {
	// var mrnasWithProperties []MRNAWithProperties
	defer close(csvOutputChannel)
	defer wg.Done()

	//
	header := append(getMRNAFields(), GetFunctionNames(properties)...)
	csvOutputChannel <- header

	for {
		select {
		case mrna, ok := <-mrnaChannel:
			if !ok {
				return
			}

			var propertyValues []string
			// var mrnaProperties map[string]interface{}
			for _, property := range properties {
				propertyValue := property.computePropertyFn(mrna)
				mrna.properties[GetFunctionName(property.computePropertyFn)] = propertyValue
				if property.includeInOutputCSV {
					propertyValues = append(propertyValues, toString(propertyValue))
				}
				// mapValue, propertyValue := property(mrna, mrnaProperties)
				// mapKey := GetFunctionName(property)
				// mrnaProperties[mapKey] = mapValue
			}

			mrnaWithProperties := MRNAWithProperties{
				mrna:           mrna,
				propertyValues: propertyValues,
			}

			csvOutputChannel <- mrnaWithProperties.ToSlice()
		}
	}
}

func writeToCSV(outputFile string, osFlags int, csvOutputChannel chan []string, wg *sync.WaitGroup) {
	wFile, err := csv_helper.OpenFile(outputFile, osFlags)
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

func (d DatasetWithProperties) getHeaders() []string {

	mrna := d.mrnaWithProperties[0].mrna

	// Add fields of mrna to headers
	valueOfMRNA := reflect.ValueOf(mrna)
	typeOfMRNA := valueOfMRNA.Type()
	var headers []string
	for i := 0; i < valueOfMRNA.NumField(); i++ {
		headers = append(headers, typeOfMRNA.Field(i).Name)
	}

	// Add rest of fields of DatasetWithProperties to headers
	for i := 0; i < len(d.propertyNames); i++ {
		headers = append(headers, d.propertyNames[i])
	}

	return headers
}

func getMRNAFields() []string {
	mrna := MRNA{}

	// Add fields of mrna to headers
	valueOfMRNA := reflect.ValueOf(mrna)
	typeOfMRNA := valueOfMRNA.Type()
	var headers []string
	for i := 0; i < valueOfMRNA.NumField(); i++ {
		if valueOfMRNA.Field(i).CanInterface() {
			headers = append(headers, typeOfMRNA.Field(i).Name)
		}
	}
	return headers
}

func (mrna MRNA) ToSlice() []string {
	// Add fields of mrna to headers
	valueOfMRNA := reflect.ValueOf(mrna)
	var values []string
	for i := 0; i < valueOfMRNA.NumField(); i++ {
		if valueOfMRNA.Field(i).CanInterface() {
			values = append(values, fmt.Sprint(valueOfMRNA.Field(i).Interface()))
		}
	}
	return values
}

func (d MRNAWithProperties) ToSlice() []string {
	mrna := d.mrna

	values := mrna.ToSlice()
	// Add rest of fields of DatasetWithProperties to headers
	values = append(values, d.propertyValues...)

	return values
}

func GetFunctionName(fn interface{}) string {
	nameWithPackagePath := runtime.FuncForPC(reflect.ValueOf(fn).Pointer()).Name()
	nameSlice := strings.Split(nameWithPackagePath, ".")
	functionName := nameSlice[len(nameSlice)-1]
	return functionName
}

func GetFunctionNames(i []Property) []string {
	var names []string
	for _, property := range i {
		if property.includeInOutputCSV {
			names = append(names, GetFunctionName(property.computePropertyFn))
		}
	}
	// for _, fn := range i {
	// }
	return names
}
