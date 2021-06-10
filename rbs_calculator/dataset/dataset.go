package dataset

import (
	"encoding/csv"
	"fmt"
	"io"
	"log"
	"os"
	"reflect"
	"runtime"
	"strconv"
	"strings"

	"github.com/TimothyStiles/poly"
	"github.com/TimothyStiles/poly/rbs_calculator/csv_helper"
)

/**********************
Step 1: Add folded structrue to dataset

Helper functions to add the folded structure of a sequence to a dataset
***********************/

func AddFoldedStructure(dataset, datasetOutputFile string, mRNAColNumber int) error {

	// csvOutputChannel := make(chan []string)

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

/**********************
Step 2: Add variables of model to visualize and explore

Functions to add variables to visualize and explore to the dataset
***********************/

// MRNA contains all the 'ground-truth' information of a mRNA sequence
type MRNA struct {
	Sequence              string
	Structure             string
	FivePrimeUTR          string
	ProteinCodingSequence string
	Output                float64
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
func createDataset(csvFile string, sequenceColNum, fivePrimeUTRColNum, cdsColNum, structureColNum, outputColNum int, mrnaChannel chan MRNA) {
	defer close(mrnaChannel)

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
		var output float64
		if output, err = strconv.ParseFloat(csvRow[outputColNum], 64); err != nil {
			panic(err)
		}
		// sequence, structure := cleanRNA(csvRow[sequenceColNum]), csvRow[structureColNum]
		// _, energyContributions, err := mfe.MinimumFreeEnergy(sequence, structure, mfe.DefaultTemperature)
		// if err != nil {
		// 	panic(err)
		// }
		mRNA := MRNA{
			Sequence:              cleanRNA(csvRow[sequenceColNum]),
			Structure:             csvRow[structureColNum],
			Output:                output,
			FivePrimeUTR:          cleanRNA(csvRow[fivePrimeUTRColNum]),
			ProteinCodingSequence: cleanRNA(csvRow[cdsColNum]),
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
func streamDataset(csvFile string, sequenceColNum, fivePrimeUTRColNum, cdsColNum, structureColNum, outputColNum, batchSize, batchNum int) (Dataset, bool) {
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
			Sequence:              cleanRNA(csvRow[sequenceColNum]),
			Structure:             csvRow[structureColNum],
			Output:                output,
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

func computeProperties(dataset Dataset, properties [](func(MRNA) string)) DatasetWithProperties {
	var mrnasWithProperties []MRNAWithProperties

	for _, mrna := range dataset {
		var propertyValues []string
		for _, property := range properties {
			propertyValues = append(propertyValues, property(mrna))
		}

		mrnaWithProperties := MRNAWithProperties{
			mrna:           mrna,
			propertyValues: propertyValues,
		}
		mrnasWithProperties = append(mrnasWithProperties, mrnaWithProperties)
	}

	return DatasetWithProperties{
		mrnaWithProperties: mrnasWithProperties,
		propertyNames:      GetFunctionNames(properties),
	}
}

func computeProperties2(mrnaChannel chan MRNA, properties [](func(MRNA) string), csvOutputChannel chan []string) {
	// var mrnasWithProperties []MRNAWithProperties
	defer close(csvOutputChannel)

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
			for _, property := range properties {
				propertyValues = append(propertyValues, property(mrna))
			}

			mrnaWithProperties := MRNAWithProperties{
				mrna:           mrna,
				propertyValues: propertyValues,
			}

			// fmt.Println(mrnaWithProperties)
			csvOutputChannel <- mrnaWithProperties.ToSlice()
			// computedPropChannel <- mrnaWithProperties
		}
	}

	// return DatasetWithProperties{
	// 	mrnaWithProperties: mrnasWithProperties,
	// 	propertyNames:      GetFunctionNames(properties),
	// }
}

func writeToCSV(outputFile string, csvOutputChannel chan []string) {
	// defer wg.Done()
	wFile, err := csv_helper.CreateFile(outputFile)
	if err != nil {
		panic(err)
	}
	defer wFile.Close()
	csvWriter := csv.NewWriter(wFile)
	defer csvWriter.Flush()

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

func writeDatasetWithPropertiesToCSV(dataset DatasetWithProperties, outputFile string) {
	wFile, err := csv_helper.CreateFile(outputFile)
	if err != nil {
		panic(err)
	}
	defer wFile.Close()
	csvWriter := csv.NewWriter(wFile)
	defer csvWriter.Flush()

	headers := dataset.getHeaders()
	if err := csvWriter.Write(headers); err != nil {
		panic(err)
	}

	for _, mrnaWithProperties := range dataset.mrnaWithProperties {
		values := mrnaWithProperties.ToSlice()
		if err := csvWriter.Write(values); err != nil {
			panic(err)
		}
	}
}

func (d DatasetWithProperties) getHeaders() []string {
	// valueOfDatasetWithProperties := reflect.ValueOf(d)
	// typeOfDatasetWithProperties := valueOfDatasetWithProperties.Type()

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
		headers = append(headers, typeOfMRNA.Field(i).Name)
	}
	return headers
}

func (d MRNAWithProperties) ToSlice() []string {
	mrna := d.mrna

	// Add fields of mrna to headers
	valueOfMRNA := reflect.ValueOf(mrna)
	var values []string
	for i := 0; i < valueOfMRNA.NumField(); i++ {
		if valueOfMRNA.Field(i).CanInterface() {
			values = append(values, fmt.Sprint(valueOfMRNA.Field(i).Interface()))
		}
	}

	// Add rest of fields of DatasetWithProperties to headers
	values = append(values, d.propertyValues...)

	fmt.Println(values)
	return values
}

func GetFunctionName(fn interface{}) string {
	nameWithPackagePath := runtime.FuncForPC(reflect.ValueOf(fn).Pointer()).Name()
	nameSlice := strings.Split(nameWithPackagePath, ".")
	functionName := nameSlice[len(nameSlice)-1]
	return functionName
}

func GetFunctionNames(i [](func(MRNA) string)) []string {
	var names []string
	for _, fn := range i {
		names = append(names, GetFunctionName(fn))
	}
	return names
}

/***********************
Functions to compute properties of mRNAs
***********************/

// The default properties to compute values of and add to DatasetWithProperties
var DefaultProperies [](func(MRNA) string) = [](func(MRNA) string){
	// same,
}

func Reverse(s string) string {
	runes := []rune(s)
	for i, j := 0, len(runes)-1; i < j; i, j = i+1, j-1 {
		runes[i], runes[j] = runes[j], runes[i]
	}
	return string(runes)
}
