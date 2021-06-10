package dataset

import (
	"encoding/csv"
	"fmt"
	"io"
	"os"
	"reflect"
	"regexp"
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

func AddFoldedStructure(dataset, datasetOutputFile string, csvStartRowIdx int, mRNAColNumber int) {
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
		csvOutputChannel <- header
	}

	for i := 1; i < csvStartRowIdx; i++ {
		csvr.Read()
	}

	csvRow, err := csvr.Read()
	for err == nil {
		mRNA := csvRow[mRNAColNumber]
		structure, _ := poly.LinearFold(mRNA, poly.DefaultBeamSize)
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
	Sequence              string
	Structure             string
	FivePrimeUTR          string
	ProteinCodingSequence string
	Output                float64
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
func createDataset(csvFile string, sequenceColNum, fivePrimeUTRColNum, cdsColNum, structureColNum, outputColNum int, mrnaChannel chan MRNA, wg *sync.WaitGroup) {
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
				propertyValue := property.computeFn(mrna)
				// p := reflect.ValueOf(mrna).MethodByName(property)
				// if !p.IsValid() {
				// 	panic(property)
				// }
				// reflectValues := p.Call([]reflect.Value{})
				// propertyValue := reflectValues[0].Interface()
				mrna.properties[GetFunctionName(property.computeFn)] = propertyValue
				if property.includeInOutput {
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
	defer wg.Done()
	wFile, err := csv_helper.OpenFile(outputFile, osFlags)
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

func GetFunctionNames(i []Property) []string {
	var names []string
	for _, property := range i {
		if property.includeInOutput {
			names = append(names, GetFunctionName(property.computeFn))
		}
	}
	// for _, fn := range i {
	// }
	return names
}

/***********************
Functions to compute properties of mRNAs
***********************/

// The default properties to compute values of and add to DatasetWithProperties
// var DefaultProperies [](func(MRNA, map[string]interface{}) (interface{}, string)) = [](func(MRNA, map[string]interface{}) (interface{}, string)){
// 	// same,
// 	lenFivePrimeUTR,
// 	sdSequence,
// }
type Property struct {
	computeFn       func(MRNA) interface{}
	includeInOutput bool
}

// var DefaultProperies map[string]bool = map[string]bool{
// 	"LenFivePrimeUTR": true,
// 	"SDSequence":      true,
// 	"LenSDSequence":   true,
// }
var DefaultProperies []Property = []Property{
	Property{MRNA.LenFivePrimeUTR, true},
	Property{MRNA.SDSequence, false},
	Property{MRNA.LenSDSequence, true},
}

func toString(i interface{}) string {
	return fmt.Sprint(i)
}

func (mrna MRNA) LenFivePrimeUTR() interface{} {
	lenFivePrimeUTR := len(mrna.FivePrimeUTR)
	return lenFivePrimeUTR
}

func (mrna MRNA) SDSequence() interface{} {
	fivePrimeUTR := mrna.FivePrimeUTR
	// fivePrimeUTR = Reverse(fivePrimeUTR)
	// antiSDSequence := "ACCUCCUUA"
	smallestSDSequence := "GAG"
	// 18 is because SD sequence can end max nine nucelotides from the start codon
	// and is at max 9 nucleotides long. We remove 2 because GAG is two away from
	// the start of the SD sequence
	smallestSDSequenceStartPos := 2
	maxDistFromStartCodonForSDSequece := 18 - smallestSDSequenceStartPos
	minDistForSDSequence := 3

	// The sequences of the mRNA to search for the SD sequence in
	sdMRNASearchSequence := fivePrimeUTR[len(fivePrimeUTR)-maxDistFromStartCodonForSDSequece : len(fivePrimeUTR)-minDistForSDSequence]
	smallestSDSequenceRegex := regexp.MustCompile(smallestSDSequence)

	matches := smallestSDSequenceRegex.FindAllStringIndex(sdMRNASearchSequence, -1)
	if matches == nil {
		// No SD sequence in 5' untranslated region
		return ""
	} else {
		maxLenSDSequence := 0
		var bestSDSequence string
		for _, match := range matches {
			smallestSDSequenceStartIdx, smallestSDSequenceEndIdx := len(fivePrimeUTR)-maxDistFromStartCodonForSDSequece+match[0], len(fivePrimeUTR)-maxDistFromStartCodonForSDSequece+match[1]
			lenSDSequence, sdSequence := getFullSDSequence(fivePrimeUTR, smallestSDSequence, smallestSDSequenceStartPos, smallestSDSequenceStartIdx, smallestSDSequenceEndIdx)
			if lenSDSequence >= maxLenSDSequence {
				maxLenSDSequence = lenSDSequence
				bestSDSequence = sdSequence
			}
		}
		return bestSDSequence
	}
}

func (mrna MRNA) LenSDSequence() interface{} {
	// panic(len(mrna.properties["SDSequence"].(string)))
	sdSequence := mrna.properties["SDSequence"]
	if sdSequence == nil {
		return 0
	} else {
		return len(mrna.properties["SDSequence"].(string))
	}
}

func getFullSDSequence(fivePrimeUTR, smallestSDSequence string, smallestSDSequenceStartPos, smallestSDSequenceStartIdx, smallestSDSequenceEndIdx int) (int, string) {
	smallestSDSequenceEndPos := smallestSDSequenceStartPos + len(smallestSDSequence)
	minSDSequenceIdx, maxSDSequenceIdx := smallestSDSequenceStartIdx, smallestSDSequenceEndIdx
	fullSdSequence := "ACCUCCUUA"
	for i := smallestSDSequenceStartPos - 1; i >= 0; i-- {
		if minSDSequenceIdx-1 == 0 {
			panic("here")
		}
		if fullSdSequence[i] == fivePrimeUTR[minSDSequenceIdx-1] {
			minSDSequenceIdx--
		} else {
			break
		}
	}

	for i := smallestSDSequenceEndPos + 1; i < len(fullSdSequence); i++ {
		if fullSdSequence[i] == fivePrimeUTR[maxSDSequenceIdx+1] {
			maxSDSequenceIdx++
		} else {
			break
		}
	}

	lenSDSequence := maxSDSequenceIdx - minSDSequenceIdx
	return lenSDSequence, fivePrimeUTR[minSDSequenceIdx : maxSDSequenceIdx+1]
}

/********************
Helper functions used to compute properties
*********************/

func Reverse(s string) string {
	runes := []rune(s)
	for i, j := 0, len(runes)-1; i < j; i, j = i+1, j-1 {
		runes[i], runes[j] = runes[j], runes[i]
	}
	return string(runes)
}
