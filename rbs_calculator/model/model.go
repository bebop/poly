package model

import (
	"embed"
	"fmt"
	"reflect"
	"runtime"
	"strconv"
	"strings"
	"sync"

	"github.com/TimothyStiles/poly/rbs_calculator/csv_helper"
)

/*****************************************************************************
Process to develop a model for the RBS Calculator

Current research provides intuition as to the process behind a ribosome binding
to a mRNA, but there is no concensus on an actual model of this interaction.

To develop a model, we need examine different primary and secondary structure
features and figure out their relationship to the end result (translation
initiation rate).

To create the model, we follow a Data -> Model Theory, Variables & Assumptions
-> Visualization & Quantitative Relationship Formalization -> Final Model
Equation workflow.


Data:
We need a high quality dataset that includes the mRNA sequence, its structure,
the organism the sequence was inserted into, the protein expressed by the CDS
and the amount of protein synthesized.
The mRNA sequence must be split into its 5' untranslated region (5' UTR) and its
protein coding sequence (CDS). This data is available in the `dataset` folder.
Please read the readme file in the folder for more information about the
datasets.


Model Theory, Variables & Assumptions:
After we have the data, we need to create a theoretical model of the
interactions of a ribosome with a mRNA strand. Once we have a model, we need to
write the code to come up with the variables required by the model. This part
of the process is specific to the theories people have in mind so nothing more
can be generalized about this part. If you'd like to create your own model,
create a subpackage under this subpackage (for example `model/salis_lab_v2_1`)
and describe your model and write the code to compute properties for your model
in a `properties.go` file in your created subpackage. Please have a look at the
`properties.go` file in the `model/salis_lab_v2_1` subpackage for an example.

Once the code is written to compute your desired properties, we have to run
the ground-truth dataset against your code to generate properties for all the
data points in the dataset. This package exports the function
`ComputeProperties` which takes each data point in the 'ground-truth' dataset,
computes the properties you've specified, and outputs a csv file in a directory
named `dataset_with_properties` in your subpackage.

Please see the `properties_test.go` file in the `model/salis_lab_v2_1`
subpackage to understand how to create a dataset with your computed properties.
You can see an example of the dataset file after properties have been computed
(`model/salis_lab_v2_1/dataset_with_properties/train.csv`).


Visualization & Quantitative Relationship Formalization:

After the above step, we have a dataset with properties, but we don't yet
quantitatively understand the relationship between a computed property and
the ground-truth protein levels.

Thus, to understand these relationships, we need to visualize the relationship
between a computed property and the ground-truth protein levels.

Once we can see a (linear, quadratic, cubic) relationship between a variable
and the outcome, we can figure out the quantitative relationship through
curve-fitting apps like MATLAB and add this quantitative relationship to our
model.

Unfortunately, guidance to do this step is not yet included in this repo as
we've based our first calculator on research papers from The Salis Lab which
have included the quantitative relationship between properties and the
resulting protein levels. We hope to include this information and a workflow
for this step in the future if / when the model of the RBS calculator is
improved.


Final Model Equation:

Although the previous step gives us the relationships between individual
variables and the outcome, it could be the case that adding all the variables
together leads to worse prediction of results.

Thus, as a final step we run regression analysis that tries out all the possible
combinations of variables in the model to see if removing some terms could
actually improve the predicted results. The least set of variables with the
highest correlation with actual results are then used in the model.

Unfortunately, as in the step above, guidance to do this is not included in
this repo as we've not had to do this yet.


Implementation details to create a model:
To create a model, create a subpackage in the `model` package with a
`properties.go` file. This file must include the following:
* Create a subpackage in the `model` package with a `properties.go` file.
* The `properties.go` file should include a `PropertiesToCompute` variable of
  type `[]Property`.
	The value returned by each `ComputePropertyFn` in `PropertiesToCompute` is
	added to the `Properties` map (the key in the map is the name of the
	`ComputePropertyFn`) of the `RibosomeBindingSite` struct. Thus, the order of
	the properties in `PropertiesToCompute` is important and properties later in
	the `PropertiesToCompute` slice can use the values of properties earlier in
	the slice by accessing the relevant value using the `ComputePropertyFn` name
	of the earlier property.
* Create a copy of `properties_test.go` from the `salis_lab_v2_1` subpackage
	into your own subpackage. You will have to update the name of the package on
	the top of the file. Then compute the properties you wrote by `cd`ing into
	the directory of your subpackage and running
	`go test -timeout 0 -run ^TestComputeProperties`.
* You will now have a csv file with your computed properties for each data
	point in the source dataset. At this point, you will have to curve-fit and
	quantitatively figure out the relationship between the values of each computed
	property and their relationship to the actual protein mean and protein std
	values. (As mentioned above, we don't provide guidance on this step currently.
	Excel and MATLAB are useful tools for this. PRs are very welcome.)
* Update your `properties.go` file to include the quantitative relationships
	from the previous step.
* Finally, create a function named `TranslationInitiationRate` of type
	`func(RibosomeBindingSite) float64` that returns the translation initiation
	rate for a ribosome binding site.
* Update `rbs_calculator.go` in the `rbs_calculator` package to use your model
	by updating the `rbs_model` import at the top of the file.

To compute statistics for your model, read the documentation of `stats.py`.

*****************************************************************************/

//go:embed datasets/*
var embeddedDatasetsDirectory embed.FS
var datasetsDirectory = "datasets"

// RibosomeBindingSite contains all the 'ground-truth' information of a mRNA
// sequence. This struct is used to create a model for the rbs calculator (
// please read the description of the `model.ComputeProperties` func for more
// information), and is used for inference (please see the code in the
// `TranslationInitiationRate` func in the `rbs_calculator` package for more
// information).
type RibosomeBindingSite struct {
	// FivePrimeUTR is the untranslated region of the mRNA sequence on the five
	// prime end
	FivePrimeUTR string
	// ProteinCodingSequence is the sequence after the 5' UTR that encodes a
	// protein
	ProteinCodingSequence string
	Temperature           float64
	// RibosomalRNA is the 16S rRNA of the organism the mRNA strand is inserted
	// into
	RibosomalRNA string
	// Properties is a map that contains the computed properties
	Properties map[string]interface{}
	// OtherInformation contains information that we'd like to output to the
	// dataset with properties, but is not used to compute properties of a
	// binding site
	OtherInformation []string
	// TranslationInitiationRate is the computed translation initiation rate for
	// this binding site. Please note that this field is only set when infering
	// in the `TranslationInitiationRate` func in the `rbs_calculator` package
	TranslationInitiationRate float64
}

// Property contains information about the properties that need to be computed
// for a binding site (and whether or not to include the computed property in
// the csv outputted to the `dataset_with_properties` directory after calling
// the `ComputeProperties` func)
type Property struct {
	ComputePropertyFn  func(*RibosomeBindingSite) interface{}
	IncludeInOutputCSV bool
}

// CleanRNA converts DNA to RNA and makes sequence upper case
func CleanRNA(rna string) string {
	rna = strings.TrimSpace(rna)
	rna = strings.ToUpper(rna)
	rna = strings.ReplaceAll(rna, "T", "U")
	return rna
}

// ComputeProperties creates a `RibosomeBindingSite` struct for each data point
// in the input dataset, computes the desired properties (the
// `propertiesToCompute` argument), and outputs each `RibosomeBindingSite`
// struct after the properties have been computed to the
// `dataset_with_properties` directory in the package this function is called
// from.
//
// The computed properties that are included in the output csv file is
// controlled by the `IncludeInOutputCSV` field of the `Property` struct.
//
// The source dataset may contain information that we'd like to include the
// output csv which may not be part of the `RibosomeBindingSite` struct.
// The `otherInformationMap` is an argument which specifies a map of column
// number to column header name which will be included in the outputted csv
// file.
func ComputeProperties(datasetName string, propertiesToCompute []Property,
	fivePrimeUTRColIdx, cdsColIdx, tempColIdx, ribosomalRNAColIdx int,
	otherInformationMap map[int]string) {
	dataset := datasetsDirectory + "/" + datasetName + ".csv"

	// get the required information from `otherInformationMap`
	var otherInformationColIdxs []int
	var otherInformationColHeaders []string
	for colNum, colHeaderName := range otherInformationMap {
		otherInformationColIdxs = append(otherInformationColIdxs, colNum)
		otherInformationColHeaders = append(otherInformationColHeaders, colHeaderName)
	}

	// Create a wait group that will only 'release' when all go subroutines
	// call `wg.Done()`
	var wg sync.WaitGroup

	// Step 1: populate a channel with `RibosomeBindingSite` structs from each row
	// of the dataset
	rbsChannel := make(chan *RibosomeBindingSite)
	// add to the wait group for the `populateMRNAChannelFromDataset` subroutine
	wg.Add(1)
	go populateRBSChannelFromDataset(dataset, rbsChannel, &wg, fivePrimeUTRColIdx, cdsColIdx, tempColIdx, ribosomalRNAColIdx, otherInformationColIdxs...)

	// Step 2: For each struct in `rbsChannel`, compute properties for it (based
	// on `propertiesToCompute`), convert the struct to `[]string` and push it
	// to `csvOutputChannel`
	csvOutputChannel := make(chan []string)
	// add to the wait group for the `computeProperties` subroutine
	wg.Add(1)
	go computeRBSProperties(rbsChannel, propertiesToCompute, csvOutputChannel, &wg, otherInformationColHeaders...)

	// Step 3: Output data from `csvOutputChannel` to a CSV file
	datasetOutputFile := "./dataset_with_properties/" + datasetName + ".csv"
	// add the current unix timestamp to the file name
	datasetOutputFile = csv_helper.FileNameWithUNIXTimestamp(datasetOutputFile)
	// add to the wait group for the `csv_helper.WriteToCSV` subroutine
	wg.Add(1)
	go csv_helper.WriteToCSV(datasetOutputFile, csv_helper.CREATE, csvOutputChannel, &wg)

	// wait till all the subroutines call `wg.Done()`
	wg.Wait()
}

// populateRBSChannelFromDataset creates a `RibosomeBindingSite` struct for each
// row in the `csvFile` and adds it to `rbsChannel`
func populateRBSChannelFromDataset(csvFile string, rbsChannel chan *RibosomeBindingSite,
	wg *sync.WaitGroup, fivePrimeUTRColIdx, cdsColIdx, tempColIdx, RibosomalRNAColIdx int,
	otherInformationColIdxs ...int) {
	// close the rbs channel when this func returns
	defer close(rbsChannel)
	// call `wg.Done()` when this func returns
	defer wg.Done()

	// a channel that will hold rows of the CSV that are scanned
	csvRowsChannel := make(chan []string)
	// add to the wait group for the `csv_helper.ReadCSV` subroutine
	wg.Add(1)
	go csv_helper.ReadCSV(embeddedDatasetsDirectory, csvFile, csvRowsChannel, wg, 1)

	// weird Go syntax to wait for a message from a channel
	for {
		select {
		case csvRow, ok := <-csvRowsChannel:
			if !ok {
				// occurs when `csvRowsChannel` is closed by `csv_helper.ReadCSV`
				// so return this func as well as we don't have any more input
				return
			}

			if strings.TrimSpace(csvRow[fivePrimeUTRColIdx]) == "" && strings.TrimSpace(csvRow[cdsColIdx]) == "" {
				panic("dataset row doesn't contain both five prime UTR and coding sequence")
			}

			// parse temperature into `float64`
			var temp float64
			var err error
			if temp, err = strconv.ParseFloat(csvRow[tempColIdx], 64); err != nil {
				panic(err)
			}

			// get the other information which don't contribute to the
			// `RibosomeBindingSite` struct, but we want to keep in the output csv
			var otherInformation []string
			for _, otherInformationColNum := range otherInformationColIdxs {
				otherInformation = append(otherInformation, csvRow[otherInformationColNum])
			}

			// finally, create a `RibosomeBindingSite` struct with the required info
			rbs := RibosomeBindingSite{
				FivePrimeUTR:          CleanRNA(csvRow[fivePrimeUTRColIdx]),
				ProteinCodingSequence: CleanRNA(csvRow[cdsColIdx]),
				Temperature:           temp,
				RibosomalRNA:          CleanRNA(csvRow[RibosomalRNAColIdx]),
				Properties:            make(map[string]interface{}),
				OtherInformation:      otherInformation,
			}

			// add the struct to the channel for use in the `computeRBSProperties` subroutine
			rbsChannel <- &rbs

		}
	}
}

// computeRBSProperties computes properties for each `RibosomeBindingSite`
// struct in `rbsChannel`, convert the struct to `[]string`, and sends the
// `[]string` to `csvOutputChannel`
func computeRBSProperties(rbsChannel chan *RibosomeBindingSite,
	properties []Property, csvOutputChannel chan []string,
	wg *sync.WaitGroup, otherInformationColHeaders ...string) {
	// close the csv output channel when this func returns
	defer close(csvOutputChannel)
	// call `wg.Done()` when this func returns
	defer wg.Done()

	// add the required `RibsomeBindingSite` sturct fields to the header
	header := ribosomeBindingSiteStructFields()
	// add the other information column headers to the output csv file header
	header = append(header, otherInformationColHeaders...)
	// add the names of the property functions (that need to be included in the
	// output csv) to the header
	header = append(header, functionNames(properties)...)
	// write the header to the csv
	csvOutputChannel <- header

	// weird Go syntax to wait for a message from a channel
	for {
		select {
		case rbs, ok := <-rbsChannel:
			if !ok {
				// occurs when `rbsChannel` is closed by `populateRBSChannelFromDataset`
				// so return this func as well as we don't have any more rbs structs
				// to process
				return
			}

			// compute properties for the rbs struct and get the list of computed
			// property values that need to be added to the output csv
			propertyValues := rbs.ComputeProperties(properties)

			// add the required values of the rbs struct to the output
			// note: this should be in the same order as the header row above
			var output []string = rbs.fieldValues()
			output = append(output, rbs.OtherInformation...)
			output = append(output, propertyValues...)

			csvOutputChannel <- output
		}
	}
}

// rbsStructFieldsToExclude contains fields that we don't want to include
// This variable is used in the `ribosomeBindingSiteStructFields` and
// `RibosomeBindingSite.fieldValues` funcs
var rbsStructFieldsToExclude map[string]bool = map[string]bool{
	"OtherInformation":          true,
	"Properties":                true,
	"TranslationInitiationRate": true,
}

// ribosomeBindingSiteStructFields returns the fields of the
// `RibosomeBindingSite` struct excluding fields from the
// `rbsStructFieldsToExclude` map
func ribosomeBindingSiteStructFields() (fieldNames []string) {
	rbs := RibosomeBindingSite{}
	valueOfRBSStruct := reflect.ValueOf(rbs)
	typeOfRBSStruct := valueOfRBSStruct.Type()

	for i := 0; i < valueOfRBSStruct.NumField(); i++ {
		if valueOfRBSStruct.Field(i).CanInterface() {
			fieldName := typeOfRBSStruct.Field(i).Name

			// exclude fields that are in the `rbsStructFieldsToExclude` map
			if rbsStructFieldsToExclude[fieldName] {
				continue
			}
			fieldNames = append(fieldNames, typeOfRBSStruct.Field(i).Name)
		}
	}
	return
}

// fieldValues returns the values of the fields of `rbs` excluding fields from
// the `rbsStructFieldsToExclude` map
func (rbs RibosomeBindingSite) fieldValues() (fieldValues []string) {
	valueOfRBSStruct := reflect.ValueOf(rbs)
	typeOfRBSStruct := valueOfRBSStruct.Type()

	for i := 0; i < valueOfRBSStruct.NumField(); i++ {
		if valueOfRBSStruct.Field(i).CanInterface() {
			fieldName := typeOfRBSStruct.Field(i).Name

			// exclude fields that are in the `rbsStructFieldsToExclude` map
			if rbsStructFieldsToExclude[fieldName] {
				continue
			}
			fieldValue := toString(valueOfRBSStruct.Field(i).Interface())
			fieldValues = append(fieldValues, fieldValue)
		}
	}
	return
}

// functionName returns a string of the name of the function `fn`
func functionName(fn interface{}) string {
	// get the full name of the function with package path
	// `nameWithPackagePath` will be "<package>.<function name>"
	fnNameWithPackagePath := runtime.FuncForPC(reflect.ValueOf(fn).Pointer()).Name()

	// get only the function name
	fnNameWithPackagePathSlice := strings.Split(fnNameWithPackagePath, ".")
	functionName := fnNameWithPackagePathSlice[len(fnNameWithPackagePathSlice)-1]

	return functionName
}

// functionNames returns the names of the function from a list of `Property`
func functionNames(properties []Property) (names []string) {
	for _, property := range properties {
		if property.IncludeInOutputCSV {
			names = append(names, functionName(property.ComputePropertyFn))
		}
	}

	return
}

// ComputeProperties computes the properties specified by the `properties`
// argument and returns a list of the computed property values if the property's
// `IncludeInOutputCSV` field is set to true
func (rbs *RibosomeBindingSite) ComputeProperties(properties []Property) (propertyValues []string) {
	for _, property := range properties {
		propertyValue := property.ComputePropertyFn(rbs)

		// add the computed property value to the `Properties` map
		rbs.Properties[functionName(property.ComputePropertyFn)] = propertyValue

		// add the computed property value to the list that will be returned (if
		// required)
		if property.IncludeInOutputCSV {
			propertyValues = append(propertyValues, toString(propertyValue))
		}
	}

	return
}

// toString returns the string value of `i`
func toString(i interface{}) string {
	return fmt.Sprint(i)
}
