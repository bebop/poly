package main

import (
	"fmt"
	"io/ioutil"
	"strings"
)

// this is a really dumb function that I had to write to get a comprehensive list of feature keys and qualifier keys for the genbank parser
// I thought it'd be good to leave here for posterity - Tim. Source Data : http://www.insdc.org/files/feature_table.html#3.3
func parseFeaturesList() {
	file, err := ioutil.ReadFile("gbfeatures.txt")
	if err != nil {
		fmt.Println(err)
	}
	splitFile := strings.Split(string(file), "\n")
	for _, line := range splitFile {
		splitLine := strings.Split(line, " ")
		if splitLine[0] == "Feature" {
			fmt.Println("\"" + strings.TrimSpace(strings.Join(splitLine[2:], " ")) + "\",")
		}

	}

}

// ^ ditto to above but for feature qualifiers : http://www.insdc.org/files/feature_table.html#3.3
func parseQualifiersList() {
	file, err := ioutil.ReadFile("gbqualifiers.txt")
	if err != nil {
		fmt.Println(err)
	}
	splitFile := strings.Split(string(file), "\n")
	for _, line := range splitFile {
		splitLine := strings.Split(line, " ")
		if strings.TrimSpace(splitLine[0]) == "Qualifier" {
			fmt.Println("\"" + strings.TrimSpace(strings.Join(splitLine[1:], " ")) + "\",")
		}

	}

}
