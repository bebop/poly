package main

import (
	"embed"
	"fmt"
	"sync"

	"github.com/TimothyStiles/poly/rbs_calculator"
	"github.com/TimothyStiles/poly/rbs_calculator/csv_helper"
)

func main() {

}

func add16SrRNA(fs embed.FS, csvFile string, outputCSVFile string) {
	var wg sync.WaitGroup

	nbHeaderRows := 1
	headerRow := csv_helper.ReadHeader(fs, csvFile, nbHeaderRows)[0]
	headerRow = append(headerRow, "16S rRNA")

	// add to wait group for new csv rows channel
	wg.Add(1)
	csvRowsChannel := make(chan []string)
	go csv_helper.ReadCSV(fs, csvFile, csvRowsChannel, &wg, nbHeaderRows)

	// create a channel for csv output
	csvOutputChannel := make(chan []string)

	// write all header row to file
	csvOutputChannel <- headerRow

	wg.Add(1)
	organismColIdx := 2
	go doAdd16SrRNA(csvRowsChannel, organismColIdx, csvOutputChannel, &wg)

	wg.Add(1)
	go csv_helper.WriteToCSV(outputCSVFile, csv_helper.CREATE, csvOutputChannel, &wg)
}

func doAdd16SrRNA(csvInputChannel chan []string, organismColIdx int, csvOutputChannel chan []string, wg *sync.WaitGroup) {
	// var mrnasWithProperties []MRNAWithProperties
	defer close(csvOutputChannel)
	defer wg.Done()

	for {
		select {
		case csvRow, ok := <-csvInputChannel:
			if !ok {
				return
			}

			organism := csvRow[organismColIdx]

			if val, ok := rbs_calculator.Organism16SrRNAMap[organism]; ok {
				csvRow = append(csvRow, val)
			} else {
				panic(fmt.Errorf("do not have 16S rRNA for organism: %v", organism))
			}

			csvOutputChannel <- csvRow
		}
	}
}
