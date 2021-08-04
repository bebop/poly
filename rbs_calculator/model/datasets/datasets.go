package model

import (
	"embed"
	"fmt"
	"sync"

	"github.com/TimothyStiles/poly/rbs_calculator"
	"github.com/TimothyStiles/poly/rbs_calculator/csv_helper"
)

/***********************************************
Functions for modifying datasets in the `datasets` directory.
Currenlty only includes a func `add16SrRNA` to add the 16S
rRNA of an organism to a dataset.
***********************************************/

func add16SrRNA(fs embed.FS, csvFile string, outputCSVFile string, organismColIdx, rrnaColIdx int) {
	var wg sync.WaitGroup

	nbHeaderRows := 1
	headerRow := csv_helper.ReadHeader(fs, csvFile, nbHeaderRows)[0]
	headerRow[rrnaColIdx] = "16S rRNA"

	// create a channel for csv output
	csvOutputChannel := make(chan []string, 1)

	// write all header row to file
	csvOutputChannel <- headerRow

	// add to wait group for new csv rows channel
	wg.Add(1)
	csvRowsChannel := make(chan []string)
	go csv_helper.ReadCSV(fs, csvFile, csvRowsChannel, &wg, nbHeaderRows)

	wg.Add(1)
	go doAdd16SrRNA(csvRowsChannel, organismColIdx, rrnaColIdx, csvOutputChannel, &wg)

	wg.Add(1)
	go csv_helper.WriteToCSV(outputCSVFile, csv_helper.CREATE, csvOutputChannel, &wg)

	wg.Wait()
}

func doAdd16SrRNA(csvInputChannel chan []string, organismColIdx, rrnaColIdx int, csvOutputChannel chan []string, wg *sync.WaitGroup) {
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
				csvRow[rrnaColIdx] = val
			} else {
				panic(fmt.Errorf("do not have 16S rRNA for organism: %v", organism))
			}

			csvOutputChannel <- csvRow
		}
	}
}
