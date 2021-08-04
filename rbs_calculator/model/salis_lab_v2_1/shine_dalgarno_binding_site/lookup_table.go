package shine_dalgarno_binding_site

import (
	"embed"
	"encoding/csv"
	"io"
	"strconv"

	"github.com/TimothyStiles/poly/rbs_calculator/csv_helper"
)

//go:embed data/*
var embeddedLookupTableDataDirectory embed.FS
var rbsLookupDataDirectory = "data"

type ShineDalgarnoBindingSite struct {
	MessengerRNA, RibosomalRNA                   string
	MessengerRNAStructure, RibosomalRNAStructure string
	ShineDalgarnoFreeEnergy                      float64
	MRNAFivePrimeIdx, MRNAThreePrimeIdx          int
	RRNAFivePrimeIdx, RRNAThreePrimeIdx          int
	AlignedSpacing                               int
	AlignedSpacingFreeEnergy                     float64
}

func parseCSVIntoLookupTable(fs embed.FS, file string, lookupTable *map[string]map[string]map[float64]ShineDalgarnoBindingSite) error {
	f, err := fs.Open(file)
	if err != nil {
		return err
	}
	defer f.Close()

	csvr := csv.NewReader(f)

	// consume header line
	csvr.Read()

	rRNAColNum, mRNAColNum, sdFreeEnergyColNum, mrnaStructureColNum, rrnaStructureColNum, mrnaFivePrimeColNum, mrnaThreePrimeColNum, rrnaFivePrimeColNum, rrnaThreePrimeColNum, temperatureColNum := 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

	for {
		row, err := csvr.Read()
		if err != nil {
			if err == io.EOF {
				err = nil
			}
			return err
		}

		ribosomalRNA, mRNA := row[rRNAColNum], row[mRNAColNum]

		if (*lookupTable)[ribosomalRNA] == nil {
			(*lookupTable)[ribosomalRNA] = make(map[string]map[float64]ShineDalgarnoBindingSite)
		}

		if (*lookupTable)[ribosomalRNA][mRNA] == nil {
			(*lookupTable)[ribosomalRNA][mRNA] = make(map[float64]ShineDalgarnoBindingSite)
		}

		var bindingSiteFreeEnergy float64 = parseFloat64(row[sdFreeEnergyColNum])

		ribosomeBindingSite := ShineDalgarnoBindingSite{
			RibosomalRNA: ribosomalRNA, MessengerRNA: mRNA,
			MessengerRNAStructure: row[mrnaStructureColNum], RibosomalRNAStructure: row[rrnaStructureColNum],
			MRNAFivePrimeIdx:        parseInt(row[mrnaFivePrimeColNum]),
			MRNAThreePrimeIdx:       parseInt(row[mrnaThreePrimeColNum]),
			RRNAFivePrimeIdx:        parseInt(row[rrnaFivePrimeColNum]),
			RRNAThreePrimeIdx:       parseInt(row[rrnaThreePrimeColNum]),
			ShineDalgarnoFreeEnergy: bindingSiteFreeEnergy,
		}

		var temperature float64 = parseFloat64(row[temperatureColNum])
		(*lookupTable)[ribosomalRNA][mRNA][temperature] = ribosomeBindingSite
	}
}

func parseInt(token string) int {
	valueInt64, err := strconv.ParseInt(token, 10, 0)
	if err != nil {
		panic(err)
	}
	return int(valueInt64)
}

func parseFloat64(token string) (ret float64) {
	ret, err := strconv.ParseFloat(token, 64)
	if err != nil {
		panic(err)
	}
	return
}

func LookupTable() map[string]map[string]map[float64]ShineDalgarnoBindingSite {
	lookupTable := make(map[string]map[string]map[float64]ShineDalgarnoBindingSite)
	csvFiles := csv_helper.CSVFilesFromEmbeddedFS(embeddedLookupTableDataDirectory, rbsLookupDataDirectory)

	for _, csv := range csvFiles {
		err := parseCSVIntoLookupTable(embeddedLookupTableDataDirectory, csv, &lookupTable)
		if err != nil {
			panic(err)
			// return nil, err
		}
	}
	return lookupTable
}
