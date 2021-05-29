package rbs_calculator

// import (
// 	"encoding/csv"
// 	"io"
// 	"io/fs"
// 	"os"
// 	"path/filepath"
// 	"regexp"
// 	"strconv"
// )

// A file that contains helper functions needed to compare poly's RBS calculator
// with other implementations of RBS calculators and with data from experiments

var defaultTranslationInitiationRateHeader = "r"

// processCSV reads a CSV (skipping the first header line), and passes mRNA sequence
// at column number `mRNAColNumber` to `RibosomeBindingSite()`. The output is then
// written to the last available column in the CSV under the header
// func processCSV(mRNAColNumber int, translationInitiationRateHeader string) {
// 	translationInitiationRateHeader := defaultTranslationInitiationRateHeader
// 	rRNA := defaultRRNA

// }
