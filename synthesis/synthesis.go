/*
Package synthesis fixes synthetic DNA molecules in preparation for synthesis.

Many synthesis companies have restrictions on the DNA they can synthesize. This
synthesis fixer takes advantage of synonymous codons in protein coding
sequences (CDS) to remove problematic sequences that either users don't want
(like restriction enzymes sites) or that would cause DNA synthesis companies to
reject a synthesis project.

This synthesis fixer is meant to cover the majority of use cases for DNA
fixing. It is not intended to cover all possible use cases, since the majority
of DNA design does not actually have these edge cases.

FixCds does not guarantee that all requested features will be removed. If you
have use case that FixCds cannot properly fix, please put an issue in the poly
github.
*/
package synthesis

import (
	"errors"
	"fmt"
	"regexp"
	"strings"
	"sync"

	"github.com/TimothyStiles/poly/checks"
	"github.com/TimothyStiles/poly/transform"
	"github.com/TimothyStiles/poly/transform/codon"
	"github.com/jmoiron/sqlx"
	_ "modernc.org/sqlite" // imports CGO-less sqlite
)

// DnaSuggestion is a suggestion of a fixer, generated by a
// problematicSequenceFunc. Bias must be `NA`, `GC`, or `AT`, with `NA`
// representing a neutral skew.
type DnaSuggestion struct {
	Start          int    `db:"start"`
	End            int    `db:"end"`
	Bias           string `db:"gcbias"`
	QuantityFixes  int    `db:"quantityfixes"`
	SuggestionType string `db:"suggestiontype"`
}

// Change is a change to a given DNA sequence. A list of changes is given as
// the output of FixCds.
type Change struct {
	Position int    `db:"position"`
	Step     int    `db:"step"`
	From     string `db:"codonfrom"`
	To       string `db:"codonto"`
	Reason   string `db:"reason"`
}

type dbDnaSuggestion struct {
	Start          int    `db:"start"`
	End            int    `db:"end"`
	Bias           string `db:"gcbias"`
	QuantityFixes  int    `db:"quantityfixes"`
	SuggestionType string `db:"suggestiontype"`
	Step           int    `db:"step"`
	ID             int    `db:"id"`
}

// RemoveSequence is a generator for a problematicSequenceFuncs for specific
// sequences.
func RemoveSequence(sequencesToRemove []string, reason string) func(string, chan DnaSuggestion, *sync.WaitGroup) {
	return func(sequence string, c chan DnaSuggestion, waitgroup *sync.WaitGroup) {
		var sequencesToRemoveForReverse []string
		for _, seq := range sequencesToRemove {
			sequencesToRemoveForReverse = []string{seq, transform.ReverseComplement(seq)}
			for _, site := range sequencesToRemoveForReverse {
				re := regexp.MustCompile(site)
				locations := re.FindAllStringIndex(sequence, -1)
				for _, location := range locations {
					codonLength := 3
					position := location[0] / codonLength
					leftover := location[0] % codonLength
					switch {
					case leftover == 0:
						c <- DnaSuggestion{position, (location[1] / codonLength), "NA", 1, reason}
					case leftover != 0:
						c <- DnaSuggestion{position, (location[1] / codonLength) - 1, "NA", 1, reason}
					}
				}
			}
		}
		waitgroup.Done()
	}
}

// RemoveRepeat is a generator to make a problematicSequenceFunc for repeats.
func RemoveRepeat(repeatLen int) func(string, chan DnaSuggestion, *sync.WaitGroup) {
	return func(sequence string, c chan DnaSuggestion, waitgroup *sync.WaitGroup) {
		// Get a kmer list
		kmers := make(map[string]bool)
		for sequencePosition := 0; sequencePosition < len(sequence)-repeatLen; sequencePosition++ {
			_, alreadyFoundForward := kmers[sequence[sequencePosition:sequencePosition+repeatLen]]
			_, alreadyFoundReverse := kmers[transform.ReverseComplement(sequence[sequencePosition:sequencePosition+repeatLen])]
			kmers[sequence[sequencePosition:sequencePosition+repeatLen]] = true
			if alreadyFoundForward || alreadyFoundReverse {
				codonLength := 3
				position := sequencePosition / codonLength
				leftover := sequencePosition % codonLength
				endPosition := (sequencePosition + repeatLen) / codonLength
				switch {
				case leftover == 0:
					c <- DnaSuggestion{position, endPosition, "NA", 1, "Repeat sequence"}
				case leftover != 0:
					c <- DnaSuggestion{position, endPosition - 1, "NA", 1, "Repeat sequence"}
				}
				sequencePosition = sequencePosition + leftover
			}
		}
		waitgroup.Done()
	}
}

// GcContentFixer is a generator to increase or decrease the overall GcContent
// of a CDS. GcContent is defined as the percentage of guanine and cytosine
// base pairs in comparison to adenine and thymine base pairs. Usually, you
// want the range to be somewhere around 50%, with a decent upperBound being
// 80% GC and a decent lowerBound being 20%.
func GcContentFixer(upperBound, lowerBound float64) func(string, chan DnaSuggestion, *sync.WaitGroup) {
	return func(sequence string, c chan DnaSuggestion, waitgroup *sync.WaitGroup) {
		gcContent := checks.GcContent(sequence)
		var numberOfChanges int
		if gcContent > upperBound {
			numberOfChanges = int((gcContent-upperBound)*float64(len(sequence))) + 1
			c <- DnaSuggestion{0, (len(sequence) / 3) - 1, "AT", numberOfChanges, "GcContent too high"}
		}
		if gcContent < lowerBound {
			numberOfChanges = int((lowerBound-gcContent)*float64(len(sequence))) + 1
			c <- DnaSuggestion{0, (len(sequence) / 3) - 1, "GC", numberOfChanges, "GcContent too low"}
		}
		waitgroup.Done()
	}

}

// findProblems is a helper function in FixCDS that concurrently runs each
// sequence check and returns a list of all the suggested changes.
func findProblems(sequence string, problematicSequenceFuncs []func(string, chan DnaSuggestion, *sync.WaitGroup)) []DnaSuggestion {
	// Run functions to get suggestions
	suggestions := make(chan DnaSuggestion, 100)
	var waitgroup sync.WaitGroup
	for _, function := range problematicSequenceFuncs {
		waitgroup.Add(1)
		go function(sequence, suggestions, &waitgroup)
	}
	waitgroup.Wait()
	close(suggestions)

	var suggestionsList []DnaSuggestion
	for suggestion := range suggestions {
		suggestionsList = append(suggestionsList, suggestion)
	}
	return suggestionsList
}

// FixCds fixes a CDS given the CDS sequence, a codon table, a list of
// functions to solve for, and a number of iterations to attempt fixing.
// 10 iterations is a reasonable default for fixIterations.
func FixCds(sqlitePath string, sequence string, codontable codon.Table, problematicSequenceFuncs []func(string, chan DnaSuggestion, *sync.WaitGroup), fixIterations int) (string, []Change, error) {

	if len(sequence)%3 != 0 {
		return "", []Change{}, errors.New("this sequence isn't a complete CDS, please try to use a CDS without interrupted codons")
	}

	db := sqlx.MustConnect("sqlite", sqlitePath)
	createMemoryDbSQL := `
	PRAGMA foreign_keys = ON;
	CREATE TABLE codon (
		codon TEXT PRIMARY KEY,
		aa TEXT
	);

	CREATE TABLE seq (
		pos INT PRIMARY KEY
	);

	-- Weights are set on a per position basis for codon harmonization at a later point
	CREATE TABLE weights (
		pos INTEGER REFERENCES seq(pos),
		codon TEXT NOT NULL,
		weight INTEGER,
		FOREIGN KEY(codon) REFERENCES codon(codon)
	);

	CREATE TABLE codonbias (
		gcbias TEXT CHECK(gcbias IN ('NA', 'GC', 'AT')),
		fromcodon TEXT NOT NULL,
		tocodon TEXT NOT NULL,
		FOREIGN KEY(fromcodon) REFERENCES codon(codon),
		FOREIGN KEY(tocodon) REFERENCES codon(codon)
	);

	CREATE TABLE suggestedfix (
		id INTEGER PRIMARY KEY AUTOINCREMENT,
		step INTEGER,
		start INT NOT NULL,
		end INT NOT NULL,
		gcbias TEXT,
		quantityfixes INTEGER,
		suggestiontype TEXT,
		FOREIGN KEY(start) REFERENCES seq(pos),
		FOREIGN KEY(end) REFERENCES seq(pos)
	);

	CREATE TABLE history (
                pos INTEGER,
                codon TEXT NOT NULL,
                step INT,
                suggestedfix INT,
                FOREIGN KEY(codon) REFERENCES codon(codon),
                FOREIGN KEY(suggestedfix) REFERENCES suggestedfix(id),
		FOREIGN KEY(pos) REFERENCES seq(pos)
        );
`
	db.MustExec(createMemoryDbSQL)
	// Insert codons
	weightTable := make(map[string]int)
	codonInsert := `INSERT INTO codon(codon, aa) VALUES (?, ?)`
	// First just insert the codons
	for _, aminoAcid := range codontable.AminoAcids {
		for _, codon := range aminoAcid.Codons {
			db.MustExec(codonInsert, codon.Triplet, aminoAcid.Letter)
		}
	}
	// Then, add in GC biases
	for _, aminoAcid := range codontable.AminoAcids {
		for _, codon := range aminoAcid.Codons {
			weightTable[codon.Triplet] = codon.Weight

			codonBias := strings.Count(codon.Triplet, "G") + strings.Count(codon.Triplet, "C")
			for _, toCodon := range aminoAcid.Codons {
				if codon.Triplet != toCodon.Triplet {
					toCodonBias := strings.Count(toCodon.Triplet, "G") + strings.Count(toCodon.Triplet, "C")
					switch {
					case codonBias == toCodonBias:
						db.MustExec(`INSERT INTO codonbias(fromcodon, tocodon, gcbias) VALUES (?, ?, ?)`, codon.Triplet, toCodon.Triplet, "NA")
					case codonBias > toCodonBias:
						db.MustExec(`INSERT INTO codonbias(fromcodon, tocodon, gcbias) VALUES (?, ?, ?)`, codon.Triplet, toCodon.Triplet, "AT")
					case codonBias < toCodonBias:
						db.MustExec(`INSERT INTO codonbias(fromcodon, tocodon, gcbias) VALUES (?, ?, ?)`, codon.Triplet, toCodon.Triplet, "GC")
					}
				}
			}
		}
	}

	// Insert seq and history
	pos := 0
	for codonPosition := 0; codonPosition < len(sequence); codonPosition = codonPosition + 3 {
		codon := sequence[codonPosition : codonPosition+3]
		db.MustExec(`INSERT INTO seq(pos) VALUES (?)`, pos)
		db.MustExec(`INSERT INTO history(pos, codon, step) VALUES (?, ?, 0)`, pos, codon)
		pos++
	}

	for key, value := range weightTable {
		db.MustExec(`INSERT INTO weights(codon, weight) VALUES (?,?)`, key, value)
	}

	// For a maximum of 100 iterations, see if we can do better. Usually sequences will be solved within 1-3 rounds,
	// so 100 just effectively acts as the max cap for iterations. Once you get to 100, you pretty much know that
	// we cannot fix the sequence.
	for i := 1; i < fixIterations; i++ {
		suggestions := findProblems(sequence, problematicSequenceFuncs)
		// If there are no suggestions, break the iteration!
		if len(suggestions) == 0 {
			// Add a historical log of changes
			var changes []Change
			// This SQL will basically select positions + steps and organize
			// the codons by firstly their positions, and then by the step. We
			// want codons from the last fix iteration. It then selects the
			// codons at each position that has been changed and returns what it
			// has been changed to and what it has been changed from (this requires
			// a subquery)
			getChangesSql := `SELECT h.pos             AS position,
					       h.step            AS step,
					       (SELECT codon
					        FROM   history
					        WHERE  pos = h.pos
					               AND step = h.step - 1
					        LIMIT  1)        AS codonfrom,
					       h.codon           AS codonto,
					       sf.suggestiontype AS reason
					FROM   history AS h
					       JOIN suggestedfix AS sf
					         ON sf.id = h.suggestedfix
					WHERE  h.suggestedfix IS NOT NULL
					ORDER  BY h.step,
					          h.pos 
					`
			_ = db.Select(&changes, getChangesSql)
			return sequence, changes, nil
		}
		for _, suggestion := range suggestions { // if you want to add overlaps, add suggestionIndex
			// Check for a valid bias direction
			var validBias bool
			switch suggestion.Bias {
			case "NA", "GC", "AT":
				validBias = true
			}
			if !validBias {
				return sequence, []Change{}, fmt.Errorf("Invalid bias. Expected NA, GC, or AT, got %s", suggestion.Bias)
			}
			// First, let's insert the suggestions that we found using our problematicSequenceFuncs
			_, err := db.Exec(`INSERT INTO suggestedfix(step, start, end, gcbias, quantityfixes, suggestiontype) VALUES (?, ?, ?, ?, ?, ?)`, i, suggestion.Start, suggestion.End, suggestion.Bias, suggestion.QuantityFixes, suggestion.SuggestionType)
			if err != nil {
				return sequence, []Change{}, err
			}
		}

		// The following statements are the magic sauce that makes this all worthwhile.
		// Parameters: step, gcbias, start, end, quantityfix
		sqlFix1 := `INSERT INTO history
		            (codon,
		             pos,
		             step,
			     suggestedfix)
		SELECT t.codon,
		       t.pos,
		       ? AS step,
		       ? AS suggestedfix
		FROM   (SELECT cb.tocodon AS codon,
		               s.pos      AS pos
		        FROM   seq AS s
		               JOIN history AS h
		                 ON h.pos = s.pos
		               JOIN codon AS c
		                 ON h.codon = c.codon
		               JOIN codonbias AS cb
		                 ON cb.fromcodon = c.codon
			       JOIN 'weights' AS w
				 ON w.codon = cb.tocodon
		        WHERE `
		sqlFix2 := ` s.pos >= ?
		               AND s.pos <= ?
		               AND h.codon != cb.tocodon
		        ORDER  BY h.step DESC, w.weight DESC) AS t
		GROUP  BY t.pos
		LIMIT  ?; `

		var independentSuggestions []dbDnaSuggestion
		_ = db.Select(&independentSuggestions, `SELECT * FROM suggestedfix WHERE step = ?`, i)

		for _, independentSuggestion := range independentSuggestions {
			switch independentSuggestion.Bias {
			case "NA":
				db.MustExec(sqlFix1+sqlFix2, i, independentSuggestion.ID, independentSuggestion.Start, independentSuggestion.End, independentSuggestion.QuantityFixes)
			case "GC":
				db.MustExec(sqlFix1+`cb.gcbias = 'GC' AND `+sqlFix2, i, independentSuggestion.ID, independentSuggestion.Start, independentSuggestion.End, independentSuggestion.QuantityFixes)
			case "AT":
				db.MustExec(sqlFix1+`cb.gcbias = 'AT' AND `+sqlFix2, i, independentSuggestion.ID, independentSuggestion.Start, independentSuggestion.End, independentSuggestion.QuantityFixes)
			}
		}
		var codons []string
		_ = db.Select(&codons, `SELECT codon FROM (SELECT codon, pos FROM history ORDER BY step DESC) GROUP BY pos`)
		sequence = strings.Join(codons, "")
	}

	var changes []Change
	_ = db.Select(&changes, `SELECT h.pos AS position, h.step AS step, (SELECT codon FROM history WHERE pos = h.pos AND step = h.step-1 LIMIT 1) AS codonfrom, h.codon AS codonto, sf.suggestiontype AS reason FROM history AS h JOIN suggestedfix AS sf ON sf.id = h.suggestedfix WHERE h.suggestedfix IS NOT NULL ORDER BY h.step, h.pos`)

	return sequence, changes, errors.New("Could not find a solution to sequence space")
}

// FixCdsSimple is FixCds with some defaults for normal usage, including
// removing of homopolymers, removing any repeat larger than 18 base pairs, and
// fixing if a CDS's gc content is above 80% or below 20%
func FixCdsSimple(sequence string, codontable codon.Table, sequencesToRemove []string) (string, []Change, error) {
	var functions []func(string, chan DnaSuggestion, *sync.WaitGroup)
	// Remove homopolymers
	functions = append(functions, RemoveSequence([]string{"AAAAAAAA", "GGGGGGGG"}, "Homopolymers"))

	// Remove user defined sequences
	functions = append(functions, RemoveSequence(sequencesToRemove, "Removal requested by user"))

	// Remove repeats
	functions = append(functions, RemoveRepeat(18))

	// Ensure normal GC range
	functions = append(functions, GcContentFixer(0.80, 0.20))

	return FixCds(":memory:", sequence, codontable, functions, 100)
}
