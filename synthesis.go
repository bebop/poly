package poly

import (
	"github.com/jmoiron/sqlx"
	"github.com/juliangruber/go-intersect"
	_ "github.com/mattn/go-sqlite3"
	"regexp"
	"strings"
	"sync"
)

/******************************************************************************
Dec, 9, 2020

Synthesis fixing stuff begins here.

Codon optimization is not good enough for synthetic DNA - there are certain
elements that synthesis doesn't really work for. For example, long homopolymers
will nearly always need removal in order to synthesize a gene efficiently.

Poly codon optimization + synthesis fixing is meant to be as efficient yet
accurate as possible.

# How do we do it?

From experience using Twist, there are a couple of major elements to fix:

1. Global GC content
2. ~50bp GC hotspots
3. Homopolymers
4. K-mer repeats
5. Various banned restriction enzyme sites

To implement this functionality, I imagine each major element will have a function
that takes in a sequence and returns 2 integers and either "GC", "AT", or "NA". The
2 integers are start and end positions of problematic sequences, and the "GC"/"AT"
are for if the sequence needs to be biased in a particular direction (we call this
the change type).

If there are two elements that have an overlap (for example, a high "GC" region
with a restriction enzyme site), fixing the outer element is tabled until the next
round of fixes. However, the nearest outer element change type will be inherited for
biased stochastic selection of a new codon. With the example of high "GC" and an
internal restriction enzyme site, the internal restriction enzyme will be changed
with an AT bias if possible.

Once all change sites are identified, changes are implemented across the entire
sequence. If any given change causes a problem which is *smaller* than the change,
that change is reverted, the entire context of the problematic region is noted, and
a change at a different position is attempted, until there is a change that fixes
the original problem and does not cause a problem which is *smaller* than that change.

The procedure above is repeated up the chain until all problems are resolved. If a problem
cannot be resolved, this is noted in the output, but it is not an error - "It is possible to
commit no mistakes and still lose. That is not a weakness. That is life."

[Removed notes about backtracking. It's hard and I don't really want to implement it right now]

Although there is lots of recursion and rechecking, I roughly estimate only 1~2 rounds will be
necessary for ~95% of sequences, which should be rather fast.

# Anything else?

This functionality should be fundamentally valuable to companies and labs. We will also
be able to improve the functionality of these fixes as we collect more data, or even
have company-by-company recommendations on things to fix. Hopefully by working in the
open, we can improve codon optimization and usage in a transparent way. One day,
this will enable us to be far better at expressing proteins at certain levels within
cells. That is the dream - working in the open will enable us to be far better at
engineering life.

Information wants to be free,

Keoni Gandall







This is Keoni coming back after a few months. Let me re-explain what this code will do.

1 First, we generate an in-memory SQLite DB. This database lets us do relational queries
  on important parts of our sequence that would be difficult to write in Golang
2 Next, we use the inputed functions to look for problematic regions (problematicRegion).
  Each problematicRegion has a list of suggested fixes (suggestedFixes) that, if fixed up,
  should resolve the problematicRegion.
3 Once we have the list of problematicRegions + suggestedFixes, we sort the suggestedFixes by
  the estimated amount of changes needed to fix that region.
3a The suggestedFix with the smallest quantity of changes necessary to fix the problem is then
   fixed in sequence space.
3b The specific codon that is changed to fix a problem is first sorted by if it has been used already.
   For example, if one codon has already been switched around, the system will avoid changing it
   again unless absolutely necessary.
3c Codons are then sorted to move towards optimal codons. For example, if one codon can be changed from
   a fairly rare codon to the most common codon for that amino acid, it will be chosen to be changed over
   a codon who is already fairly optimal.
4. All problematicRegions *without overlaps* are repeated with [GOTO 3] in order of quantity of changes
5. After all independent problematicRegions have been checked, repeat search with [GOTO 2]
6. If we end up with a sequence without any problematicRegions, we return the sequence!



******************************************************************************/

type DnaSuggestion struct {
	Start          int    `db:"start"`
	End            int    `db:"end"`
	Bias           string `db:"gcbias"`
	QuantityFixes  int    `db:"quantityfixes"`
	SuggestionType string `db:"suggestiontype"`
	Step           int    `db:"step"`
	Id             int    `db:"id"`
}

func FindBsaI(sequence string, c chan DnaSuggestion, wg *sync.WaitGroup) {
	re := regexp.MustCompile(`GGTCTC`)
	locs := re.FindAllStringIndex(sequence, -1)
	for _, loc := range locs {
		position := loc[0] / 3
		leftover := loc[0] % 3
		switch {
		case leftover == 0:
			c <- DnaSuggestion{position, (loc[1] / 3), "NA", 1, "BsaI removal", 0, 0}
		case leftover != 0:
			c <- DnaSuggestion{position, (loc[1] / 3) - 1, "NA", 1, "BsaI removal", 0, 0}
		}
	}
	wg.Done()
}

func findProblems(sequence string, problematicSequenceFuncs []func(string, chan DnaSuggestion, *sync.WaitGroup)) []DnaSuggestion {
	// Run functions to get suggestions
	suggestions := make(chan DnaSuggestion, 100)
	var wg sync.WaitGroup
	for _, f := range problematicSequenceFuncs {
		wg.Add(1)
		go f(sequence, suggestions, &wg)
	}
	wg.Wait()
	close(suggestions)

	var suggestionsList []DnaSuggestion
	for suggestion := range suggestions {
		suggestionsList = append(suggestionsList, suggestion)
	}
	return suggestionsList
}

func FixCds(sequence string, codontable CodonTable, problematicSequenceFuncs []func(string, chan DnaSuggestion, *sync.WaitGroup)) (string, error) {
	var db *sqlx.DB
	db = sqlx.MustConnect("sqlite3", ":memory:")
	createMemoryDbSql := `
	CREATE TABLE codon (
		codon TEXT PRIMARY KEY,
		aa TEXT
	);

	CREATE TABLE seq (
		pos INT PRIMARY KEY
	);

	CREATE TABLE history (
		pos INTEGER REFERENCES seq(pos),
		codon TEXT NOT NULL REFERENCES codon(codon),
		step INT
	);

	-- Weights are set on a per position basis for codon harmonization at a later point
	CREATE TABLE weights (
		pos INTEGER REFERENCES seq(pos),
		codon TEXT NOT NULL REFERENCES codon(codon),
		weight INTEGER
	);

	CREATE TABLE codonbias (
		fromcodon TEXT REFERENCES codon(codon),
		tocodon TEXT REFERENCES codon(codon),
		gcbias TEXT
	);

	CREATE TABLE suggestedfix (
		id INTEGER PRIMARY KEY AUTOINCREMENT,
		step INTEGER,
		start INTEGER REFERENCES seq(pos),
		end INTEGER REFERENCES seq(pos),
		gcbias TEXT,
		quantityfixes INTEGER,
		suggestiontype TEXT
	);
`
	db.MustExec(createMemoryDbSql)
	// Insert codons
	weightTable := make(map[string]int)
	codonInsert := `INSERT INTO codon(codon, aa) VALUES (?, ?)`
	for _, aminoAcid := range codontable.AminoAcids {
		for _, codon := range aminoAcid.Codons {
			db.MustExec(codonInsert, codon.Triplet, aminoAcid.Letter)
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
	for i := 0; i < len(sequence)-3; i = i + 3 {
		codon := sequence[i : i+3]
		db.MustExec(`INSERT INTO seq(pos) VALUES (?)`, pos)
		db.MustExec(`INSERT INTO history(pos, codon, step) VALUES (?, ?, 0)`, pos, codon)
		db.MustExec(`INSERT INTO weights(pos, codon, weight) VALUES (?,?,?)`, pos, codon, weightTable[codon])
		pos++
	}

	// For a maximum of 1000 iterations, see if we can do better
	for i := 1; i < 1000; i++ {
		suggestions := findProblems(sequence, problematicSequenceFuncs)
		// If there are no suggestions, break the iteration!
		if len(suggestions) == 0 {
			break
		}
		for suggestionIndex, suggestion := range suggestions {
			// First, let's insert the suggestions that we found using our problematicSequenceFuncs
			db.MustExec(`INSERT INTO suggestedfix(step, start, end, gcbias, quantityfixes, suggestiontype) VALUES (?, ?, ?, ?, ?, ?)`, i, suggestion.Start, suggestion.End, suggestion.Bias, suggestion.QuantityFixes, suggestion.SuggestionType)
			// Second, lets look for if there are any overlapping suggestions. This is equivalent to a spot where we can "kill two birds with one stone"
			// TODO: make this into a function so you can get overlapping overlaps (for when 3 problems are overlapped on top of each other)
			for secondSuggestionIndex, secondSuggestion := range suggestions {
				if suggestionIndex != secondSuggestionIndex {
					// makeRange is a helper function that will allow us to generate lists for intersection
					makeRange := func(min, max int) []int {
						a := make([]int, max-min+1)
						for i := range a {
							a[i] = min + i
						}
						return a
					}

					// The overlapping regions would need to have the same bias applied to them. So we check first if they have no preference or the same bias
					if (suggestion.Bias == "NA") || (secondSuggestion.Bias == "NA") || (suggestion.Bias == secondSuggestion.Bias) {
						// Each suggestion will have a start and end. If we take the intersection of the numbers between the start and the end of both suggestions,
						// we should be able to find a region where an edit would affect both suggested changes.
						overlap := intersect.Sorted(makeRange(suggestion.Start, suggestion.End), makeRange(secondSuggestion.Start, secondSuggestion.End)).([]int)
						if len(overlap) != 0 {
							// Each overlap that is successfully found will be inserted as a suggestedfix. First, we need to find its bias
							var overlapBias string
							if suggestion.Bias == "NA" {
								overlapBias = secondSuggestion.Bias
							} else {
								overlapBias = suggestion.Bias
							}

							// Each overlap will have a number of fixes that need to be hit to fix the outer objects
							var overlapQuantityFixes int
							if suggestion.QuantityFixes >= secondSuggestion.QuantityFixes {
								overlapQuantityFixes = suggestion.QuantityFixes
							} else {
								overlapQuantityFixes = secondSuggestion.QuantityFixes
							}

							// There can't be more desired fixes than there are positions
							if overlapQuantityFixes > len(overlap) {
								overlapQuantityFixes = len(overlap)
							}

							// Lets check if there is already an overlap at this position during this step
							var id int
							err := db.Get(&id, `SELECT id FROM suggestedfix WHERE step = ? AND start = ? AND end = ? AND suggestiontype = "overlap"`, i, overlap[0], overlap[len(overlap)-1])
							// If err is not nil, then we couldn't scan a single instance of the above query
							if err != nil {
								db.MustExec(`INSERT INTO suggestedfix(step, start, end, gcbias, quantityfixes, suggestiontype) VALUES (?, ?, ?, ?, ?, "overlap")`, i, overlap[0], overlap[len(overlap)-1], overlapBias, overlapQuantityFixes)
							}

						}
					}
				}
			}
		}

		// The following statements are the magic sauce that makes this all worthwhile.
		// Parameters: step, gcbias, start, end, quantityfix
		sqlFix := `INSERT INTO history
		            (codon,
		             pos,
		             step)
		SELECT t.codon,
		       t.pos,
		       ? AS step
		FROM   (SELECT cb.tocodon AS codon,
		               s.pos      AS pos
		        FROM   seq AS s
		               JOIN history AS h
		                 ON h.pos = s.pos
		               JOIN weights AS w
		                 ON w.pos = s.pos
		               JOIN codon AS c
		                 ON h.codon = c.codon
		               JOIN codonbias AS cb
		                 ON cb.fromcodon = c.codon
		        WHERE  cb.gcbias = ?
		               AND s.pos >= ?
		               AND s.pos <= ?
		               AND h.codon != cb.tocodon
		        ORDER  BY w.weight) AS t
		GROUP  BY t.pos
		LIMIT  ?; `

		// During this step, first see if there are any overlapping suggestedfixes
		var overlapSuggestions []DnaSuggestion
		db.Select(&overlapSuggestions, `SELECT * FROM suggestedfix WHERE suggestiontype = "overlap" AND step = ?`, i)
		// If there are overlapping suggestedfixes, fix those and iterate along
		if len(overlapSuggestions) > 0 {
			for _, overlapSuggestion := range overlapSuggestions {
				db.MustExec(sqlFix, i, overlapSuggestion.Bias, overlapSuggestion.Start, overlapSuggestion.End, overlapSuggestion.QuantityFixes)
			}
		} else { // If there aren't any overlapping suggestedfixes, fix the the remaining independent suggestedfixes
			independentSuggestions := []DnaSuggestion{}
			db.Select(&independentSuggestions, `SELECT * FROM suggestedfix WHERE step = ?`, i)

			for _, independentSuggestion := range independentSuggestions {
				db.MustExec(sqlFix, i, independentSuggestion.Bias, independentSuggestion.Start, independentSuggestion.End, independentSuggestion.QuantityFixes)
			}
		}
		var codons []string
		db.Select(&codons, `SELECT codon FROM (SELECT codon, pos FROM history ORDER BY step DESC) GROUP BY pos`)
		sequence = strings.Join(codons, "")
	}
	return sequence, nil
}
