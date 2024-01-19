/*
Package align is a package for aligning (comparing) DNA, RNA, and protein sequences.

Biology is fickle and full of quirks that make it hard to do even the most basic of tasks
which we would normally take for granted when working with other kinds of data.

Comparing two biological sequences to see if they're roughly equivalent is one of those tasks.

Essentially two almost identical sequences with almost identical functionality can contain
small insertions or deletions that shift the entire string such that a meaningful comparison via
hamming distance or levenshtein distance becomes impossible.

For example:

Timothy Stiles
||||||| ||||||
Timothy Stiles

is an easy match with hamming or levenshtein distance.

However, say we introduce a new character to the beginning of the sequence.

Timothy Stiles
xxxxxxxxxxxxxxxx
A Timothy Stiles

Now our edit distance via levenshtein is maximized at 16 and we wouldn't
be able to tell that semantically these strings are almost identical.

This frame shifting seen above is incredibly common within biological sequences and alignment
algorithms are designed in part to deal with these shifts so that when we compare two sequences
like the two below we can get a more useful edit distance.

GAAAAAAT
GAA----T

As of writing this package includes the two most basic algorithms for alignment,
Needleman-Wunsch and Smith-Waterman. Needleman-Wunsch is used when you are
looking for global alignment between two full-length sequences, while
Smith-Waterman is better at smaller sequences with local similarities and
handling sequences with long non-homologous regions. BLAST, on the other hand,
takes advantage of more heuristic techniques to speed up alignment, and is better
at finding similar sequences in large database, sacrificing precision for faster
results.

Both are "dynamic programming algorithms" which is a fancy 1980's term for they use
matrices. If you're familiar with kernel operations, linear filters, or whatever term
ML researchers are using nowadays for, "slide a window over a matrix and determine that
entry's values using its neighbor's values", then this should be pretty easy to grok.

If not these algorithms essentially compare every character in one sequence with another
sequence and create an edit distance along with human readable string to show gaps like the
previous example.

I'm not really an expert on alignment so if you want to learn more about this class of algorithms
wikipedia has a decent overview.

https://en.wikipedia.org/wiki/Sequence_alignment

Even if I may not know the answer to your alignment questions please ask and I'll do my best
to help!

TTFN,
Tim
*/
package align

import (
	"github.com/bebop/poly/search/align/matrix"
)

// Scoring is a struct that holds the scoring matrix for match, mismatch, and gap penalties.
type Scoring struct {
	SubstitutionMatrix *matrix.SubstitutionMatrix
	GapPenalty         int
}

// NewScoring returns a new Scoring struct with default values for DNA.
func NewScoring(substitutionMatrix *matrix.SubstitutionMatrix, gapPenalty int) (Scoring, error) {
	if substitutionMatrix == nil {
		substitutionMatrix = matrix.Default
	}
	return Scoring{
		SubstitutionMatrix: substitutionMatrix,
		GapPenalty:         gapPenalty,
	}, nil
}

func (s Scoring) Score(a, b byte) (int, error) {
	matchScore, err := s.SubstitutionMatrix.Score(string(a), string(b))
	if err != nil {
		return 0, err
	}
	return matchScore, nil
}

// NeedlemanWunsch performs global alignment between two strings using the Needleman-Wunsch algorithm.
// It returns the final score and the optimal alignments of the two strings in O(nm) time and O(nm) space.
// https://en.wikipedia.org/wiki/Needleman-Wunsch_algorithm
func NeedlemanWunsch(stringA string, stringB string, scoring Scoring) (int, string, string, error) {
	// Get the M and N dimensions of the matrix. The M x N matrix is standard linear algebra notation.
	// But I added columns and rows to the variable name to make it more clear what the dimensions are.
	columnLengthM, rowLengthN := len(stringA), len(stringB)

	// Initialize the matrix columns
	matrix := make([][]int, columnLengthM+1) // matrix is a 2D slice of ints.

	// Initialize the matrix rows.
	for columnM := range matrix {
		matrix[columnM] = make([]int, rowLengthN+1)
	}

	// Fill in the first column with gap penalties.
	for columnM := 1; columnM <= columnLengthM; columnM++ {
		matrix[columnM][0] = matrix[columnM-1][0] + scoring.GapPenalty
	}

	// Fill in the first row with gap penalties.
	for rowN := 1; rowN <= rowLengthN; rowN++ {
		matrix[0][rowN] = matrix[0][rowN-1] + scoring.GapPenalty
	}

	// Fill in the rest of the matrix.
	for columnM := 1; columnM <= columnLengthM; columnM++ {
		for rowN := 1; rowN <= rowLengthN; rowN++ {
			// Calculate the scores for scoring.Match/mismatch and gap.
			var score, err = scoring.Score(stringA[columnM-1], stringB[rowN-1])
			if err != nil {
				return 0, "", "", err
			}
			matrix[columnM][rowN] = max(
				matrix[columnM-1][rowN-1]+score,
				max(matrix[columnM-1][rowN]+scoring.GapPenalty, matrix[columnM][rowN-1]+scoring.GapPenalty),
			)
		}
	}

	// Traceback to find the optimal alignment.
	var alignA, alignB []rune
	columnM, rowN := columnLengthM, rowLengthN
	for columnM > 0 && rowN > 0 {
		var matchScore, err = scoring.Score(stringA[columnM-1], stringB[rowN-1])
		if err != nil {
			return 0, "", "", err
		}
		if matrix[columnM][rowN] == matrix[columnM-1][rowN-1]+matchScore {
			alignA = append(alignA, rune(stringA[columnM-1]))
			alignB = append(alignB, rune(stringB[rowN-1]))
			columnM--
			rowN--
		} else if matrix[columnM][rowN] == matrix[columnM-1][rowN]+scoring.GapPenalty {
			alignA = append(alignA, rune(stringA[columnM-1]))
			alignB = append(alignB, '-')
			columnM--
		} else {
			alignA = append(alignA, '-')
			alignB = append(alignB, rune(stringB[rowN-1]))
			rowN--
		}
	}

	// Reverse the alignments to get the optimal alignment.
	alignA = reverseRuneArray(alignA)
	alignB = reverseRuneArray(alignB)
	return matrix[columnLengthM][rowLengthN], string(alignA), string(alignB), nil
}

// SmithWaterman performs local alignment between two strings using the Smith-Waterman algorithm.
// It returns the max score and optimal local alignments between two strings alignments of the two strings in O(nm) time and O(nm) space.
// https://en.wikipedia.org/wiki/Smith-Waterman_algorithm
func SmithWaterman(stringA string, stringB string, scoring Scoring) (int, string, string, error) {
	columnLengthM, rowLengthN := len(stringA), len(stringB)

	// Initialize the alignment matrix
	matrix := make([][]int, columnLengthM+1)
	for columnM := 0; columnM <= columnLengthM; columnM++ {
		matrix[columnM] = make([]int, rowLengthN+1)
	}

	// Initialize variables to keep track of the maximum score and its position
	maxScore := 0
	maxScoreRow := 0
	maxScoreCol := 0

	// Fill the alignment matrix
	for columnM := 1; columnM <= columnLengthM; columnM++ {
		for rowN := 1; rowN <= rowLengthN; rowN++ {
			var matchScore, err = scoring.Score(stringA[columnM-1], stringB[rowN-1])
			if err != nil {
				return 0, "", "", err
			}
			diagScore := matrix[columnM-1][rowN-1] + matchScore
			upScore := matrix[columnM-1][rowN] + scoring.GapPenalty
			leftScore := matrix[columnM][rowN-1] + scoring.GapPenalty
			matrix[columnM][rowN] = max(0, max(diagScore, max(upScore, leftScore)))

			if matrix[columnM][rowN] > maxScore {
				maxScore = matrix[columnM][rowN]
				maxScoreRow = columnM
				maxScoreCol = rowN
			}
		}
	}

	// Traceback to construct the aligned strings
	alignA := ""
	alignB := ""
	columnM := maxScoreRow
	rowN := maxScoreCol
	for matrix[columnM][rowN] > 0 {
		var matchScore, err = scoring.Score(stringA[columnM-1], stringB[rowN-1])
		if err != nil {
			return 0, "", "", err
		}
		if matrix[columnM][rowN] == matrix[columnM-1][rowN-1]+matchScore {
			alignA = string(stringA[columnM-1]) + alignA
			alignB = string(stringB[rowN-1]) + alignB
			columnM--
			rowN--
		} else if matrix[columnM][rowN] == matrix[columnM-1][rowN]+scoring.GapPenalty {
			alignA = string(stringA[columnM-1]) + alignA
			alignB = "-" + alignB
			columnM--
		} else if matrix[columnM][rowN] == matrix[columnM][rowN-1]+scoring.GapPenalty {
			alignA = "-" + alignA
			alignB = string(stringB[rowN-1]) + alignB
			rowN--
		}
	}

	return maxScore, alignA, alignB, nil
}

func reverseRuneArray(runes []rune) []rune { // wasn't able to find a built-in reverse function for runes
	length := len(runes)
	for index := 0; index < length/2; index++ {
		reverseIndex := length - index - 1
		runes[index], runes[reverseIndex] = runes[reverseIndex], runes[index]
	}
	return runes
}

func max(a, b int) int { // funny enough Go's built-in max only handles floats?
	if a > b {
		return a
	}
	return b
}
