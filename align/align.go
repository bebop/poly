package align

// Scoring is a struct that holds the scoring matrix for match, mismatch, and gap penalties.
type Scoring struct {
	Match      int
	Mismatch   int
	GapPenalty int
}

// NewScoring returns a new Scoring struct with default values.
func NewScoring() Scoring {
	return Scoring{
		Match:      1,
		Mismatch:   -1,
		GapPenalty: -1,
	}
}

// NeedlemanWunsch performs global alignment between two strings using the Needleman-Wunsch algorithm.
// It returns the final score and the optimal alignments of the two strings in O(nm) time and O(nm) space.
// https://en.wikipedia.org/wiki/Needleman-Wunsch_algorithm
func NeedlemanWunsch(stringA string, stringB string, scoring Scoring) (int, string, string) {

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
			var score int = score(stringA[columnM-1], stringB[rowN-1], scoring)

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
		if matrix[columnM][rowN] == matrix[columnM-1][rowN-1]+score(stringA[columnM-1], stringB[rowN-1], scoring) {
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
	return matrix[columnLengthM][rowLengthN], string(alignA), string(alignB)
}

// SmithWaterman
func SmithWaterman(stringA string, stringB string, scoring Scoring) (int, string, string) {

	columnLengthM, rowLengthN := len(stringA), len(stringB)

	// Initialize the alignment matrix
	matrix := make([][]int, columnLengthM+1)
	for i := 0; i <= columnLengthM; i++ {
		matrix[i] = make([]int, rowLengthN+1)
	}

	// Initialize variables to keep track of the maximum score and its position
	maxScore := 0
	maxScoreRow := 0
	maxScoreCol := 0

	// Fill the alignment matrix
	for i := 1; i <= columnLengthM; i++ {
		for j := 1; j <= rowLengthN; j++ {
			diagScore := matrix[i-1][j-1] + score(stringA[i-1], stringB[j-1], scoring)
			upScore := matrix[i-1][j] + scoring.GapPenalty
			leftScore := matrix[i][j-1] + scoring.GapPenalty
			matrix[i][j] = max(0, max(diagScore, max(upScore, leftScore)))

			if matrix[i][j] > maxScore {
				maxScore = matrix[i][j]
				maxScoreRow = i
				maxScoreCol = j
			}
		}
	}

	// Traceback to construct the aligned strings
	alignA := ""
	alignB := ""
	i := maxScoreRow
	j := maxScoreCol
	for matrix[i][j] > 0 {
		if matrix[i][j] == matrix[i-1][j-1]+score(stringA[i-1], stringB[j-1], scoring) {
			alignA = string(stringA[i-1]) + alignA
			alignB = string(stringB[j-1]) + alignB
			i--
			j--
		} else if matrix[i][j] == matrix[i-1][j]+scoring.GapPenalty {
			alignA = string(stringA[i-1]) + alignA
			alignB = "-" + alignB
			i--
		} else if matrix[i][j] == matrix[i][j-1]+scoring.GapPenalty {
			alignA = "-" + alignA
			alignB = string(stringB[j-1]) + alignB
			j--
		} else {
			// Should not happen
			panic("Unexpected case in traceback")
		}
	}

	return maxScore, alignA, alignB
}

func score(a, b byte, scoring Scoring) int {
	if a == b {
		return scoring.Match
	} else {
		return scoring.Mismatch
	}
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
