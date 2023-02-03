package align

// Scoring matrix for match and mismatch.
var match = 1
var mismatch = -1

// Penalty for gap.
var gapPenalty = -1

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}

// NeedlemanWunsch performs global alignment between two strings using the Needleman-Wunsch algorithm.
// It returns the final score and the optimal alignments of the two strings in O(nm) time and O(nm) space.
func NeedlemanWunsch(stringA, stringB string) (int, string, string) {

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
		matrix[columnM][0] = matrix[columnM-1][0] + gapPenalty
	}

	// Fill in the first row with gap penalties.
	for rowN := 1; rowN <= rowLengthN; rowN++ {
		matrix[0][rowN] = matrix[0][rowN-1] + gapPenalty
	}

	// Fill in the rest of the matrix.
	for columnM := 1; columnM <= columnLengthM; columnM++ {
		for rowN := 1; rowN <= rowLengthN; rowN++ {
			// Calculate the scores for match/mismatch and gap.
			var score int
			if stringA[columnM-1] == stringB[rowN-1] {
				score = match
			} else {
				score = mismatch
			}
			matrix[columnM][rowN] = max(
				matrix[columnM-1][rowN-1]+score,
				max(matrix[columnM-1][rowN]+gapPenalty, matrix[columnM][rowN-1]+gapPenalty),
			)
		}
	}

	// Traceback to find the optimal alignment.
	var alignA, alignB []rune
	columnM, rowN := columnLengthM, rowLengthN
	for columnM > 0 && rowN > 0 {
		if stringA[columnM-1] == stringB[rowN-1] {
			alignA = append(alignA, rune(stringA[columnM-1]))
			alignB = append(alignB, rune(stringB[rowN-1]))
			columnM--
			rowN--
		} else if matrix[columnM][rowN] == matrix[columnM-1][rowN-1]+mismatch {
			alignA = append(alignA, rune(stringA[columnM-1]))
			alignB = append(alignB, rune(stringB[rowN-1]))
			columnM--
			rowN--
		} else if matrix[columnM][rowN] == matrix[columnM-1][rowN]+gapPenalty {
			alignA = append(alignA, rune(stringA[columnM-1]))
			alignB = append(alignB, '-')
			columnM--
		} else {
			alignA = append(alignA, '-')
			alignB = append(alignB, rune(stringB[rowN-1]))
			rowN--
		}
	}

	// Add remaining characters from A.
	for columnM > 0 {
		alignA = append(alignA, rune(stringA[columnM-1]))
		alignB = append(alignB, '-')
		columnM--
	}
	// Add remaining characters from B.
	for rowN > 0 {
		alignA = append(alignA, '-')
		alignB = append(alignB, rune(stringB[rowN-1]))
		rowN--
	}

	// Reverse the alignments to get the optimal alignment.
	alignA = reverse(alignA)
	alignB = reverse(alignB)
	return matrix[columnLengthM][rowLengthN], string(alignA), string(alignB)
}

func reverse(r []rune) []rune {
	for i := 0; i < len(r)/2; i++ {
		j := len(r) - i - 1
		r[i], r[j] = r[j], r[i]
	}
	return r
}
