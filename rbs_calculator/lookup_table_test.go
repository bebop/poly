package rbs_calculator

import (
	"fmt"
	"testing"
)

func ExampleLookup() {
	dG_rRNA_mRNA, err := Initalize()
	if err != nil {
		fmt.Printf("Failed to initialize lookup table: %s", err)
		return
	}
	fmt.Println(dG_rRNA_mRNA["ACCTCCTTA"]["CAAGGAGGGTG"])
	// Output: -11.9
}

// Test consists of looking up one randomly selected value from each csv file
func TestLookup(t *testing.T) {
	dG_rRNA_mRNA, err := Initalize()
	if err != nil {
		t.Errorf(fmt.Sprintf("Failed to initialize lookup table: %s", err))
	}
	var dG float64
	var expected_dG float64
	rRNA := "ACCTCCTTA"

	// 1.AAAAAAAAATG-AAATTGGGGGG
	dG = dG_rRNA_mRNA[rRNA]["AAATTAGGGTG"]
	expected_dG = -2.9
	if dG != expected_dG {
		t.Errorf(fmt.Sprintf("Failed to get correct lookup value. Expected: %v, got: %v", expected_dG, dG))
	}

	// 2.AAATTGGGGGG-AATCGGGGGTG
	dG = dG_rRNA_mRNA[rRNA]["AATATTGGACT"]
	expected_dG = 2.7
	if dG != expected_dG {
		t.Errorf(fmt.Sprintf("Failed to get correct lookup value. Expected: %v, got: %v", expected_dG, dG))
	}

	// 3.AATCGGGGGTG-AAGTATCCATG
	dG = dG_rRNA_mRNA[rRNA]["AACTGAATTGC"]
	expected_dG = 0.4
	if dG != expected_dG {
		t.Errorf(fmt.Sprintf("Failed to get correct lookup value. Expected: %v, got: %v", expected_dG, dG))
	}

	// 4.AAGTATCCTTG-ATAGAGTATGG
	dG = dG_rRNA_mRNA[rRNA]["AAGCGTTTGGA"]
	expected_dG = 2.0
	if dG != expected_dG {
		t.Errorf(fmt.Sprintf("Failed to get correct lookup value. Expected: %v, got: %v", expected_dG, dG))
	}

	// 5.ATAGAGTAGTG-ATTCTTAGGTG
	dG = dG_rRNA_mRNA[rRNA]["ATTAGCTGTGC"]
	expected_dG = 1.3
	if dG != expected_dG {
		t.Errorf(fmt.Sprintf("Failed to get correct lookup value. Expected: %v, got: %v", expected_dG, dG))
	}

	// 6.ATTCTTAGGTG-ATCGGGGGGTG
	dG = dG_rRNA_mRNA[rRNA]["ATCCGGGAATG"]
	expected_dG = -4.3
	if dG != expected_dG {
		t.Errorf(fmt.Sprintf("Failed to get correct lookup value. Expected: %v, got: %v", expected_dG, dG))
	}

	// 7.ATCGGGGGGTG-AGGGGGGGGGG
	dG = dG_rRNA_mRNA[rRNA]["AGGAGTTGCAC"]
	expected_dG = -6.4
	if dG != expected_dG {
		t.Errorf(fmt.Sprintf("Failed to get correct lookup value. Expected: %v, got: %v", expected_dG, dG))
	}

	// 8.AGGGGGGGGGG-TTTTTTTTTTT
	dG = dG_rRNA_mRNA[rRNA]["TAGTCGGGTGC"]
	expected_dG = 1.3
	if dG != expected_dG {
		t.Errorf(fmt.Sprintf("Failed to get correct lookup value. Expected: %v, got: %v", expected_dG, dG))
	}

	// 9.TTTTTTTTTTT-TCCCCCCCCCC
	dG = dG_rRNA_mRNA[rRNA]["TCCCAGGTGAA"]
	expected_dG = 0.4
	if dG != expected_dG {
		t.Errorf(fmt.Sprintf("Failed to get correct lookup value. Expected: %v, got: %v", expected_dG, dG))
	}

	// 10.TCCCCCCCCCC-TGGGGGGGGGG
	dG = dG_rRNA_mRNA[rRNA]["TCCCCGAGGTG"]
	expected_dG = -3.0
	if dG != expected_dG {
		t.Errorf(fmt.Sprintf("Failed to get correct lookup value. Expected: %v, got: %v", expected_dG, dG))
	}

	// 11.TGGGGGGGGGG-CATAATTGTAT
	dG = dG_rRNA_mRNA[rRNA]["CAATAGTATTG"]
	expected_dG = 0.7
	if dG != expected_dG {
		t.Errorf(fmt.Sprintf("Failed to get correct lookup value. Expected: %v, got: %v", expected_dG, dG))
	}

}
