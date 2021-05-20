package poly

import (
	"bytes"
	"log"
	"os"
	"strings"
	"testing"
)

// func ExampleCalculateMfe() {
// 	mfe, _ := CalculateMfe("ACGAUCAGAGAUCAGAGCAUACGACAGCAG", "..((((...))))...((........))..")
// 	fmt.Println(mfe)
// 	// Output:
// 	// -2.9
// }

func TestCalculateMfe(t *testing.T) {
	test("ACGAUCAGAGAUCAGAGCAUACGACAGCAG",
		"..((((...))))...((........))..",
		`External loop                           :  -300
		Interior loop (  3, 13) GC; (  4, 12) AU:  -240
		Interior loop (  4, 12) AU; (  5, 11) UA:  -110
		Interior loop (  5, 11) UA; (  6, 10) CG:  -240
		Hairpin  loop (  6, 10) CG              :   540
		Interior loop ( 17, 28) GC; ( 18, 27) CG:  -340
		Hairpin  loop ( 18, 27) CG              :   400
		-2.9`, t)

	test("AAAACGGUCCUUAUCAGGACCAAACA",
		".....((((((....)))))).....",
		`External loop                           :  -150
		Interior loop (  6, 21) GC; (  7, 20) GC:  -330
		Interior loop (  7, 20) GC; (  8, 19) UA:  -220
		Interior loop (  8, 19) UA; (  9, 18) CG:  -240
		Interior loop (  9, 18) CG; ( 10, 17) CG:  -330
		Interior loop ( 10, 17) CG; ( 11, 16) UA:  -210
		Hairpin  loop ( 11, 16) UA              :   550
		-9.3
		`, t)

	test("AUUCUUGCUUCAACAGUGUUUGAACGGAAU",
		"..............................",
		`External loop                           :     0
		0`, t)

	test("UCGGCCACAAACACACAAUCUACUGUUGGUCGA",
		"(((((((...................)))))))",
		`External loop                           :    50
		Interior loop (  1, 33) UA; (  2, 32) CG:  -240
		Interior loop (  2, 32) CG; (  3, 31) GC:  -240
		Interior loop (  3, 31) GC; (  4, 30) GU:  -150
		Interior loop (  4, 30) GU; (  5, 29) CG:  -250
		Interior loop (  5, 29) CG; (  6, 28) CG:  -330
		Interior loop (  6, 28) CG; (  7, 27) AU:  -210
		Hairpin  loop (  7, 27) AU              :   700
		-6.7`, t)

	test("GUUUUUAUCUUACACACGCUUGUGUAAGAUAGUUA",
		".....(((((((((((....)))))))))))....",
		`External loop                           :   -50
		Interior loop (  6, 31) UA; (  7, 30) AU:  -130
		Interior loop (  7, 30) AU; (  8, 29) UA:  -110
		Interior loop (  8, 29) UA; (  9, 28) CG:  -240
		Interior loop (  9, 28) CG; ( 10, 27) UA:  -210
		Interior loop ( 10, 27) UA; ( 11, 26) UA:   -90
		Interior loop ( 11, 26) UA; ( 12, 25) AU:  -130
		Interior loop ( 12, 25) AU; ( 13, 24) CG:  -220
		Interior loop ( 13, 24) CG; ( 14, 23) AU:  -210
		Interior loop ( 14, 23) AU; ( 15, 22) CG:  -220
		Interior loop ( 15, 22) CG; ( 16, 21) AU:  -210
		Hairpin  loop ( 16, 21) AU              :   540
		-12.8`, t)

	test("GGGCUCGUAGAUCAGCGGUAGAUCGCUUCCUUCGCAAGGAAGCCCUGGGUUCAAAUCCCAGCGAGUCCACCA",
		"(((((((..((((.......))))(((((((.....))))))).(((((.......))))))))))))....",
		`External loop                           :  -170
		Interior loop (  1, 68) GC; (  2, 67) GC:  -330
		Interior loop (  2, 67) GC; (  3, 66) GU:  -150
		Interior loop (  3, 66) GU; (  4, 65) CG:  -250
		Interior loop (  4, 65) CG; (  5, 64) UA:  -210
		Interior loop (  5, 64) UA; (  6, 63) CG:  -240
		Interior loop (  6, 63) CG; (  7, 62) GC:  -240
		Interior loop ( 10, 24) GC; ( 11, 23) AU:  -240
		Interior loop ( 11, 23) AU; ( 12, 22) UA:  -110
		Interior loop ( 12, 22) UA; ( 13, 21) CG:  -240
		Hairpin  loop ( 13, 21) CG              :   450
		Interior loop ( 25, 43) GC; ( 26, 42) CG:  -340
		Interior loop ( 26, 42) CG; ( 27, 41) UA:  -210
		Interior loop ( 27, 41) UA; ( 28, 40) UA:   -90
		Interior loop ( 28, 40) UA; ( 29, 39) CG:  -240
		Interior loop ( 29, 39) CG; ( 30, 38) CG:  -330
		Interior loop ( 30, 38) CG; ( 31, 37) UA:  -210
		Hairpin  loop ( 31, 37) UA              :   550
		Interior loop ( 45, 61) CG; ( 46, 60) UA:  -210
		Interior loop ( 46, 60) UA; ( 47, 59) GC:  -210
		Interior loop ( 47, 59) GC; ( 48, 58) GC:  -330
		Interior loop ( 48, 58) GC; ( 49, 57) GC:  -330
		Hairpin  loop ( 49, 57) GC              :   440
		Multi    loop (  7, 62) GC              :   140
		-31`, t)

	// test("AUGAAACAAUACCAAGAUUUAAUUAAAGACAUUUUUGAAAAUGGUUAUGAAACCGAUGAUCGUACAGGCACAGGAACAAUUGCUCUGUUCGGAUCUAAAUUACGCUGGGAUUUAACUAAAGGUUUUCCUGCGGUAACAACUAAGAAGCUCGCCUGGAAAGCUUGCAUUGCUGAGCUAAUAUGGUUUUUAUCAGGAAGCACAAAUGUCAAUGAUUUACGAUUAAUUCAACACGAUUCGUUAAUCCAAGGCAAAACAGUCUGGGAUGAAAAUUACGAAAAUCAAGCAAAAGAUUUAGGAUACCAUAGCGGUGAACUUGGUCCAAUUUAUGGAAAACAGUGGCGUGAUUUUGGUGGUGUAGACCAAAUUAUAGAAGUUAUUGAUCGUAUUAAAAAACUGCCAAAUGAUAGGCGUCAAAUUGUUUCUGCAUGGAAUCCAGCUGAACUUAAAUAUAUGGCAUUACCGCCUUGUCAUAUGUUCUAUCAGUUUAAUGUGCGUAAUGGCUAUUUGGAUUUGCAGUGGUAUCAACGCUCAGUAGAUGUUUUCUUGGGUCUACCGUUUAAUAUUGCGUCAUAUGCUACGUUAGUUCAUAUUGUAGCUAAGAUGUGUAAUCUUAUUCCAGGGGAUUUGAUAUUUUCUGGUGGUAAUACUCAUAUCUAUAUGAAUCACGUAGAACAAUGUAAAGAAAUUUUGAGGCGUGAACCUAAAGAGCUUUGUGAGCUGGUAAUAAGUGGUCUACCUUAUAAAUUCCGAUAUCUUUCUACUAAAGAACAAUUAAAAUAUGUUCUUAAACUUAGGCCUAAAGAUUUCGUUCUUAACAACUAUGUAUCACACCCUCCUAUUAAAGGAAAGAUGGCGGUGUAA",
	// 	"........(((((...((((.......(((((((.((....((((...((((((...((((..(((((((..........))).))))...))))...((((.(((.((..........(((((((((.(((((............)).)))..)))))))))......)).)))))))..)))))).)))).....)).)))))))..(((.(((....))).))).......((((..((((.((((......)))).)))).......))))(((((........)))))..(((((.((((..((((((....((((.....)))).....((((((((((.((((....(((((((............(((((((.(((........(((((.(((((((((.(((...........((((((......((((((.....(((((((((((.........)))))))))))...))))))...))))))..))).)))))...))))))))).......))).)))))))...........)))))))...)))).))))))))))...((((((..((((((((........(((((((.....((((((((((((((.....))))))))).)))))....))))))))))))))).))))))...........((((((((.((((.(((........(((((...)))))((((..(((..(((...............)))..)))..)))).(((((((.........)))))))....))).)))).))))))))))))...))..)))))))))......((((.....)))).))))...)))))..",
	// 	-197.7, t)
}

func test(sequence, structure string, expected_output string, t *testing.T) {
	output := captureOutput(func() {
		mfe, _ := CalculateMfe(sequence, structure)
		log.Printf("%v", mfe)
	})
	output = stripWhiteSpace(output)
	expected_output = stripWhiteSpace(expected_output)
	if output != expected_output {
		t.Errorf("\n\nFailed to calcualte mfe for '%v'. \nExpected: \n'%v'\n\nGot: \n'%v'\n\n",
			sequence, expected_output, output)
	}
}

func captureOutput(f func()) string {
	var buf bytes.Buffer
	log.SetFlags(0)
	log.SetOutput(&buf)
	f()
	log.SetOutput(os.Stderr)
	return buf.String()
}

func stripWhiteSpace(s string) string {
	return strings.Join(strings.Fields(s), "")
}
