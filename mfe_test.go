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

	test("AUGAAACAAUACCAAGAUUUAAUUAAAGACAUUUUUGAAAAUGGUUAUGAAACCGAUGAUCGUACAGGCACAGGAACAAUUGCUCUGUUCGGAUCUAAAUUACGCUGGGAUUUAACUAAAGGUUUUCCUGCGGUAACAACUAAGAAGCUCGCCUGGAAAGCUUGCAUUGCUGAGCUAAUAUGGUUUUUAUCAGGAAGCACAAAUGUCAAUGAUUUACGAUUAAUUCAACACGAUUCGUUAAUCCAAGGCAAAACAGUCUGGGAUGAAAAUUACGAAAAUCAAGCAAAAGAUUUAGGAUACCAUAGCGGUGAACUUGGUCCAAUUUAUGGAAAACAGUGGCGUGAUUUUGGUGGUGUAGACCAAAUUAUAGAAGUUAUUGAUCGUAUUAAAAAACUGCCAAAUGAUAGGCGUCAAAUUGUUUCUGCAUGGAAUCCAGCUGAACUUAAAUAUAUGGCAUUACCGCCUUGUCAUAUGUUCUAUCAGUUUAAUGUGCGUAAUGGCUAUUUGGAUUUGCAGUGGUAUCAACGCUCAGUAGAUGUUUUCUUGGGUCUACCGUUUAAUAUUGCGUCAUAUGCUACGUUAGUUCAUAUUGUAGCUAAGAUGUGUAAUCUUAUUCCAGGGGAUUUGAUAUUUUCUGGUGGUAAUACUCAUAUCUAUAUGAAUCACGUAGAACAAUGUAAAGAAAUUUUGAGGCGUGAACCUAAAGAGCUUUGUGAGCUGGUAAUAAGUGGUCUACCUUAUAAAUUCCGAUAUCUUUCUACUAAAGAACAAUUAAAAUAUGUUCUUAAACUUAGGCCUAAAGAUUUCGUUCUUAACAACUAUGUAUCACACCCUCCUAUUAAAGGAAAGAUGGCGGUGUAA",
		"........(((((...((((.......(((((((.((....((((...((((((...((((..(((((((..........))).))))...))))...((((.(((.((..........(((((((((.(((((............)).)))..)))))))))......)).)))))))..)))))).)))).....)).)))))))..(((.(((....))).))).......((((..((((.((((......)))).)))).......))))(((((........)))))..(((((.((((..((((((....((((.....)))).....((((((((((.((((....(((((((............(((((((.(((........(((((.(((((((((.(((...........((((((......((((((.....(((((((((((.........)))))))))))...))))))...))))))..))).)))))...))))))))).......))).)))))))...........)))))))...)))).))))))))))...((((((..((((((((........(((((((.....((((((((((((((.....))))))))).)))))....))))))))))))))).))))))...........((((((((.((((.(((........(((((...)))))((((..(((..(((...............)))..)))..)))).(((((((.........)))))))....))).)))).))))))))))))...))..)))))))))......((((.....)))).))))...)))))..",
		`External loop                           :   -50
		Interior loop (  9,859) AU; ( 10,858) UG:  -140
		Interior loop ( 10,858) UG; ( 11,857) AU:  -100
		Interior loop ( 11,857) AU; ( 12,856) CG:  -220
		Interior loop ( 12,856) CG; ( 13,855) CG:  -330
		Interior loop ( 13,855) CG; ( 17,851) AU:   170
		Interior loop ( 17,851) AU; ( 18,850) UA:  -110
		Interior loop ( 18,850) UA; ( 19,849) UG:  -130
		Interior loop ( 19,849) UG; ( 20,848) UA:   -60
		Interior loop ( 28,207) GC; ( 29,206) AU:  -240
		Interior loop ( 29,206) AU; ( 30,205) CG:  -220
		Interior loop ( 30,205) CG; ( 31,204) AU:  -210
		Interior loop ( 31,204) AU; ( 32,203) UA:  -110
		Interior loop ( 32,203) UA; ( 33,202) UA:   -90
		Interior loop ( 33,202) UA; ( 34,201) UA:   -90
		Interior loop ( 34,201) UA; ( 36,199) UA:   190
		Interior loop ( 36,199) UA; ( 37,198) GC:  -210
		Interior loop ( 37,198) GC; ( 42,192) UA:   190
		Interior loop ( 42,192) UA; ( 43,191) GC:  -210
		Interior loop ( 43,191) GC; ( 44,190) GU:  -150
		Interior loop ( 44,190) GU; ( 45,189) UA:  -140
		Interior loop ( 45,189) UA; ( 49,187) GU:   370
		Interior loop ( 49,187) GU; ( 50,186) AU:  -130
		Interior loop ( 50,186) AU; ( 51,185) AU:   -90
		Interior loop ( 51,185) AU; ( 52,184) AU:   -90
		Interior loop ( 52,184) AU; ( 53,183) CG:  -220
		Interior loop ( 53,183) CG; ( 54,182) CG:  -330
		Interior loop ( 58, 95) GC; ( 59, 94) AU:  -240
		Interior loop ( 59, 94) AU; ( 60, 93) UA:  -110
		Interior loop ( 60, 93) UA; ( 61, 92) CG:  -240
		Interior loop ( 61, 92) CG; ( 64, 88) AU:   230
		Interior loop ( 64, 88) AU; ( 65, 87) CG:  -220
		Interior loop ( 65, 87) CG; ( 66, 86) AU:  -210
		Interior loop ( 66, 86) AU; ( 67, 85) GC:  -210
		Interior loop ( 67, 85) GC; ( 68, 83) GC:    50
		Interior loop ( 68, 83) GC; ( 69, 82) CG:  -340
		Interior loop ( 69, 82) CG; ( 70, 81) AU:  -210
		Hairpin  loop ( 70, 81) AU              :   630
		Interior loop ( 99,179) AU; (100,178) UA:  -110
		Interior loop (100,178) UA; (101,177) UA:   -90
		Interior loop (101,177) UA; (102,176) AU:  -130
		Interior loop (102,176) AU; (104,175) GC:   170
		Interior loop (104,175) GC; (105,174) CG:  -340
		Interior loop (105,174) CG; (106,173) UA:  -210
		Interior loop (106,173) UA; (108,171) GU:   -70
		Interior loop (108,171) GU; (109,170) GC:  -210
		Interior loop (109,170) GC; (120,163) AU:   430
		Interior loop (120,163) AU; (121,162) GU:   -60
		Interior loop (121,162) GU; (122,161) GC:  -210
		Interior loop (122,161) GC; (123,160) UG:  -250
		Interior loop (123,160) UG; (124,159) UA:   -60
		Interior loop (124,159) UA; (125,158) UA:   -90
		Interior loop (125,158) UA; (126,157) UA:   -90
		Interior loop (126,157) UA; (127,156) CG:  -240
		Interior loop (127,156) CG; (128,155) CG:  -330
		Interior loop (128,155) CG; (130,152) GC:   150
		Interior loop (130,152) GC; (131,151) CG:  -340
		Interior loop (131,151) CG; (132,150) GC:  -240
		Interior loop (132,150) GC; (133,148) GC:    50
		Interior loop (133,148) GC; (134,147) UG:  -250
		Hairpin  loop (134,147) UG              :   620
		Multi    loop ( 54,182) CG              :   340
		Interior loop (210,227) UA; (211,226) GC:  -210
		Interior loop (211,226) GC; (212,225) AU:  -240
		Interior loop (212,225) AU; (214,223) UA:   120
		Interior loop (214,223) UA; (215,222) UA:   -90
		Interior loop (215,222) UA; (216,221) AU:  -130
		Hairpin  loop (216,221) AU              :   540
		Interior loop (235,275) UA; (236,274) CG:  -240
		Interior loop (236,274) CG; (237,273) GC:  -240
		Interior loop (237,273) GC; (238,272) UA:  -220
		Interior loop (238,272) UA; (241,264) AU:   520
		Interior loop (241,264) AU; (242,263) UA:  -110
		Interior loop (242,263) UA; (243,262) CG:  -240
		Interior loop (243,262) CG; (244,261) CG:  -330
		Interior loop (244,261) CG; (246,259) AU:   120
		Interior loop (246,259) AU; (247,258) GC:  -210
		Interior loop (247,258) GC; (248,257) GU:  -150
		Interior loop (248,257) GU; (249,256) CG:  -250
		Hairpin  loop (249,256) CG              :   390
		Interior loop (276,293) AU; (277,292) AU:   -90
		Interior loop (277,292) AU; (278,291) AU:   -90
		Interior loop (278,291) AU; (279,290) UA:  -110
		Interior loop (279,290) UA; (280,289) CG:  -240
		Hairpin  loop (280,289) CG              :   400
		Interior loop (296,827) GC; (297,826) AU:  -240
		Interior loop (297,826) AU; (298,825) UA:  -110
		Interior loop (298,825) UA; (299,824) AU:  -130
		Interior loop (299,824) AU; (300,823) CG:  -220
		Interior loop (300,823) CG; (302,822) AU:   170
		Interior loop (302,822) AU; (303,821) UA:  -110
		Interior loop (303,821) UA; (304,820) AU:  -130
		Interior loop (304,820) AU; (305,819) GC:  -210
		Interior loop (305,819) GC; (308,816) GC:   -70
		Interior loop (308,816) GC; (309,815) UA:  -220
		Interior loop (309,815) UA; (310,811) GC:   370
		Interior loop (310,811) GC; (311,810) AU:  -240
		Interior loop (311,810) AU; (312,809) AU:   -90
		Interior loop (312,809) AU; (313,808) CG:  -220
		Interior loop (318,330) UA; (319,329) CG:  -240
		Interior loop (319,329) CG; (320,328) CG:  -330
		Interior loop (320,328) CG; (321,327) AU:  -210
		Hairpin  loop (321,327) AU              :   540
		Interior loop (336,571) GU; (337,570) UA:  -140
		Interior loop (337,570) UA; (338,569) GC:  -210
		Interior loop (338,569) GC; (339,568) GU:  -150
		Interior loop (339,568) GU; (340,567) CG:  -250
		Interior loop (340,567) CG; (341,566) GC:  -240
		Interior loop (341,566) GC; (342,565) UG:  -250
		Interior loop (342,565) UG; (343,564) GU:    30
		Interior loop (343,564) GU; (344,563) AU:  -130
		Interior loop (344,563) AU; (345,562) UA:  -110
		Interior loop (345,562) UA; (347,560) UA:   150
		Interior loop (347,560) UA; (348,559) UA:   -90
		Interior loop (348,559) UA; (349,558) GU:  -100
		Interior loop (349,558) GU; (350,557) GU:   -50
		Interior loop (350,557) GU; (355,553) GC:   280
		Interior loop (355,553) GC; (356,552) UA:  -220
		Interior loop (356,552) UA; (357,551) AU:  -130
		Interior loop (357,551) AU; (358,550) GC:  -210
		Interior loop (358,550) GC; (359,549) AU:  -240
		Interior loop (359,549) AU; (360,548) CG:  -220
		Interior loop (360,548) CG; (361,547) CG:  -330
		Interior loop (361,547) CG; (374,535) UG:   310
		Interior loop (374,535) UG; (375,534) UA:   -60
		Interior loop (375,534) UA; (376,533) AU:  -130
		Interior loop (376,533) AU; (377,532) UG:  -140
		Interior loop (377,532) UG; (378,531) UA:   -60
		Interior loop (378,531) UA; (379,530) GC:  -210
		Interior loop (379,530) GC; (380,529) AU:  -240
		Interior loop (380,529) AU; (382,527) CG:   120
		Interior loop (382,527) CG; (383,526) GC:  -240
		Interior loop (383,526) GC; (384,525) UA:  -220
		Interior loop (384,525) UA; (393,517) AU:   390
		Interior loop (393,517) AU; (394,516) CG:  -220
		Interior loop (394,516) CG; (395,515) UA:  -210
		Interior loop (395,515) UA; (396,514) GC:  -210
		Interior loop (396,514) GC; (397,513) CG:  -340
		Interior loop (397,513) CG; (399,512) AU:   170
		Interior loop (399,512) AU; (400,511) AU:   -90
		Interior loop (400,511) AU; (401,510) AU:   -90
		Interior loop (401,510) AU; (402,509) UA:  -110
		Interior loop (402,509) UA; (403,505) GU:   420
		Interior loop (403,505) GU; (404,504) AU:  -130
		Interior loop (404,504) AU; (405,503) UA:  -110
		Interior loop (405,503) UA; (406,502) AU:  -130
		Interior loop (406,502) AU; (407,501) GC:  -210
		Interior loop (407,501) GC; (409,499) CG:  -230
		Interior loop (409,499) CG; (410,498) GU:  -140
		Interior loop (410,498) GU; (411,497) UA:  -140
		Interior loop (411,497) UA; (423,494) UG:   720
		Interior loop (423,494) UG; (424,493) GC:  -140
		Interior loop (424,493) GC; (425,492) CG:  -340
		Interior loop (425,492) CG; (426,491) AU:  -210
		Interior loop (426,491) AU; (427,490) UG:  -140
		Interior loop (427,490) UG; (428,489) GU:    30
		Interior loop (428,489) GU; (435,485) AU:   460
		Interior loop (435,485) AU; (436,484) GU:   -60
		Interior loop (436,484) GU; (437,483) CG:  -250
		Interior loop (437,483) CG; (438,482) UA:  -210
		Interior loop (438,482) UA; (439,481) GC:  -210
		Interior loop (439,481) GC; (440,480) AU:  -240
		Interior loop (440,480) AU; (446,476) AU:   490
		Interior loop (446,476) AU; (447,475) AU:   -90
		Interior loop (447,475) AU; (448,474) UG:  -140
		Interior loop (448,474) UG; (449,473) AU:  -100
		Interior loop (449,473) AU; (450,472) UA:  -110
		Interior loop (450,472) UA; (451,471) AU:  -130
		Interior loop (451,471) AU; (452,470) UA:  -110
		Interior loop (452,470) UA; (453,469) GC:  -210
		Interior loop (453,469) GC; (454,468) GU:  -150
		Interior loop (454,468) GU; (455,467) CG:  -250
		Interior loop (455,467) CG; (456,466) AU:  -210
		Hairpin  loop (456,466) AU              :   520
		Interior loop (575,670) CG; (576,669) UA:  -210
		Interior loop (576,669) UA; (577,668) AU:  -130
		Interior loop (577,668) AU; (578,667) CG:  -220
		Interior loop (578,667) CG; (579,666) GC:  -240
		Interior loop (579,666) GC; (580,665) UA:  -220
		Interior loop (580,665) UA; (583,663) GU:   370
		Interior loop (583,663) GU; (584,662) UA:  -140
		Interior loop (584,662) UA; (585,661) UA:   -90
		Interior loop (585,661) UA; (586,660) CG:  -240
		Interior loop (586,660) CG; (587,659) AU:  -210
		Interior loop (587,659) AU; (588,658) UA:  -110
		Interior loop (588,658) UA; (589,657) AU:  -130
		Interior loop (589,657) AU; (590,656) UA:  -110
		Interior loop (590,656) UA; (599,655) AU:   570
		Interior loop (599,655) AU; (600,654) GC:  -210
		Interior loop (600,654) GC; (601,653) AU:  -240
		Interior loop (601,653) AU; (602,652) UA:  -110
		Interior loop (602,652) UA; (603,651) GU:  -100
		Interior loop (603,651) GU; (604,650) UA:  -140
		Interior loop (604,650) UA; (605,649) GC:  -210
		Interior loop (605,649) GC; (611,644) UA:   310
		Interior loop (611,644) UA; (612,643) UA:   -90
		Interior loop (612,643) UA; (613,642) AU:  -130
		Interior loop (613,642) AU; (614,641) UG:  -140
		Interior loop (614,641) UG; (615,640) UG:   -50
		Interior loop (615,640) UG; (616,638) CG:   230
		Interior loop (616,638) CG; (617,637) CG:  -330
		Interior loop (617,637) CG; (618,636) AU:  -210
		Interior loop (618,636) AU; (619,635) GC:  -210
		Interior loop (619,635) GC; (620,634) GU:  -150
		Interior loop (620,634) GU; (621,633) GU:   -50
		Interior loop (621,633) GU; (622,632) GU:   -50
		Interior loop (622,632) GU; (623,631) AU:  -130
		Interior loop (623,631) AU; (624,630) UA:  -110
		Hairpin  loop (624,630) UA              :   480
		Interior loop (682,807) GC; (683,806) AU:  -240
		Interior loop (683,806) AU; (684,805) AU:   -90
		Interior loop (684,805) AU; (685,804) AU:   -90
		Interior loop (685,804) AU; (686,803) UA:  -110
		Interior loop (686,803) UA; (687,802) UG:  -130
		Interior loop (687,802) UG; (688,801) UA:   -60
		Interior loop (688,801) UA; (689,800) UA:   -90
		Interior loop (689,800) UA; (691,798) AU:   190
		Interior loop (691,798) AU; (692,797) GC:  -210
		Interior loop (692,797) GC; (693,796) GC:  -330
		Interior loop (693,796) GC; (694,795) CG:  -340
		Interior loop (694,795) CG; (696,793) UA:  -140
		Interior loop (696,793) UA; (697,792) GU:  -100
		Interior loop (697,792) GU; (698,791) AU:  -130
		Interior loop (707,719) AU; (708,718) GC:  -210
		Interior loop (708,718) GC; (709,717) CG:  -340
		Interior loop (709,717) CG; (710,716) UA:  -210
		Interior loop (710,716) UA; (711,715) UG:  -130
		Hairpin  loop (711,715) UG              :   590
		Interior loop (720,762) GU; (721,761) GC:  -210
		Interior loop (721,761) GC; (722,760) UA:  -220
		Interior loop (722,760) UA; (723,759) AU:  -130
		Interior loop (723,759) AU; (726,756) AU:   180
		Interior loop (726,756) AU; (727,755) AU:   -90
		Interior loop (727,755) AU; (728,754) GC:  -210
		Interior loop (728,754) GC; (731,751) GU:   140
		Interior loop (731,751) GU; (732,750) UA:  -140
		Interior loop (732,750) UA; (733,749) CG:  -240
		Hairpin  loop (733,749) CG              :   550
		Interior loop (764,786) AU; (765,785) AU:   -90
		Interior loop (765,785) AU; (766,784) GC:  -210
		Interior loop (766,784) GC; (767,783) AU:  -240
		Interior loop (767,783) AU; (768,782) AU:   -90
		Interior loop (768,782) AU; (769,781) CG:  -220
		Interior loop (769,781) CG; (770,780) AU:  -210
		Hairpin  loop (770,780) AU              :   610
		Multi    loop (698,791) AU              :   370
		Multi    loop (313,808) CG              :   -10
		Interior loop (834,846) UA; (835,845) CG:  -240
		Interior loop (835,845) CG; (836,844) CG:  -330
		Interior loop (836,844) CG; (837,843) UA:  -210
		Hairpin  loop (837,843) UA              :   520
		Multi    loop ( 20,848) UA              :  -220
		-197.7`, t)

	// test("")
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
