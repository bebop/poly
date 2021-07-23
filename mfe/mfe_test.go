package mfe

import (
	"bytes"
	"fmt"
	"log"
	"os"
	"strings"
	"testing"

	"github.com/TimothyStiles/poly/energy_params"
	. "github.com/TimothyStiles/poly/secondary_structure"
)

func ExampleMinimumFreeEnergy() {
	mfe, _, _ := MinimumFreeEnergy("ACGAUCAGAGAUCAGAGCAUACGACAGCAG", "..((((...))))...((........))..", DefaultTemperature, energy_params.Turner2004, DefaultDanglingEndsModel)
	fmt.Println(mfe)
	// Output:
	// -2.9
}

func TestAndronescu2007ParamsNoDanglingEnergies(t *testing.T) {
	compareMFEOutputToViennaRNA("ACGAUCAGAGAUCAGAGCAUACGACAGCAG",
		"..((((...))))...((........))..",
		DefaultTemperature, energy_params.Andronescu2007, NoDanglingEnds,
		`External loop                           :     0
		Interior loop (  3, 13) GC; (  4, 12) AU:  -211
		Interior loop (  4, 12) AU; (  5, 11) UA:   -99
		Interior loop (  5, 11) UA; (  6, 10) CG:  -211
		Hairpin  loop (  6, 10) CG              :   475
		Interior loop ( 17, 28) GC; ( 18, 27) CG:  -300
		Hairpin  loop ( 18, 27) CG              :   324
		-0.22`, t)

	compareMFEOutputToViennaRNA(
		"AUGAAACAAUACCAAGAUUUAAUUAAAGACAUUUUUGAAAAUGGUUAUGAAACCGAUGAUCGUACAGGCACAGGAACAAUUGCUCUGUUCGGAUCUAAAUUACGCUGGGAUUUAACUAAAGGUUUUCCUGCGGUAACAACUAAGAAGCUCGCCUGGAAAGCUUGCAUUGCUGAGCUAAUAUGGUUUUUAUCAGGAAGCACAAAUGUCAAUGAUUUACGAUUAAUUCAACACGAUUCGUUAAUCCAAGGCAAAACAGUCUGGGAUGAAAAUUACGAAAAUCAAGCAAAAGAUUUAGGAUACCAUAGCGGUGAACUUGGUCCAAUUUAUGGAAAACAGUGGCGUGAUUUUGGUGGUGUAGACCAAAUUAUAGAAGUUAUUGAUCGUAUUAAAAAACUGCCAAAUGAUAGGCGUCAAAUUGUUUCUGCAUGGAAUCCAGCUGAACUUAAAUAUAUGGCAUUACCGCCUUGUCAUAUGUUCUAUCAGUUUAAUGUGCGUAAUGGCUAUUUGGAUUUGCAGUGGUAUCAACGCUCAGUAGAUGUUUUCUUGGGUCUACCGUUUAAUAUUGCGUCAUAUGCUACGUUAGUUCAUAUUGUAGCUAAGAUGUGUAAUCUUAUUCCAGGGGAUUUGAUAUUUUCUGGUGGUAAUACUCAUAUCUAUAUGAAUCACGUAGAACAAUGUAAAGAAAUUUUGAGGCGUGAACCUAAAGAGCUUUGUGAGCUGGUAAUAAGUGGUCUACCUUAUAAAUUCCGAUAUCUUUCUACUAAAGAACAAUUAAAAUAUGUUCUUAAACUUAGGCCUAAAGAUUUCGUUCUUAACAACUAUGUAUCACACCCUCCUAUUAAAGGAAAGAUGGCGGUGUAA",
		"........(((((...((((.......(((((((.((....((((...((((((...((((..(((((((..........))).))))...))))...((((.(((.((..........(((((((((.(((((............)).)))..)))))))))......)).)))))))..)))))).)))).....)).)))))))..(((.(((....))).))).......((((..((((.((((......)))).)))).......))))(((((........)))))..(((((.((((..((((((....((((.....)))).....((((((((((.((((....(((((((............(((((((.(((........(((((.(((((((((.(((...........((((((......((((((.....(((((((((((.........)))))))))))...))))))...))))))..))).)))))...))))))))).......))).)))))))...........)))))))...)))).))))))))))...((((((..((((((((........(((((((.....((((((((((((((.....))))))))).)))))....))))))))))))))).))))))...........((((((((.((((.(((........(((((...)))))((((..(((..(((...............)))..)))..)))).(((((((.........)))))))....))).)))).))))))))))))...))..)))))))))......((((.....)))).))))...)))))..",
		DefaultTemperature, energy_params.Andronescu2007, NoDanglingEnds,
		`External loop                           :    11
		Interior loop (  9,859) AU; ( 10,858) UG:   -88
		Interior loop ( 10,858) UG; ( 11,857) AU:     0
		Interior loop ( 11,857) AU; ( 12,856) CG:  -189
		Interior loop ( 12,856) CG; ( 13,855) CG:  -271
		Interior loop ( 13,855) CG; ( 17,851) AU:   166
		Interior loop ( 17,851) AU; ( 18,850) UA:   -99
		Interior loop ( 18,850) UA; ( 19,849) UG:   -47
		Interior loop ( 19,849) UG; ( 20,848) UA:    -1
		Interior loop ( 28,207) GC; ( 29,206) AU:  -211
		Interior loop ( 29,206) AU; ( 30,205) CG:  -189
		Interior loop ( 30,205) CG; ( 31,204) AU:  -178
		Interior loop ( 31,204) AU; ( 32,203) UA:   -99
		Interior loop ( 32,203) UA; ( 33,202) UA:   -69
		Interior loop ( 33,202) UA; ( 34,201) UA:   -69
		Interior loop ( 34,201) UA; ( 36,199) UA:   128
		Interior loop ( 36,199) UA; ( 37,198) GC:  -178
		Interior loop ( 37,198) GC; ( 42,192) UA:   108
		Interior loop ( 42,192) UA; ( 43,191) GC:  -178
		Interior loop ( 43,191) GC; ( 44,190) GU:  -127
		Interior loop ( 44,190) GU; ( 45,189) UA:   -88
		Interior loop ( 45,189) UA; ( 49,187) GU:   166
		Interior loop ( 49,187) GU; ( 50,186) AU:   -47
		Interior loop ( 50,186) AU; ( 51,185) AU:   -69
		Interior loop ( 51,185) AU; ( 52,184) AU:   -69
		Interior loop ( 52,184) AU; ( 53,183) CG:  -189
		Interior loop ( 53,183) CG; ( 54,182) CG:  -271
		Interior loop ( 58, 95) GC; ( 59, 94) AU:  -211
		Interior loop ( 59, 94) AU; ( 60, 93) UA:   -99
		Interior loop ( 60, 93) UA; ( 61, 92) CG:  -211
		Interior loop ( 61, 92) CG; ( 64, 88) AU:   167
		Interior loop ( 64, 88) AU; ( 65, 87) CG:  -189
		Interior loop ( 65, 87) CG; ( 66, 86) AU:  -178
		Interior loop ( 66, 86) AU; ( 67, 85) GC:  -199
		Interior loop ( 67, 85) GC; ( 68, 83) GC:    40
		Interior loop ( 68, 83) GC; ( 69, 82) CG:  -300
		Interior loop ( 69, 82) CG; ( 70, 81) AU:  -178
		Hairpin  loop ( 70, 81) AU              :   483
		Interior loop ( 99,179) AU; (100,178) UA:   -99
		Interior loop (100,178) UA; (101,177) UA:   -69
		Interior loop (101,177) UA; (102,176) AU:   -91
		Interior loop (102,176) AU; (104,175) GC:   112
		Interior loop (104,175) GC; (105,174) CG:  -300
		Interior loop (105,174) CG; (106,173) UA:  -199
		Interior loop (106,173) UA; (108,171) GU:    62
		Interior loop (108,171) GU; (109,170) GC:  -178
		Interior loop (109,170) GC; (120,163) AU:   320
		Interior loop (120,163) AU; (121,162) GU:    -1
		Interior loop (121,162) GU; (122,161) GC:  -178
		Interior loop (122,161) GC; (123,160) UG:  -193
		Interior loop (123,160) UG; (124,159) UA:    -1
		Interior loop (124,159) UA; (125,158) UA:   -69
		Interior loop (125,158) UA; (126,157) UA:   -69
		Interior loop (126,157) UA; (127,156) CG:  -211
		Interior loop (127,156) CG; (128,155) CG:  -271
		Interior loop (128,155) CG; (130,152) GC:   175
		Interior loop (130,152) GC; (131,151) CG:  -300
		Interior loop (131,151) CG; (132,150) GC:  -203
		Interior loop (132,150) GC; (133,148) GC:    40
		Interior loop (133,148) GC; (134,147) UG:  -193
		Hairpin  loop (134,147) UG              :   503
		Multi    loop ( 54,182) CG              :   492
		Interior loop (210,227) UA; (211,226) GC:  -178
		Interior loop (211,226) GC; (212,225) AU:  -211
		Interior loop (212,225) AU; (214,223) UA:   137
		Interior loop (214,223) UA; (215,222) UA:   -69
		Interior loop (215,222) UA; (216,221) AU:   -91
		Hairpin  loop (216,221) AU              :   464
		Interior loop (235,275) UA; (236,274) CG:  -211
		Interior loop (236,274) CG; (237,273) GC:  -203
		Interior loop (237,273) GC; (238,272) UA:  -189
		Interior loop (238,272) UA; (241,264) AU:   382
		Interior loop (241,264) AU; (242,263) UA:   -99
		Interior loop (242,263) UA; (243,262) CG:  -211
		Interior loop (243,262) CG; (244,261) CG:  -271
		Interior loop (244,261) CG; (246,259) AU:    62
		Interior loop (246,259) AU; (247,258) GC:  -199
		Interior loop (247,258) GC; (248,257) GU:  -127
		Interior loop (248,257) GU; (249,256) CG:  -193
		Hairpin  loop (249,256) CG              :   338
		Interior loop (276,293) AU; (277,292) AU:   -69
		Interior loop (277,292) AU; (278,291) AU:   -69
		Interior loop (278,291) AU; (279,290) UA:   -99
		Interior loop (279,290) UA; (280,289) CG:  -211
		Hairpin  loop (280,289) CG              :   324
		Interior loop (296,827) GC; (297,826) AU:  -211
		Interior loop (297,826) AU; (298,825) UA:   -99
		Interior loop (298,825) UA; (299,824) AU:   -91
		Interior loop (299,824) AU; (300,823) CG:  -189
		Interior loop (300,823) CG; (302,822) AU:   133
		Interior loop (302,822) AU; (303,821) UA:   -99
		Interior loop (303,821) UA; (304,820) AU:   -91
		Interior loop (304,820) AU; (305,819) GC:  -199
		Interior loop (305,819) GC; (308,816) GC:    22
		Interior loop (308,816) GC; (309,815) UA:  -189
		Interior loop (309,815) UA; (310,811) GC:   277
		Interior loop (310,811) GC; (311,810) AU:  -211
		Interior loop (311,810) AU; (312,809) AU:   -69
		Interior loop (312,809) AU; (313,808) CG:  -189
		Interior loop (318,330) UA; (319,329) CG:  -211
		Interior loop (319,329) CG; (320,328) CG:  -271
		Interior loop (320,328) CG; (321,327) AU:  -178
		Hairpin  loop (321,327) AU              :   447
		Interior loop (336,571) GU; (337,570) UA:   -88
		Interior loop (337,570) UA; (338,569) GC:  -178
		Interior loop (338,569) GC; (339,568) GU:  -127
		Interior loop (339,568) GU; (340,567) CG:  -193
		Interior loop (340,567) CG; (341,566) GC:  -203
		Interior loop (341,566) GC; (342,565) UG:  -193
		Interior loop (342,565) UG; (343,564) GU:   -70
		Interior loop (343,564) GU; (344,563) AU:   -47
		Interior loop (344,563) AU; (345,562) UA:   -99
		Interior loop (345,562) UA; (347,560) UA:   128
		Interior loop (347,560) UA; (348,559) UA:   -69
		Interior loop (348,559) UA; (349,558) GU:     0
		Interior loop (349,558) GU; (350,557) GU:   -71
		Interior loop (350,557) GU; (355,553) GC:   165
		Interior loop (355,553) GC; (356,552) UA:  -189
		Interior loop (356,552) UA; (357,551) AU:   -91
		Interior loop (357,551) AU; (358,550) GC:  -199
		Interior loop (358,550) GC; (359,549) AU:  -211
		Interior loop (359,549) AU; (360,548) CG:  -189
		Interior loop (360,548) CG; (361,547) CG:  -271
		Interior loop (361,547) CG; (374,535) UG:   209
		Interior loop (374,535) UG; (375,534) UA:    -1
		Interior loop (375,534) UA; (376,533) AU:   -91
		Interior loop (376,533) AU; (377,532) UG:   -88
		Interior loop (377,532) UG; (378,531) UA:    -1
		Interior loop (378,531) UA; (379,530) GC:  -178
		Interior loop (379,530) GC; (380,529) AU:  -211
		Interior loop (380,529) AU; (382,527) CG:    62
		Interior loop (382,527) CG; (383,526) GC:  -203
		Interior loop (383,526) GC; (384,525) UA:  -189
		Interior loop (384,525) UA; (393,517) AU:   305
		Interior loop (393,517) AU; (394,516) CG:  -189
		Interior loop (394,516) CG; (395,515) UA:  -199
		Interior loop (395,515) UA; (396,514) GC:  -178
		Interior loop (396,514) GC; (397,513) CG:  -300
		Interior loop (397,513) CG; (399,512) AU:   133
		Interior loop (399,512) AU; (400,511) AU:   -69
		Interior loop (400,511) AU; (401,510) AU:   -69
		Interior loop (401,510) AU; (402,509) UA:   -99
		Interior loop (402,509) UA; (403,505) GU:   288
		Interior loop (403,505) GU; (404,504) AU:   -47
		Interior loop (404,504) AU; (405,503) UA:   -99
		Interior loop (405,503) UA; (406,502) AU:   -91
		Interior loop (406,502) AU; (407,501) GC:  -199
		Interior loop (407,501) GC; (409,499) CG:  -195
		Interior loop (409,499) CG; (410,498) GU:   -85
		Interior loop (410,498) GU; (411,497) UA:   -88
		Interior loop (411,497) UA; (423,494) UG:   615
		Interior loop (423,494) UG; (424,493) GC:   -85
		Interior loop (424,493) GC; (425,492) CG:  -300
		Interior loop (425,492) CG; (426,491) AU:  -178
		Interior loop (426,491) AU; (427,490) UG:   -88
		Interior loop (427,490) UG; (428,489) GU:   -70
		Interior loop (428,489) GU; (435,485) AU:   350
		Interior loop (435,485) AU; (436,484) GU:    -1
		Interior loop (436,484) GU; (437,483) CG:  -193
		Interior loop (437,483) CG; (438,482) UA:  -199
		Interior loop (438,482) UA; (439,481) GC:  -178
		Interior loop (439,481) GC; (440,480) AU:  -211
		Interior loop (440,480) AU; (446,476) AU:   363
		Interior loop (446,476) AU; (447,475) AU:   -69
		Interior loop (447,475) AU; (448,474) UG:   -88
		Interior loop (448,474) UG; (449,473) AU:     0
		Interior loop (449,473) AU; (450,472) UA:   -99
		Interior loop (450,472) UA; (451,471) AU:   -91
		Interior loop (451,471) AU; (452,470) UA:   -99
		Interior loop (452,470) UA; (453,469) GC:  -178
		Interior loop (453,469) GC; (454,468) GU:  -127
		Interior loop (454,468) GU; (455,467) CG:  -193
		Interior loop (455,467) CG; (456,466) AU:  -178
		Hairpin  loop (456,466) AU              :   383
		Interior loop (575,670) CG; (576,669) UA:  -199
		Interior loop (576,669) UA; (577,668) AU:   -91
		Interior loop (577,668) AU; (578,667) CG:  -189
		Interior loop (578,667) CG; (579,666) GC:  -203
		Interior loop (579,666) GC; (580,665) UA:  -189
		Interior loop (580,665) UA; (583,663) GU:   315
		Interior loop (583,663) GU; (584,662) UA:   -88
		Interior loop (584,662) UA; (585,661) UA:   -69
		Interior loop (585,661) UA; (586,660) CG:  -211
		Interior loop (586,660) CG; (587,659) AU:  -178
		Interior loop (587,659) AU; (588,658) UA:   -99
		Interior loop (588,658) UA; (589,657) AU:   -91
		Interior loop (589,657) AU; (590,656) UA:   -99
		Interior loop (590,656) UA; (599,655) AU:   393
		Interior loop (599,655) AU; (600,654) GC:  -199
		Interior loop (600,654) GC; (601,653) AU:  -211
		Interior loop (601,653) AU; (602,652) UA:   -99
		Interior loop (602,652) UA; (603,651) GU:     0
		Interior loop (603,651) GU; (604,650) UA:   -88
		Interior loop (604,650) UA; (605,649) GC:  -178
		Interior loop (605,649) GC; (611,644) UA:   192
		Interior loop (611,644) UA; (612,643) UA:   -69
		Interior loop (612,643) UA; (613,642) AU:   -91
		Interior loop (613,642) AU; (614,641) UG:   -88
		Interior loop (614,641) UG; (615,640) UG:   -71
		Interior loop (615,640) UG; (616,638) CG:   184
		Interior loop (616,638) CG; (617,637) CG:  -271
		Interior loop (617,637) CG; (618,636) AU:  -178
		Interior loop (618,636) AU; (619,635) GC:  -199
		Interior loop (619,635) GC; (620,634) GU:  -127
		Interior loop (620,634) GU; (621,633) GU:   -71
		Interior loop (621,633) GU; (622,632) GU:   -71
		Interior loop (622,632) GU; (623,631) AU:   -47
		Interior loop (623,631) AU; (624,630) UA:   -99
		Hairpin  loop (624,630) UA              :   418
		Interior loop (682,807) GC; (683,806) AU:  -211
		Interior loop (683,806) AU; (684,805) AU:   -69
		Interior loop (684,805) AU; (685,804) AU:   -69
		Interior loop (685,804) AU; (686,803) UA:   -99
		Interior loop (686,803) UA; (687,802) UG:   -47
		Interior loop (687,802) UG; (688,801) UA:    -1
		Interior loop (688,801) UA; (689,800) UA:   -69
		Interior loop (689,800) UA; (691,798) AU:   128
		Interior loop (691,798) AU; (692,797) GC:  -199
		Interior loop (692,797) GC; (693,796) GC:  -271
		Interior loop (693,796) GC; (694,795) CG:  -300
		Interior loop (694,795) CG; (696,793) UA:    -4
		Interior loop (696,793) UA; (697,792) GU:     0
		Interior loop (697,792) GU; (698,791) AU:   -47
		Interior loop (707,719) AU; (708,718) GC:  -199
		Interior loop (708,718) GC; (709,717) CG:  -300
		Interior loop (709,717) CG; (710,716) UA:  -199
		Interior loop (710,716) UA; (711,715) UG:   -47
		Hairpin  loop (711,715) UG              :   486
		Interior loop (720,762) GU; (721,761) GC:  -178
		Interior loop (721,761) GC; (722,760) UA:  -189
		Interior loop (722,760) UA; (723,759) AU:   -91
		Interior loop (723,759) AU; (726,756) AU:   184
		Interior loop (726,756) AU; (727,755) AU:   -69
		Interior loop (727,755) AU; (728,754) GC:  -199
		Interior loop (728,754) GC; (731,751) GU:   159
		Interior loop (731,751) GU; (732,750) UA:   -88
		Interior loop (732,750) UA; (733,749) CG:  -211
		Hairpin  loop (733,749) CG              :   445
		Interior loop (764,786) AU; (765,785) AU:   -69
		Interior loop (765,785) AU; (766,784) GC:  -199
		Interior loop (766,784) GC; (767,783) AU:  -211
		Interior loop (767,783) AU; (768,782) AU:   -69
		Interior loop (768,782) AU; (769,781) CG:  -189
		Interior loop (769,781) CG; (770,780) AU:  -178
		Hairpin  loop (770,780) AU              :   463
		Multi    loop (698,791) AU              :   548
		Multi    loop (313,808) CG              :   569
		Interior loop (834,846) UA; (835,845) CG:  -211
		Interior loop (835,845) CG; (836,844) CG:  -271
		Interior loop (836,844) CG; (837,843) UA:  -199
		Hairpin  loop (837,843) UA              :   450
		Multi    loop ( 20,848) UA              :   616
		-144.00`, t)
}

func TestAndronescu2007Params(t *testing.T) {
	compareMFEOutputToViennaRNA("ACGAUCAGAGAUCAGAGCAUACGACAGCAG",
		"..((((...))))...((........))..",
		DefaultTemperature, energy_params.Andronescu2007, DoubleDanglingEnds,
		`External loop                           :  -350
		Interior loop (  3, 13) GC; (  4, 12) AU:  -211
		Interior loop (  4, 12) AU; (  5, 11) UA:   -99
		Interior loop (  5, 11) UA; (  6, 10) CG:  -211
		Hairpin  loop (  6, 10) CG              :   475
		Interior loop ( 17, 28) GC; ( 18, 27) CG:  -300
		Hairpin  loop ( 18, 27) CG              :   324
		-3.72`, t)

	compareMFEOutputToViennaRNA("GGGCUCGUAGAUCAGCGGUAGAUCGCUUCCUUCGCAAGGAAGCCCUGGGUUCAAAUCCCAGCGAGUCCACCA",
		"(((((((..((((.......))))(((((((.....))))))).(((((.......))))))))))))....",
		DefaultTemperature, energy_params.Andronescu2007, DoubleDanglingEnds,
		`External loop                           :  -139
		Interior loop (  1, 68) GC; (  2, 67) GC:  -271
		Interior loop (  2, 67) GC; (  3, 66) GU:  -127
		Interior loop (  3, 66) GU; (  4, 65) CG:  -193
		Interior loop (  4, 65) CG; (  5, 64) UA:  -199
		Interior loop (  5, 64) UA; (  6, 63) CG:  -211
		Interior loop (  6, 63) CG; (  7, 62) GC:  -203
		Interior loop ( 10, 24) GC; ( 11, 23) AU:  -211
		Interior loop ( 11, 23) AU; ( 12, 22) UA:   -99
		Interior loop ( 12, 22) UA; ( 13, 21) CG:  -211
		Hairpin  loop ( 13, 21) CG              :   342
		Interior loop ( 25, 43) GC; ( 26, 42) CG:  -300
		Interior loop ( 26, 42) CG; ( 27, 41) UA:  -199
		Interior loop ( 27, 41) UA; ( 28, 40) UA:   -69
		Interior loop ( 28, 40) UA; ( 29, 39) CG:  -211
		Interior loop ( 29, 39) CG; ( 30, 38) CG:  -271
		Interior loop ( 30, 38) CG; ( 31, 37) UA:  -199
		Hairpin  loop ( 31, 37) UA              :   394
		Interior loop ( 45, 61) CG; ( 46, 60) UA:  -199
		Interior loop ( 46, 60) UA; ( 47, 59) GC:  -178
		Interior loop ( 47, 59) GC; ( 48, 58) GC:  -271
		Interior loop ( 48, 58) GC; ( 49, 57) GC:  -271
		Hairpin  loop ( 49, 57) GC              :   274
		Multi    loop (  7, 62) GC              :  -178
		-32.00`, t)

	compareMFEOutputToViennaRNA(
		"AUGAAACAAUACCAAGAUUUAAUUAAAGACAUUUUUGAAAAUGGUUAUGAAACCGAUGAUCGUACAGGCACAGGAACAAUUGCUCUGUUCGGAUCUAAAUUACGCUGGGAUUUAACUAAAGGUUUUCCUGCGGUAACAACUAAGAAGCUCGCCUGGAAAGCUUGCAUUGCUGAGCUAAUAUGGUUUUUAUCAGGAAGCACAAAUGUCAAUGAUUUACGAUUAAUUCAACACGAUUCGUUAAUCCAAGGCAAAACAGUCUGGGAUGAAAAUUACGAAAAUCAAGCAAAAGAUUUAGGAUACCAUAGCGGUGAACUUGGUCCAAUUUAUGGAAAACAGUGGCGUGAUUUUGGUGGUGUAGACCAAAUUAUAGAAGUUAUUGAUCGUAUUAAAAAACUGCCAAAUGAUAGGCGUCAAAUUGUUUCUGCAUGGAAUCCAGCUGAACUUAAAUAUAUGGCAUUACCGCCUUGUCAUAUGUUCUAUCAGUUUAAUGUGCGUAAUGGCUAUUUGGAUUUGCAGUGGUAUCAACGCUCAGUAGAUGUUUUCUUGGGUCUACCGUUUAAUAUUGCGUCAUAUGCUACGUUAGUUCAUAUUGUAGCUAAGAUGUGUAAUCUUAUUCCAGGGGAUUUGAUAUUUUCUGGUGGUAAUACUCAUAUCUAUAUGAAUCACGUAGAACAAUGUAAAGAAAUUUUGAGGCGUGAACCUAAAGAGCUUUGUGAGCUGGUAAUAAGUGGUCUACCUUAUAAAUUCCGAUAUCUUUCUACUAAAGAACAAUUAAAAUAUGUUCUUAAACUUAGGCCUAAAGAUUUCGUUCUUAACAACUAUGUAUCACACCCUCCUAUUAAAGGAAAGAUGGCGGUGUAA",
		"........(((((...((((.......(((((((.((....((((...((((((...((((..(((((((..........))).))))...))))...((((.(((.((..........(((((((((.(((((............)).)))..)))))))))......)).)))))))..)))))).)))).....)).)))))))..(((.(((....))).))).......((((..((((.((((......)))).)))).......))))(((((........)))))..(((((.((((..((((((....((((.....)))).....((((((((((.((((....(((((((............(((((((.(((........(((((.(((((((((.(((...........((((((......((((((.....(((((((((((.........)))))))))))...))))))...))))))..))).)))))...))))))))).......))).)))))))...........)))))))...)))).))))))))))...((((((..((((((((........(((((((.....((((((((((((((.....))))))))).)))))....))))))))))))))).))))))...........((((((((.((((.(((........(((((...)))))((((..(((..(((...............)))..)))..)))).(((((((.........)))))))....))).)))).))))))))))))...))..)))))))))......((((.....)))).))))...)))))..",
		DefaultTemperature, energy_params.Andronescu2007, DoubleDanglingEnds,
		`External loop                           :   -16
			Interior loop (  9,859) AU; ( 10,858) UG:   -88
			Interior loop ( 10,858) UG; ( 11,857) AU:     0
			Interior loop ( 11,857) AU; ( 12,856) CG:  -189
			Interior loop ( 12,856) CG; ( 13,855) CG:  -271
			Interior loop ( 13,855) CG; ( 17,851) AU:   166
			Interior loop ( 17,851) AU; ( 18,850) UA:   -99
			Interior loop ( 18,850) UA; ( 19,849) UG:   -47
			Interior loop ( 19,849) UG; ( 20,848) UA:    -1
			Interior loop ( 28,207) GC; ( 29,206) AU:  -211
			Interior loop ( 29,206) AU; ( 30,205) CG:  -189
			Interior loop ( 30,205) CG; ( 31,204) AU:  -178
			Interior loop ( 31,204) AU; ( 32,203) UA:   -99
			Interior loop ( 32,203) UA; ( 33,202) UA:   -69
			Interior loop ( 33,202) UA; ( 34,201) UA:   -69
			Interior loop ( 34,201) UA; ( 36,199) UA:   128
			Interior loop ( 36,199) UA; ( 37,198) GC:  -178
			Interior loop ( 37,198) GC; ( 42,192) UA:   108
			Interior loop ( 42,192) UA; ( 43,191) GC:  -178
			Interior loop ( 43,191) GC; ( 44,190) GU:  -127
			Interior loop ( 44,190) GU; ( 45,189) UA:   -88
			Interior loop ( 45,189) UA; ( 49,187) GU:   166
			Interior loop ( 49,187) GU; ( 50,186) AU:   -47
			Interior loop ( 50,186) AU; ( 51,185) AU:   -69
			Interior loop ( 51,185) AU; ( 52,184) AU:   -69
			Interior loop ( 52,184) AU; ( 53,183) CG:  -189
			Interior loop ( 53,183) CG; ( 54,182) CG:  -271
			Interior loop ( 58, 95) GC; ( 59, 94) AU:  -211
			Interior loop ( 59, 94) AU; ( 60, 93) UA:   -99
			Interior loop ( 60, 93) UA; ( 61, 92) CG:  -211
			Interior loop ( 61, 92) CG; ( 64, 88) AU:   167
			Interior loop ( 64, 88) AU; ( 65, 87) CG:  -189
			Interior loop ( 65, 87) CG; ( 66, 86) AU:  -178
			Interior loop ( 66, 86) AU; ( 67, 85) GC:  -199
			Interior loop ( 67, 85) GC; ( 68, 83) GC:    40
			Interior loop ( 68, 83) GC; ( 69, 82) CG:  -300
			Interior loop ( 69, 82) CG; ( 70, 81) AU:  -178
			Hairpin  loop ( 70, 81) AU              :   483
			Interior loop ( 99,179) AU; (100,178) UA:   -99
			Interior loop (100,178) UA; (101,177) UA:   -69
			Interior loop (101,177) UA; (102,176) AU:   -91
			Interior loop (102,176) AU; (104,175) GC:   112
			Interior loop (104,175) GC; (105,174) CG:  -300
			Interior loop (105,174) CG; (106,173) UA:  -199
			Interior loop (106,173) UA; (108,171) GU:    62
			Interior loop (108,171) GU; (109,170) GC:  -178
			Interior loop (109,170) GC; (120,163) AU:   320
			Interior loop (120,163) AU; (121,162) GU:    -1
			Interior loop (121,162) GU; (122,161) GC:  -178
			Interior loop (122,161) GC; (123,160) UG:  -193
			Interior loop (123,160) UG; (124,159) UA:    -1
			Interior loop (124,159) UA; (125,158) UA:   -69
			Interior loop (125,158) UA; (126,157) UA:   -69
			Interior loop (126,157) UA; (127,156) CG:  -211
			Interior loop (127,156) CG; (128,155) CG:  -271
			Interior loop (128,155) CG; (130,152) GC:   175
			Interior loop (130,152) GC; (131,151) CG:  -300
			Interior loop (131,151) CG; (132,150) GC:  -203
			Interior loop (132,150) GC; (133,148) GC:    40
			Interior loop (133,148) GC; (134,147) UG:  -193
			Hairpin  loop (134,147) UG              :   503
			Multi    loop ( 54,182) CG              :   106
			Interior loop (210,227) UA; (211,226) GC:  -178
			Interior loop (211,226) GC; (212,225) AU:  -211
			Interior loop (212,225) AU; (214,223) UA:   137
			Interior loop (214,223) UA; (215,222) UA:   -69
			Interior loop (215,222) UA; (216,221) AU:   -91
			Hairpin  loop (216,221) AU              :   464
			Interior loop (235,275) UA; (236,274) CG:  -211
			Interior loop (236,274) CG; (237,273) GC:  -203
			Interior loop (237,273) GC; (238,272) UA:  -189
			Interior loop (238,272) UA; (241,264) AU:   382
			Interior loop (241,264) AU; (242,263) UA:   -99
			Interior loop (242,263) UA; (243,262) CG:  -211
			Interior loop (243,262) CG; (244,261) CG:  -271
			Interior loop (244,261) CG; (246,259) AU:    62
			Interior loop (246,259) AU; (247,258) GC:  -199
			Interior loop (247,258) GC; (248,257) GU:  -127
			Interior loop (248,257) GU; (249,256) CG:  -193
			Hairpin  loop (249,256) CG              :   338
			Interior loop (276,293) AU; (277,292) AU:   -69
			Interior loop (277,292) AU; (278,291) AU:   -69
			Interior loop (278,291) AU; (279,290) UA:   -99
			Interior loop (279,290) UA; (280,289) CG:  -211
			Hairpin  loop (280,289) CG              :   324
			Interior loop (296,827) GC; (297,826) AU:  -211
			Interior loop (297,826) AU; (298,825) UA:   -99
			Interior loop (298,825) UA; (299,824) AU:   -91
			Interior loop (299,824) AU; (300,823) CG:  -189
			Interior loop (300,823) CG; (302,822) AU:   133
			Interior loop (302,822) AU; (303,821) UA:   -99
			Interior loop (303,821) UA; (304,820) AU:   -91
			Interior loop (304,820) AU; (305,819) GC:  -199
			Interior loop (305,819) GC; (308,816) GC:    22
			Interior loop (308,816) GC; (309,815) UA:  -189
			Interior loop (309,815) UA; (310,811) GC:   277
			Interior loop (310,811) GC; (311,810) AU:  -211
			Interior loop (311,810) AU; (312,809) AU:   -69
			Interior loop (312,809) AU; (313,808) CG:  -189
			Interior loop (318,330) UA; (319,329) CG:  -211
			Interior loop (319,329) CG; (320,328) CG:  -271
			Interior loop (320,328) CG; (321,327) AU:  -178
			Hairpin  loop (321,327) AU              :   447
			Interior loop (336,571) GU; (337,570) UA:   -88
			Interior loop (337,570) UA; (338,569) GC:  -178
			Interior loop (338,569) GC; (339,568) GU:  -127
			Interior loop (339,568) GU; (340,567) CG:  -193
			Interior loop (340,567) CG; (341,566) GC:  -203
			Interior loop (341,566) GC; (342,565) UG:  -193
			Interior loop (342,565) UG; (343,564) GU:   -70
			Interior loop (343,564) GU; (344,563) AU:   -47
			Interior loop (344,563) AU; (345,562) UA:   -99
			Interior loop (345,562) UA; (347,560) UA:   128
			Interior loop (347,560) UA; (348,559) UA:   -69
			Interior loop (348,559) UA; (349,558) GU:     0
			Interior loop (349,558) GU; (350,557) GU:   -71
			Interior loop (350,557) GU; (355,553) GC:   165
			Interior loop (355,553) GC; (356,552) UA:  -189
			Interior loop (356,552) UA; (357,551) AU:   -91
			Interior loop (357,551) AU; (358,550) GC:  -199
			Interior loop (358,550) GC; (359,549) AU:  -211
			Interior loop (359,549) AU; (360,548) CG:  -189
			Interior loop (360,548) CG; (361,547) CG:  -271
			Interior loop (361,547) CG; (374,535) UG:   209
			Interior loop (374,535) UG; (375,534) UA:    -1
			Interior loop (375,534) UA; (376,533) AU:   -91
			Interior loop (376,533) AU; (377,532) UG:   -88
			Interior loop (377,532) UG; (378,531) UA:    -1
			Interior loop (378,531) UA; (379,530) GC:  -178
			Interior loop (379,530) GC; (380,529) AU:  -211
			Interior loop (380,529) AU; (382,527) CG:    62
			Interior loop (382,527) CG; (383,526) GC:  -203
			Interior loop (383,526) GC; (384,525) UA:  -189
			Interior loop (384,525) UA; (393,517) AU:   305
			Interior loop (393,517) AU; (394,516) CG:  -189
			Interior loop (394,516) CG; (395,515) UA:  -199
			Interior loop (395,515) UA; (396,514) GC:  -178
			Interior loop (396,514) GC; (397,513) CG:  -300
			Interior loop (397,513) CG; (399,512) AU:   133
			Interior loop (399,512) AU; (400,511) AU:   -69
			Interior loop (400,511) AU; (401,510) AU:   -69
			Interior loop (401,510) AU; (402,509) UA:   -99
			Interior loop (402,509) UA; (403,505) GU:   288
			Interior loop (403,505) GU; (404,504) AU:   -47
			Interior loop (404,504) AU; (405,503) UA:   -99
			Interior loop (405,503) UA; (406,502) AU:   -91
			Interior loop (406,502) AU; (407,501) GC:  -199
			Interior loop (407,501) GC; (409,499) CG:  -195
			Interior loop (409,499) CG; (410,498) GU:   -85
			Interior loop (410,498) GU; (411,497) UA:   -88
			Interior loop (411,497) UA; (423,494) UG:   615
			Interior loop (423,494) UG; (424,493) GC:   -85
			Interior loop (424,493) GC; (425,492) CG:  -300
			Interior loop (425,492) CG; (426,491) AU:  -178
			Interior loop (426,491) AU; (427,490) UG:   -88
			Interior loop (427,490) UG; (428,489) GU:   -70
			Interior loop (428,489) GU; (435,485) AU:   350
			Interior loop (435,485) AU; (436,484) GU:    -1
			Interior loop (436,484) GU; (437,483) CG:  -193
			Interior loop (437,483) CG; (438,482) UA:  -199
			Interior loop (438,482) UA; (439,481) GC:  -178
			Interior loop (439,481) GC; (440,480) AU:  -211
			Interior loop (440,480) AU; (446,476) AU:   363
			Interior loop (446,476) AU; (447,475) AU:   -69
			Interior loop (447,475) AU; (448,474) UG:   -88
			Interior loop (448,474) UG; (449,473) AU:     0
			Interior loop (449,473) AU; (450,472) UA:   -99
			Interior loop (450,472) UA; (451,471) AU:   -91
			Interior loop (451,471) AU; (452,470) UA:   -99
			Interior loop (452,470) UA; (453,469) GC:  -178
			Interior loop (453,469) GC; (454,468) GU:  -127
			Interior loop (454,468) GU; (455,467) CG:  -193
			Interior loop (455,467) CG; (456,466) AU:  -178
			Hairpin  loop (456,466) AU              :   383
			Interior loop (575,670) CG; (576,669) UA:  -199
			Interior loop (576,669) UA; (577,668) AU:   -91
			Interior loop (577,668) AU; (578,667) CG:  -189
			Interior loop (578,667) CG; (579,666) GC:  -203
			Interior loop (579,666) GC; (580,665) UA:  -189
			Interior loop (580,665) UA; (583,663) GU:   315
			Interior loop (583,663) GU; (584,662) UA:   -88
			Interior loop (584,662) UA; (585,661) UA:   -69
			Interior loop (585,661) UA; (586,660) CG:  -211
			Interior loop (586,660) CG; (587,659) AU:  -178
			Interior loop (587,659) AU; (588,658) UA:   -99
			Interior loop (588,658) UA; (589,657) AU:   -91
			Interior loop (589,657) AU; (590,656) UA:   -99
			Interior loop (590,656) UA; (599,655) AU:   393
			Interior loop (599,655) AU; (600,654) GC:  -199
			Interior loop (600,654) GC; (601,653) AU:  -211
			Interior loop (601,653) AU; (602,652) UA:   -99
			Interior loop (602,652) UA; (603,651) GU:     0
			Interior loop (603,651) GU; (604,650) UA:   -88
			Interior loop (604,650) UA; (605,649) GC:  -178
			Interior loop (605,649) GC; (611,644) UA:   192
			Interior loop (611,644) UA; (612,643) UA:   -69
			Interior loop (612,643) UA; (613,642) AU:   -91
			Interior loop (613,642) AU; (614,641) UG:   -88
			Interior loop (614,641) UG; (615,640) UG:   -71
			Interior loop (615,640) UG; (616,638) CG:   184
			Interior loop (616,638) CG; (617,637) CG:  -271
			Interior loop (617,637) CG; (618,636) AU:  -178
			Interior loop (618,636) AU; (619,635) GC:  -199
			Interior loop (619,635) GC; (620,634) GU:  -127
			Interior loop (620,634) GU; (621,633) GU:   -71
			Interior loop (621,633) GU; (622,632) GU:   -71
			Interior loop (622,632) GU; (623,631) AU:   -47
			Interior loop (623,631) AU; (624,630) UA:   -99
			Hairpin  loop (624,630) UA              :   418
			Interior loop (682,807) GC; (683,806) AU:  -211
			Interior loop (683,806) AU; (684,805) AU:   -69
			Interior loop (684,805) AU; (685,804) AU:   -69
			Interior loop (685,804) AU; (686,803) UA:   -99
			Interior loop (686,803) UA; (687,802) UG:   -47
			Interior loop (687,802) UG; (688,801) UA:    -1
			Interior loop (688,801) UA; (689,800) UA:   -69
			Interior loop (689,800) UA; (691,798) AU:   128
			Interior loop (691,798) AU; (692,797) GC:  -199
			Interior loop (692,797) GC; (693,796) GC:  -271
			Interior loop (693,796) GC; (694,795) CG:  -300
			Interior loop (694,795) CG; (696,793) UA:    -4
			Interior loop (696,793) UA; (697,792) GU:     0
			Interior loop (697,792) GU; (698,791) AU:   -47
			Interior loop (707,719) AU; (708,718) GC:  -199
			Interior loop (708,718) GC; (709,717) CG:  -300
			Interior loop (709,717) CG; (710,716) UA:  -199
			Interior loop (710,716) UA; (711,715) UG:   -47
			Hairpin  loop (711,715) UG              :   486
			Interior loop (720,762) GU; (721,761) GC:  -178
			Interior loop (721,761) GC; (722,760) UA:  -189
			Interior loop (722,760) UA; (723,759) AU:   -91
			Interior loop (723,759) AU; (726,756) AU:   184
			Interior loop (726,756) AU; (727,755) AU:   -69
			Interior loop (727,755) AU; (728,754) GC:  -199
			Interior loop (728,754) GC; (731,751) GU:   159
			Interior loop (731,751) GU; (732,750) UA:   -88
			Interior loop (732,750) UA; (733,749) CG:  -211
			Hairpin  loop (733,749) CG              :   445
			Interior loop (764,786) AU; (765,785) AU:   -69
			Interior loop (765,785) AU; (766,784) GC:  -199
			Interior loop (766,784) GC; (767,783) AU:  -211
			Interior loop (767,783) AU; (768,782) AU:   -69
			Interior loop (768,782) AU; (769,781) CG:  -189
			Interior loop (769,781) CG; (770,780) AU:  -178
			Hairpin  loop (770,780) AU              :   463
			Multi    loop (698,791) AU              :   310
			Multi    loop (313,808) CG              :    45
			Interior loop (834,846) UA; (835,845) CG:  -211
			Interior loop (835,845) CG; (836,844) CG:  -271
			Interior loop (836,844) CG; (837,843) UA:  -199
			Hairpin  loop (837,843) UA              :   450
			Multi    loop ( 20,848) UA              :    89
			-161.02`, t)
}

func TestMinimumFreeEnergy(t *testing.T) {
	compareMFEOutputToViennaRNA("ACGAUCAGAGAUCAGAGCAUACGACAGCAG",
		"..((((...))))...((........))..",
		DefaultTemperature, energy_params.Turner2004, DoubleDanglingEnds,
		`External loop                           :  -300
		Interior loop (  3, 13) GC; (  4, 12) AU:  -240
		Interior loop (  4, 12) AU; (  5, 11) UA:  -110
		Interior loop (  5, 11) UA; (  6, 10) CG:  -240
		Hairpin  loop (  6, 10) CG              :   540
		Interior loop ( 17, 28) GC; ( 18, 27) CG:  -340
		Hairpin  loop ( 18, 27) CG              :   400
		-2.90`, t)

	compareMFEOutputToViennaRNA("ACGAUCAGAGAUCAGAGCAUACGACAGCAG",
		"..((((...))))...((........))..",
		4, energy_params.Turner2004, DoubleDanglingEnds,
		`External loop                           :  -423
		Interior loop (  3, 13) GC; (  4, 12) AU:  -346
		Interior loop (  4, 12) AU; (  5, 11) UA:  -198
		Interior loop (  5, 11) UA; (  6, 10) CG:  -346
		Hairpin  loop (  6, 10) CG              :   496
		Interior loop ( 17, 28) GC; ( 18, 27) CG:  -462
		Hairpin  loop ( 18, 27) CG              :   230
		-10.49`, t)

	compareMFEOutputToViennaRNA("AAAACGGUCCUUAUCAGGACCAAACA",
		".....((((((....)))))).....",
		DefaultTemperature, energy_params.Turner2004, DoubleDanglingEnds,
		`External loop                           :  -150
		Interior loop (  6, 21) GC; (  7, 20) GC:  -330
		Interior loop (  7, 20) GC; (  8, 19) UA:  -220
		Interior loop (  8, 19) UA; (  9, 18) CG:  -240
		Interior loop (  9, 18) CG; ( 10, 17) CG:  -330
		Interior loop ( 10, 17) CG; ( 11, 16) UA:  -210
		Hairpin  loop ( 11, 16) UA              :   550
		-9.30
		`, t)

	compareMFEOutputToViennaRNA("AUUCUUGCUUCAACAGUGUUUGAACGGAAU",
		"..............................",
		DefaultTemperature, energy_params.Turner2004, DoubleDanglingEnds,
		`External loop                           :     0
		0.00`, t)

	compareMFEOutputToViennaRNA("UCGGCCACAAACACACAAUCUACUGUUGGUCGA",
		"(((((((...................)))))))",
		DefaultTemperature, energy_params.Turner2004, DoubleDanglingEnds,
		`External loop                           :    50
		Interior loop (  1, 33) UA; (  2, 32) CG:  -240
		Interior loop (  2, 32) CG; (  3, 31) GC:  -240
		Interior loop (  3, 31) GC; (  4, 30) GU:  -150
		Interior loop (  4, 30) GU; (  5, 29) CG:  -250
		Interior loop (  5, 29) CG; (  6, 28) CG:  -330
		Interior loop (  6, 28) CG; (  7, 27) AU:  -210
		Hairpin  loop (  7, 27) AU              :   700
		-6.70`, t)

	compareMFEOutputToViennaRNA("UCGGCCACAAACACACAAUCUACUGUUGGUCGA",
		"(((((((...................)))))))",
		58.0, energy_params.Turner2004, DoubleDanglingEnds,
		`External loop                           :    28
		Interior loop (  1, 33) UA; (  2, 32) CG:  -172
		Interior loop (  2, 32) CG; (  3, 31) GC:  -184
		Interior loop (  3, 31) GC; (  4, 30) GU:  -103
		Interior loop (  4, 30) GU; (  5, 29) CG:  -181
		Interior loop (  5, 29) CG; (  6, 28) CG:  -261
		Interior loop (  6, 28) CG; (  7, 27) AU:  -153
		Hairpin  loop (  7, 27) AU              :   690
		-3.36`, t)

	compareMFEOutputToViennaRNA("GUUUUUAUCUUACACACGCUUGUGUAAGAUAGUUA",
		".....(((((((((((....)))))))))))....",
		DefaultTemperature, energy_params.Turner2004, DoubleDanglingEnds,
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
		-12.80`, t)

	compareMFEOutputToViennaRNA("GGGCUCGUAGAUCAGCGGUAGAUCGCUUCCUUCGCAAGGAAGCCCUGGGUUCAAAUCCCAGCGAGUCCACCA",
		"(((((((..((((.......))))(((((((.....))))))).(((((.......))))))))))))....",
		DefaultTemperature, energy_params.Turner2004, DoubleDanglingEnds,
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
		-31.00`, t)

	compareMFEOutputToViennaRNA("GGGCUCGUAGAUCAGCGGUAGAUCGCUUCCUUCGCAAGGAAGCCCUGGGUUCAAAUCCCAGCGAGUCCACCA",
		"(((((((..((((.......))))(((((((.....))))))).(((((.......))))))))))))....",
		22.19, energy_params.Turner2004, DoubleDanglingEnds,
		`External loop                           :  -204
		Interior loop (  1, 68) GC; (  2, 67) GC:  -378
		Interior loop (  2, 67) GC; (  3, 66) GU:  -182
		Interior loop (  3, 66) GU; (  4, 65) CG:  -298
		Interior loop (  4, 65) CG; (  5, 64) UA:  -250
		Interior loop (  5, 64) UA; (  6, 63) CG:  -287
		Interior loop (  6, 63) CG; (  7, 62) GC:  -279
		Interior loop ( 10, 24) GC; ( 11, 23) AU:  -287
		Interior loop ( 11, 23) AU; ( 12, 22) UA:  -149
		Interior loop ( 12, 22) UA; ( 13, 21) CG:  -287
		Hairpin  loop ( 13, 21) CG              :   391
		Interior loop ( 25, 43) GC; ( 26, 42) CG:  -394
		Interior loop ( 26, 42) CG; ( 27, 41) UA:  -250
		Interior loop ( 27, 41) UA; ( 28, 40) UA:  -118
		Interior loop ( 28, 40) UA; ( 29, 39) CG:  -287
		Interior loop ( 29, 39) CG; ( 30, 38) CG:  -378
		Interior loop ( 30, 38) CG; ( 31, 37) UA:  -250
		Hairpin  loop ( 31, 37) UA              :   538
		Interior loop ( 45, 61) CG; ( 46, 60) UA:  -250
		Interior loop ( 46, 60) UA; ( 47, 59) GC:  -249
		Interior loop ( 47, 59) GC; ( 48, 58) GC:  -378
		Interior loop ( 48, 58) GC; ( 49, 57) GC:  -378
		Hairpin  loop ( 49, 57) GC              :   374
		Multi    loop (  7, 62) GC              :   149
		-40.81`, t)

	compareMFEOutputToViennaRNA(
		"AUGAAACAAUACCAAGAUUUAAUUAAAGACAUUUUUGAAAAUGGUUAUGAAACCGAUGAUCGUACAGGCACAGGAACAAUUGCUCUGUUCGGAUCUAAAUUACGCUGGGAUUUAACUAAAGGUUUUCCUGCGGUAACAACUAAGAAGCUCGCCUGGAAAGCUUGCAUUGCUGAGCUAAUAUGGUUUUUAUCAGGAAGCACAAAUGUCAAUGAUUUACGAUUAAUUCAACACGAUUCGUUAAUCCAAGGCAAAACAGUCUGGGAUGAAAAUUACGAAAAUCAAGCAAAAGAUUUAGGAUACCAUAGCGGUGAACUUGGUCCAAUUUAUGGAAAACAGUGGCGUGAUUUUGGUGGUGUAGACCAAAUUAUAGAAGUUAUUGAUCGUAUUAAAAAACUGCCAAAUGAUAGGCGUCAAAUUGUUUCUGCAUGGAAUCCAGCUGAACUUAAAUAUAUGGCAUUACCGCCUUGUCAUAUGUUCUAUCAGUUUAAUGUGCGUAAUGGCUAUUUGGAUUUGCAGUGGUAUCAACGCUCAGUAGAUGUUUUCUUGGGUCUACCGUUUAAUAUUGCGUCAUAUGCUACGUUAGUUCAUAUUGUAGCUAAGAUGUGUAAUCUUAUUCCAGGGGAUUUGAUAUUUUCUGGUGGUAAUACUCAUAUCUAUAUGAAUCACGUAGAACAAUGUAAAGAAAUUUUGAGGCGUGAACCUAAAGAGCUUUGUGAGCUGGUAAUAAGUGGUCUACCUUAUAAAUUCCGAUAUCUUUCUACUAAAGAACAAUUAAAAUAUGUUCUUAAACUUAGGCCUAAAGAUUUCGUUCUUAACAACUAUGUAUCACACCCUCCUAUUAAAGGAAAGAUGGCGGUGUAA",
		"........(((((...((((.......(((((((.((....((((...((((((...((((..(((((((..........))).))))...))))...((((.(((.((..........(((((((((.(((((............)).)))..)))))))))......)).)))))))..)))))).)))).....)).)))))))..(((.(((....))).))).......((((..((((.((((......)))).)))).......))))(((((........)))))..(((((.((((..((((((....((((.....)))).....((((((((((.((((....(((((((............(((((((.(((........(((((.(((((((((.(((...........((((((......((((((.....(((((((((((.........)))))))))))...))))))...))))))..))).)))))...))))))))).......))).)))))))...........)))))))...)))).))))))))))...((((((..((((((((........(((((((.....((((((((((((((.....))))))))).)))))....))))))))))))))).))))))...........((((((((.((((.(((........(((((...)))))((((..(((..(((...............)))..)))..)))).(((((((.........)))))))....))).)))).))))))))))))...))..)))))))))......((((.....)))).))))...)))))..",
		DefaultTemperature, energy_params.Turner2004, DoubleDanglingEnds,
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
		-197.70`, t)

	compareMFEOutputToViennaRNA(
		"GAGAUACCUACAGCGUGAGCUAUGAGAAAGCGCCACGCUUCCCGAAGGGAGAAAGGCGGACAGGUAUCCGGUAAGCGGCAGGGUCGGAACAGGAGAGCGCACGAGGGAGCUUCCAGGGGGAAACGCCUGGUAUCUUUAUAGUCCUGUCGGGUUUCGCCACCUCUGACUUGAGCGUCGAUUUUUGUGAUGCUCGUCAGGGGGGCGGAGCCUAUGGAAAAACGCCAGCAACGCGGCCUUUUUACGGUUCCUGGCCUUUUGCUGGCCUUUUGCUCACAUGUUCUUUCCUGCGUUAUCCCCUGAUUCUGUGGAUAACCGUAUUACCGCCUUUGAGUGAGCUGAUACCGCUCGCCGCAGCCGAACGACCGAGCGCAGCGAGUCAGUGAGCGAGGAAGCGGAAGAGCGCCCAAUACGCAAACCGCCUCUCCCCGCGCGUUGGCCGAUUCAUUAAUGCAGCUGGCACGACAGGUUUCCCGACUGGAAAGCGGGCAGUGAGCGCAACGCAAUUAAUGUGAGUUAGCUCACUCAUUAGGCACCCCAGGCUUUACACUUUAUGCUUCCGGCUCGUAUGUUGUGUGGAAUUGUGAGCGGAUAACAAUUUCACACAGGAAACAGCUAUGACCAUGAUUACGCCAAGCUUGCAUGCCUGCAGGUCGACUCUAGAGGAUCCCCGGGUACCGAGCUCGAAUUCACUGGCCGUCGUUUUACAACGUCGUGACUGGGAAAACCCUGGCGUUACCCAACUUAAUCGCCUUGCAGCACAUCCCCCUUUCGCCAGCUGGCGUAAUAGCGAAGAGGCCCGCACCGAUCGCCCUUCCCAACAGUUGCGCAGCCUGAAUGGCGAAUGGCGCCUGAUGCGGUAUUUUCUCCUUACGCAUCUGUGCGGUAUUUCACACCGCAUAUGGUGCACUCUCAGUACAAUCUGCUCUGAUGCCGCAUAGUUAAGCCAGCCCCGACACCCGCCAACACCCGCUGACGCGCCCUGACGGGCUUGUCUGCUCCCGGCAUCCGCUUACAGACAAGCUGUGACCGUCUCCGGGAGCUGCAUGUGUCAGAGGUUUUCACCGUCAUCACCGAAACGCGCGAGACGAAAGGGCCUCGUGAUACGCCUAUUUUUAUAGGUUAAUGUCAUGAUAAUAAUGGUUUCUUAGACGUCAGGUGGCACUUUUCGGGGAAAUGUGCGCGGAACCCCUAUUUGUUUAUUUUUCUAAAUACAUUCAAAUAUGUAUCCGCUCAUGAGACAAUAACCCUGAUAAAUGCUUCAAUAAUAUUGAAAAAGGAAGAGUAUGAGUAUUCAACAUUUCCGUGUCGCCCUUAUUCCCUUUUUUGCGGCAUUUUGCCUUCCUGUUUUUGCUCACCCAGAAACGCUGGUGAAAGUAAAAGAUGCUGAAGAUCAGUUGGGUGCACGAGUGGGUUACAUCGAACUGGAUCUCAACAGCGGUAAGAUCCUUGAGAGUUUUCGCCCCGAAGAACGUUUUCCAAUGAUGAGCACUUUUAAAGUUCUGCUAUGUGGCGCGGUAUUAUCCCGUAUUGACGCCGGGCAAGAGCAACUCGGUCGCCGCAUACACUAUUCUCAGAAUGACUUGGUUGAGUACUCACCAGUCACAGAAAAGCAUCUUACGGAUGGCAUGACAGUAAGAGAAUUAUGCAGUGCUGCCAUAACCAUGAGUGAUAACACUGCGGCCAACUUACUUCUGACAACGAUCGGAGGACCGAAGGAGCUAACCGCUUUUUUGCACAACAUGGGGGAUCAUGUAACUCGCCUUGAUCGUUGGGAACCGGAGCUGAAUGAAGCCAUACCAAACGACGAGCGUGACACCACGAUGCCUGUAGCAAUGGCAACAACGUUGCGCAAACUAUUAACUGGCGAACUACUUACUCUAGCUUCCCGGCAACAAUUAAUAGACUGGAUGGAGGCGGAUAAAGUUGCAGGACCACUUCUGCGCUCGGCCCUUCCGGCUGGCUGGUUUAUUGCUGAUAAAUCUGGAGCCGGUGAGCGUGGGUCUCGCGGUAUCAUUGCAGCACUGGGGCCAGAUGGUAAGCCCUCCCGUAUCGUAGUUAUCUACACGACGGGGAGUCAGGCAACUAUGGAUGAACGAAAUAGACAGAUCGCUGAGAUAGGUGCCUCACUGAUUAAGCAUUGGUAACUGUCAGACCAAGUUUACUCAUAUAUACUUUAGAUUGAUUUAAAACUUCAUUUUUAAUUUAAAAGGAUCUAGGUGAAGAUCCUUUUUGAUAAUCUCAUGACCAAAAUCCCUUAACGUGAGUUUUCGUUCCACUGAGCGUCAGACCCCGUAGAAAAGAUCAAAGGAUCUUCUUGAGAUCCUUUUUUUCUGCGCGUAAUCUGCUGCUUGCAAACAAAAAAACCACCGCUACCAGCGGUGGUUUGUUUGCCGGAUCAAGAGCUACCAACUCUUUUUCCGAAGGUAACUGGCUUCAGCAGAGCGCAGAUACCAAAUACUGUUCUUCUAGUGUAGCCGUAGUUAGGCCACCACUUCAAGAACUCUGUAGCACCGCCUACAUACCUCGCUCUGCUAAUCCUGUUACCAGUGGCUGCUGCCAGUGGCGAUAAGUCGUGUCUUACCGGGUUGGACUCAAGACGAUAGUUACCGGAUAAGGCGCAGCGGUCGGGCUGAACGGGGGGUUCGUGCACACAGCCCAGCUUGGAGCGAACGACCUACACCGAACU",
		"(((((......((((((.(((.......)))..))))))((((...)))).((((((((...((((..((((.....((((((((((....(((.((((((.((((((((..(((((........)))))...........(((....((((((((((.((((((((..((((((((.......))))))))))))))))))))))))))..)))((((.((((((((...((((..............))))..)))))))).))))........)))))))).)))))).))).)))))))))).....)))))))).))))))))((((((((((((...((((((.((.((((......)).)))).))))))......(((.(((.((((....(((.......)))...))))...))).)))((((((..((.((((((......(((.((....(((.........)))....)).)))))))))))))))))...........))))))))))))..((((((...(((((((......(((((.((..((((.((...((((((((((((((.........))))))))))))))...))))))..)).))))).......)))))))..)))))).(((((((.(((((((......(((((.....))))).......(((((((.((.(((..((.(((((...(((....))).(((((((.......))).))))..(((((((.....((.(((((((((((.((((((...((((.(((.((...))...)))))))......)))))))))).....))))))).))(((.....)))((((.(((((.(((((((((((((((((.......)))))))).))))).....(((......)))....(((((((((.((............((((((..((........)).((((.((((....)))).)))).((((((((....((.((((((.....)))))).))...)))))))).((((.((((((.((((......(((.((.((((((.((((..((((.(((((((((((((((.((((((....))))))...)))))))).......)).)))))...))))..)).))...))))))......(((.(((((((...((((((..(((((....)))))....)))))).)).))))).))))))))...)))))))))).))))((((((...))))))...(((((.((.(((....))))).))))))))))).............)).))))))))).((((.(((.....(((((((((.((.....)))))))..))))((((((((..(((..(((((((((..((((((((((..((((..(((........)))))))..)))))((((((((.((((....(((((...)))))......)))).)))))))).((((((((...(((((((((((......)))))).)))))))).))))).)))))..)))).))...)))..))).....(((((.(((.((....))))))))))......))))))))..))).)))).....)))))))))..)))).)))))))......)))))))..))).)).))))))).(((((((((..(((((((((.(((..(((((((((.....)))))))))....((((((....))))))......))))))))))))((..((((.(((((((((.(((..............(((.(((.(((((.....(((.((((((.......)))))))))..(((((((((((.((((((.......))).)))))))......))))))).(((((.((.(((.((....((.((((((....)))))))))).))))))))))((((((((((((((...))))))))...))))))....)))))))).)))))).))))).)))).))))..))))).))))))...((((((...((((....))))...))))))(((.(((((.((((.......(((..........))).....)))).))))).)))((.(((((.(((((.........)))))))))).)).......))))))))))))))...............(((.((((((((((......)))))))))).))).)))))...................(((.((((((((((((..(((((.((((...(((((((((((....((((((((...))))))))))))))))))).....)))).)))))((((((((....(((((((((...)))))))))))))))))((((..(((((.......)))))..))))..((((((.((.((.(((((((((.((............((((((..((((..(((.......)))...))))..))))))..(((((......)))))...))))))))))))).)).))))))...(((((((((....(((.....)))..((((((((((((((.........)))).....)))).)))))))))))))))((((((.((((.....)))).....))))))....))))))))).))).)))........",
		DefaultTemperature, energy_params.Turner2004, DoubleDanglingEnds,
		`External loop                           :  -320
		Interior loop (  1,2240) GC; (  2,2239) AU:  -240
		Interior loop (  2,2239) AU; (  3,2238) GC:  -210
		Interior loop (  3,2238) GC; (  4,2237) AU:  -240
		Interior loop (  4,2237) AU; (  5,2236) UA:  -110
		Interior loop ( 12, 39) AU; ( 13, 38) GC:  -210
		Interior loop ( 13, 38) GC; ( 14, 37) CG:  -340
		Interior loop ( 14, 37) CG; ( 15, 36) GC:  -240
		Interior loop ( 15, 36) GC; ( 16, 35) UA:  -220
		Interior loop ( 16, 35) UA; ( 17, 34) GC:  -210
		Interior loop ( 17, 34) GC; ( 19, 31) GC:   110
		Interior loop ( 19, 31) GC; ( 20, 30) CG:  -340
		Interior loop ( 20, 30) CG; ( 21, 29) UA:  -210
		Hairpin  loop ( 21, 29) UA              :   550
		Interior loop ( 40, 50) UA; ( 41, 49) CG:  -240
		Interior loop ( 41, 49) CG; ( 42, 48) CG:  -330
		Interior loop ( 42, 48) CG; ( 43, 47) CG:  -330
		Hairpin  loop ( 43, 47) CG              :   540
		Interior loop ( 52,328) AU; ( 53,327) AU:   -90
		Interior loop ( 53,327) AU; ( 54,326) AU:   -90
		Interior loop ( 54,326) AU; ( 55,325) GC:  -210
		Interior loop ( 55,325) GC; ( 56,324) GC:  -330
		Interior loop ( 56,324) GC; ( 57,323) CG:  -340
		Interior loop ( 57,323) CG; ( 58,322) GC:  -240
		Interior loop ( 58,322) GC; ( 59,321) GC:  -330
		Interior loop ( 59,321) GC; ( 63,319) GU:   300
		Interior loop ( 63,319) GU; ( 64,318) GU:   -50
		Interior loop ( 64,318) GU; ( 65,317) UA:  -140
		Interior loop ( 65,317) UA; ( 66,316) AU:  -130
		Interior loop ( 66,316) AU; ( 69,315) CG:   330
		Interior loop ( 69,315) CG; ( 70,314) GC:  -240
		Interior loop ( 70,314) GC; ( 71,313) GC:  -330
		Interior loop ( 71,313) GC; ( 72,312) UA:  -220
		Interior loop ( 72,312) UA; ( 78,306) GU:   290
		Interior loop ( 78,306) GU; ( 79,305) CG:  -250
		Interior loop ( 79,305) CG; ( 80,304) AU:  -210
		Interior loop ( 80,304) AU; ( 81,303) GC:  -210
		Interior loop ( 81,303) GC; ( 82,302) GU:  -150
		Interior loop ( 82,302) GU; ( 83,301) GU:   -50
		Interior loop ( 83,301) GU; ( 84,300) UA:  -140
		Interior loop ( 84,300) UA; ( 85,299) CG:  -240
		Interior loop ( 85,299) CG; ( 86,298) GU:  -140
		Interior loop ( 86,298) GU; ( 87,297) GC:  -210
		Interior loop ( 87,297) GC; ( 92,295) GC:   380
		Interior loop ( 92,295) GC; ( 93,294) GC:  -330
		Interior loop ( 93,294) GC; ( 94,293) AU:  -240
		Interior loop ( 94,293) AU; ( 96,291) AU:   190
		Interior loop ( 96,291) AU; ( 97,290) GU:   -60
		Interior loop ( 97,290) GU; ( 98,289) CG:  -250
		Interior loop ( 98,289) CG; ( 99,288) GC:  -240
		Interior loop ( 99,288) GC; (100,287) CG:  -340
		Interior loop (100,287) CG; (101,286) AU:  -210
		Interior loop (101,286) AU; (103,284) GC:   120
		Interior loop (103,284) GC; (104,283) AU:  -240
		Interior loop (104,283) AU; (105,282) GU:   -60
		Interior loop (105,282) GU; (106,281) GU:   -50
		Interior loop (106,281) GU; (107,280) GC:  -210
		Interior loop (107,280) GC; (108,279) AU:  -240
		Interior loop (108,279) AU; (109,278) GU:   -60
		Interior loop (109,278) GU; (110,277) CG:  -250
		Interior loop (113,130) CG; (114,129) CG:  -330
		Interior loop (114,129) CG; (115,128) AU:  -210
		Interior loop (115,128) AU; (116,127) GC:  -210
		Interior loop (116,127) GC; (117,126) GC:  -330
		Hairpin  loop (117,126) GC              :   330
		Interior loop (142,215) UA; (143,214) CG:  -240
		Interior loop (143,214) CG; (144,213) CG:  -330
		Interior loop (144,213) CG; (149,210) GU:   330
		Interior loop (149,210) GU; (150,209) GC:  -210
		Interior loop (150,209) GC; (151,208) GC:  -330
		Interior loop (151,208) GC; (152,207) UG:  -250
		Interior loop (152,207) UG; (153,206) UA:   -60
		Interior loop (153,206) UA; (154,205) UG:  -130
		Interior loop (154,205) UG; (155,204) CG:  -150
		Interior loop (155,204) CG; (156,203) GC:  -240
		Interior loop (156,203) GC; (157,202) CG:  -340
		Interior loop (157,202) CG; (158,201) CG:  -330
		Interior loop (158,201) CG; (160,200) CG:    50
		Interior loop (160,200) CG; (161,199) CG:  -330
		Interior loop (161,199) CG; (162,198) UG:  -210
		Interior loop (162,198) UG; (163,197) CG:  -150
		Interior loop (163,197) CG; (164,196) UA:  -210
		Interior loop (164,196) UA; (165,195) GC:  -210
		Interior loop (165,195) GC; (166,194) AU:  -240
		Interior loop (166,194) AU; (167,193) CG:  -220
		Interior loop (167,193) CG; (170,192) GC:   280
		Interior loop (170,192) GC; (171,191) AU:  -240
		Interior loop (171,191) AU; (172,190) GC:  -210
		Interior loop (172,190) GC; (173,189) CG:  -340
		Interior loop (173,189) CG; (174,188) GU:  -140
		Interior loop (174,188) GU; (175,187) UA:  -140
		Interior loop (175,187) UA; (176,186) CG:  -240
		Interior loop (176,186) CG; (177,185) GU:  -140
		Hairpin  loop (177,185) GU              :   570
		Interior loop (216,268) AU; (217,267) AU:   -90
		Interior loop (217,267) AU; (218,266) AU:   -90
		Interior loop (218,266) AU; (219,265) AU:   -90
		Interior loop (219,265) AU; (221,263) GC:   120
		Interior loop (221,263) GC; (222,262) CG:  -340
		Interior loop (222,262) CG; (223,261) CG:  -330
		Interior loop (223,261) CG; (224,260) AU:  -210
		Interior loop (224,260) AU; (225,259) GC:  -210
		Interior loop (225,259) GC; (226,258) CG:  -340
		Interior loop (226,258) CG; (227,257) AU:  -210
		Interior loop (227,257) AU; (228,256) AU:   -90
		Interior loop (228,256) AU; (232,253) GC:   330
		Interior loop (232,253) GC; (233,252) GC:  -330
		Interior loop (233,252) GC; (234,251) CG:  -340
		Interior loop (234,251) CG; (235,250) CG:  -330
		Hairpin  loop (235,250) CG              :   480
		Multi    loop (110,277) CG              :   290
		Interior loop (329,524) GC; (330,523) AU:  -240
		Interior loop (330,523) AU; (331,522) GC:  -210
		Interior loop (331,522) GC; (332,521) UA:  -220
		Interior loop (332,521) UA; (333,520) GC:  -210
		Interior loop (333,520) GC; (334,519) AU:  -240
		Interior loop (334,519) AU; (335,518) GC:  -210
		Interior loop (335,518) GC; (336,517) CG:  -340
		Interior loop (336,517) CG; (337,516) UA:  -210
		Interior loop (337,516) UA; (338,515) GU:  -100
		Interior loop (338,515) GU; (339,514) AU:  -130
		Interior loop (339,514) AU; (340,513) UG:  -140
		Interior loop (344,377) GU; (345,376) CG:  -250
		Interior loop (345,376) CG; (346,375) UA:  -210
		Interior loop (346,375) UA; (347,374) CG:  -240
		Interior loop (347,374) CG; (348,373) GC:  -240
		Interior loop (348,373) GC; (349,372) CG:  -340
		Interior loop (349,372) CG; (351,370) GC:    50
		Interior loop (351,370) GC; (352,369) CG:  -340
		Interior loop (352,369) CG; (354,368) GC:   140
		Interior loop (354,368) GC; (355,367) CG:  -340
		Interior loop (355,367) CG; (356,365) CG:    50
		Interior loop (356,365) CG; (357,364) GC:  -240
		Hairpin  loop (357,364) GC              :   390
		Interior loop (384,429) GC; (385,428) CG:  -340
		Interior loop (385,428) CG; (386,427) GC:  -240
		Interior loop (386,427) GC; (388,425) GC:    30
		Interior loop (388,425) GC; (389,424) GC:  -330
		Interior loop (389,424) GC; (390,423) AU:  -240
		Interior loop (390,423) AU; (392,419) GC:   300
		Interior loop (392,419) GC; (393,418) CG:  -340
		Interior loop (393,418) CG; (394,417) GC:  -240
		Interior loop (394,417) GC; (395,416) GC:  -330
		Interior loop (395,416) GC; (400,412) GC:   270
		Interior loop (400,412) GC; (401,411) CG:  -340
		Interior loop (401,411) CG; (402,410) GC:  -240
		Hairpin  loop (402,410) GC              :   490
		Interior loop (430,501) GC; (431,500) CG:  -340
		Interior loop (431,500) CG; (432,499) GC:  -240
		Interior loop (432,499) GC; (433,498) UA:  -220
		Interior loop (433,498) UA; (434,497) UA:   -90
		Interior loop (434,497) UA; (435,496) GC:  -210
		Interior loop (435,496) GC; (438,495) CG:   280
		Interior loop (438,495) CG; (439,494) GC:  -240
		Interior loop (439,494) GC; (441,493) UG:   130
		Interior loop (441,493) UG; (442,492) UA:   -60
		Interior loop (442,492) UA; (443,491) CG:  -240
		Interior loop (443,491) CG; (444,490) AU:  -210
		Interior loop (444,490) AU; (445,489) UG:  -140
		Interior loop (445,489) UG; (446,488) UA:   -60
		Interior loop (446,488) UA; (453,487) GC:   490
		Interior loop (453,487) GC; (454,486) CG:  -340
		Interior loop (454,486) CG; (455,485) UG:  -210
		Interior loop (455,485) UG; (457,483) GC:  -140
		Interior loop (457,483) GC; (458,482) CG:  -340
		Interior loop (458,482) CG; (463,477) CG:   130
		Interior loop (463,477) CG; (464,476) AU:  -210
		Interior loop (464,476) AU; (465,475) GC:  -210
		Hairpin  loop (465,475) GC              :   390
		Multi    loop (340,513) UG              :   220
		Interior loop (527,646) UG; (528,645) AU:  -100
		Interior loop (528,645) AU; (529,644) GC:  -210
		Interior loop (529,644) GC; (530,643) GC:  -330
		Interior loop (530,643) GC; (531,642) CG:  -340
		Interior loop (531,642) CG; (532,641) AU:  -210
		Interior loop (532,641) AU; (536,638) CG:   330
		Interior loop (536,638) CG; (537,637) AU:  -210
		Interior loop (537,637) AU; (538,636) GU:   -60
		Interior loop (538,636) GU; (539,635) GC:  -210
		Interior loop (539,635) GC; (540,634) CG:  -340
		Interior loop (540,634) CG; (541,633) UA:  -210
		Interior loop (541,633) UA; (542,632) UA:   -90
		Interior loop (542,632) UA; (549,624) UA:   420
		Interior loop (549,624) UA; (550,623) UG:  -130
		Interior loop (550,623) UG; (551,622) AU:  -100
		Interior loop (551,622) AU; (552,621) UA:  -110
		Interior loop (552,621) UA; (553,620) GC:  -210
		Interior loop (553,620) GC; (555,618) UA:   120
		Interior loop (555,618) UA; (556,617) UG:  -130
		Interior loop (556,617) UG; (559,614) GU:   280
		Interior loop (559,614) GU; (560,613) GC:  -210
		Interior loop (560,613) GC; (561,612) CG:  -340
		Interior loop (561,612) CG; (562,611) UA:  -210
		Interior loop (562,611) UA; (564,610) GC:   170
		Interior loop (564,610) GC; (565,609) UA:  -220
		Interior loop (565,609) UA; (569,605) UG:   240
		Interior loop (569,605) UG; (570,604) UA:   -60
		Interior loop (570,604) UA; (571,603) GC:  -210
		Interior loop (571,603) GC; (572,602) UA:  -220
		Interior loop (572,602) UA; (573,601) GC:  -210
		Interior loop (573,601) GC; (574,600) UA:  -220
		Interior loop (574,600) UA; (575,599) GC:  -210
		Interior loop (575,599) GC; (576,598) GU:  -150
		Interior loop (576,598) GU; (577,597) AU:  -130
		Interior loop (577,597) AU; (578,596) AU:   -90
		Interior loop (578,596) AU; (579,595) UA:  -110
		Interior loop (579,595) UA; (580,594) UA:   -90
		Interior loop (580,594) UA; (581,593) GC:  -210
		Interior loop (581,593) GC; (582,592) UA:  -220
		Hairpin  loop (582,592) UA              :   490
		Interior loop (648,2185) AU; (649,2184) GU:   -60
		Interior loop (649,2184) GU; (650,2183) GU:   -50
		Interior loop (650,2183) GU; (651,2182) UA:  -140
		Interior loop (651,2182) UA; (652,2181) CG:  -240
		Interior loop (652,2181) CG; (653,2180) GU:  -140
		Interior loop (653,2180) GU; (654,2179) AU:  -130
		Interior loop (654,2179) AU; (656,2178) UA:   270
		Interior loop (656,2178) UA; (657,2177) CG:  -240
		Interior loop (657,2177) CG; (658,2176) UA:  -210
		Interior loop (658,2176) UA; (659,2175) AU:  -130
		Interior loop (659,2175) AU; (660,2174) GU:   -60
		Interior loop (660,2174) GU; (661,2173) AU:  -130
		Interior loop (661,2173) AU; (662,2172) GC:  -210
		Interior loop (669,683) CG; (670,682) GC:  -240
		Interior loop (670,682) GC; (671,681) GU:  -150
		Interior loop (671,681) GU; (672,680) GC:  -210
		Interior loop (672,680) GC; (673,679) UG:  -250
		Hairpin  loop (673,679) UG              :   520
		Interior loop (691,1683) UA; (692,1682) GC:  -210
		Interior loop (692,1682) GC; (693,1681) GC:  -330
		Interior loop (693,1681) GC; (694,1680) CG:  -340
		Interior loop (694,1680) CG; (695,1679) CG:  -330
		Interior loop (695,1679) CG; (696,1678) GC:  -240
		Interior loop (696,1678) GC; (697,1677) UG:  -250
		Interior loop (697,1677) UG; (699,1675) GC:   120
		Interior loop (699,1675) GC; (700,1674) UA:  -220
		Interior loop (700,1674) UA; (702,1672) UA:   190
		Interior loop (702,1672) UA; (703,1671) UA:   -90
		Interior loop (703,1671) UA; (704,1670) AU:  -130
		Interior loop (704,1670) AU; (707,1667) AU:   150
		Interior loop (707,1667) AU; (708,1666) CG:  -220
		Interior loop (708,1666) CG; (710,1665) UA:   170
		Interior loop (710,1665) UA; (711,1664) CG:  -240
		Interior loop (711,1664) CG; (712,1663) GU:  -140
		Interior loop (712,1663) GU; (713,1662) UA:  -140
		Interior loop (713,1662) UA; (714,1661) GC:  -210
		Interior loop (718,727) GC; (719,726) GC:  -330
		Interior loop (719,726) GC; (720,725) GC:  -330
		Hairpin  loop (720,725) GC              :   450
		Interior loop (729,750) GC; (730,749) GC:  -330
		Interior loop (730,749) GC; (731,748) CG:  -340
		Interior loop (731,748) CG; (732,747) GC:  -240
		Interior loop (732,747) GC; (733,745) UA:   160
		Interior loop (733,745) UA; (734,744) UA:   -90
		Interior loop (734,744) UA; (735,743) AU:  -130
		Hairpin  loop (735,743) AU              :   580
		Interior loop (753,1654) GC; (754,1653) CG:  -340
		Interior loop (754,1653) CG; (755,1652) AU:  -210
		Interior loop (755,1652) AU; (756,1651) GC:  -210
		Interior loop (756,1651) GC; (757,1650) CG:  -340
		Interior loop (757,1650) CG; (758,1649) AU:  -210
		Interior loop (758,1649) AU; (759,1648) CG:  -220
		Interior loop (765,845) CG; (766,844) CG:  -330
		Interior loop (766,844) CG; (768,842) UA:    80
		Interior loop (768,842) UA; (769,841) UA:   -90
		Interior loop (769,841) UA; (770,840) CG:  -240
		Interior loop (770,840) CG; (771,839) GC:  -240
		Interior loop (771,839) GC; (772,838) CG:  -340
		Interior loop (772,838) CG; (773,837) CG:  -330
		Interior loop (773,837) CG; (774,836) AU:  -210
		Interior loop (774,836) AU; (775,830) GC:   450
		Interior loop (775,830) GC; (776,829) CG:  -340
		Interior loop (776,829) CG; (777,828) UA:  -210
		Interior loop (777,828) UA; (778,827) GC:  -210
		Interior loop (778,827) GC; (780,826) CG:    40
		Interior loop (780,826) CG; (781,825) GC:  -240
		Interior loop (781,825) GC; (782,824) UG:  -250
		Interior loop (782,824) UG; (783,823) AU:  -100
		Interior loop (783,823) AU; (784,822) AU:   -90
		Interior loop (784,822) AU; (785,821) UG:  -140
		Interior loop (785,821) UG; (789,814) GC:   490
		Interior loop (789,814) GC; (790,813) AU:  -240
		Interior loop (790,813) AU; (791,812) AU:   -90
		Interior loop (791,812) AU; (792,811) GC:  -210
		Interior loop (792,811) GC; (794,810) GC:    50
		Interior loop (794,810) GC; (795,809) GC:  -330
		Interior loop (795,809) GC; (796,808) CG:  -340
		Interior loop (796,808) CG; (798,804) CG:   230
		Interior loop (798,804) CG; (799,803) GC:  -240
		Hairpin  loop (799,803) GC              :   540
		Interior loop (846,856) CG; (847,855) GC:  -240
		Interior loop (847,855) GC; (848,854) CG:  -340
		Hairpin  loop (848,854) CG              :   490
		Interior loop (857,1646) GC; (858,1645) UG:  -250
		Interior loop (858,1645) UG; (859,1644) AU:  -100
		Interior loop (859,1644) AU; (860,1643) UA:  -110
		Interior loop (860,1643) UA; (862,1640) UA:   300
		Interior loop (862,1640) UA; (863,1639) UA:   -90
		Interior loop (863,1639) UA; (864,1638) CG:  -240
		Interior loop (864,1638) CG; (865,1637) UA:  -210
		Interior loop (865,1637) UA; (866,1636) CG:  -240
		Interior loop (866,1636) CG; (868,1635) UA:   170
		Interior loop (868,1635) UA; (869,1634) UA:   -90
		Interior loop (869,1634) UA; (870,1633) AU:  -130
		Interior loop (870,1633) AU; (871,1632) CG:  -220
		Interior loop (872,905) GC; (873,904) CG:  -340
		Interior loop (873,904) CG; (874,903) AU:  -210
		Interior loop (874,903) AU; (875,902) UG:  -140
		Interior loop (875,902) UG; (876,901) CG:  -150
		Interior loop (876,901) CG; (877,899) UA:   170
		Interior loop (877,899) UA; (878,898) GU:  -100
		Interior loop (878,898) GU; (879,897) UA:  -140
		Interior loop (879,897) UA; (880,896) GC:  -210
		Interior loop (880,896) GC; (881,895) CG:  -340
		Interior loop (881,895) CG; (882,894) GC:  -240
		Interior loop (882,894) GC; (883,893) GC:  -330
		Interior loop (883,893) GC; (884,892) UA:  -220
		Hairpin  loop (884,892) UA              :   570
		Interior loop (911,922) CG; (912,921) AU:  -210
		Interior loop (912,921) AU; (913,920) GC:  -210
		Hairpin  loop (913,920) GC              :   380
		Interior loop (927,1333) GU; (928,1332) AU:  -130
		Interior loop (928,1332) AU; (929,1331) UA:  -110
		Interior loop (929,1331) UA; (930,1330) GC:  -210
		Interior loop (930,1330) GC; (931,1329) CG:  -340
		Interior loop (931,1329) CG; (932,1328) CG:  -330
		Interior loop (932,1328) CG; (933,1327) GC:  -240
		Interior loop (933,1327) GC; (934,1326) CG:  -340
		Interior loop (934,1326) CG; (935,1325) AU:  -210
		Interior loop (935,1325) AU; (937,1323) AU:   150
		Interior loop (937,1323) AU; (938,1322) GU:   -60
		Interior loop (938,1322) GU; (951,1308) CG:   420
		Interior loop (951,1308) CG; (952,1307) GC:  -240
		Interior loop (952,1307) GC; (953,1306) AU:  -240
		Interior loop (953,1306) AU; (954,1305) CG:  -220
		Interior loop (954,1305) CG; (955,1304) AU:  -210
		Interior loop (955,1304) AU; (956,1303) CG:  -220
		Interior loop (959,970) GC; (960,969) CG:  -340
		Hairpin  loop (960,969) CG              :   440
		Interior loop (972,993) GC; (973,992) AU:  -240
		Interior loop (973,992) AU; (974,991) CG:  -220
		Interior loop (974,991) CG; (975,990) GU:  -140
		Interior loop (975,990) GU; (977,988) GC:   100
		Interior loop (977,988) GC; (978,987) CG:  -340
		Interior loop (978,987) CG; (979,986) CG:  -330
		Interior loop (979,986) CG; (980,985) CG:  -330
		Hairpin  loop (980,985) CG              :   420
		Interior loop (995,1040) GC; (996,1039) CG:  -340
		Interior loop (996,1039) CG; (997,1038) UA:  -210
		Interior loop (997,1038) UA; (998,1037) CG:  -240
		Interior loop (998,1037) CG; (999,1036) CG:  -330
		Interior loop (999,1036) CG; (1000,1035) CG:  -330
		Interior loop (1000,1035) CG; (1001,1034) GC:  -240
		Interior loop (1001,1034) GC; (1002,1033) GC:  -330
		Interior loop (1002,1033) GC; (1007,1029) CG:   270
		Interior loop (1007,1029) CG; (1008,1028) GC:  -240
		Interior loop (1008,1028) GC; (1010,1026) UA:   120
		Interior loop (1010,1026) UA; (1011,1025) UG:  -130
		Interior loop (1011,1025) UG; (1012,1024) AU:  -100
		Interior loop (1012,1024) AU; (1013,1023) CG:  -220
		Interior loop (1013,1023) CG; (1014,1022) AU:  -210
		Interior loop (1014,1022) AU; (1015,1021) GC:  -210
		Hairpin  loop (1015,1021) GC              :   440
		Interior loop (1042,1257) GC; (1043,1256) CG:  -340
		Interior loop (1043,1256) CG; (1044,1255) AU:  -210
		Interior loop (1044,1255) AU; (1045,1254) UA:  -110
		Interior loop (1045,1254) UA; (1047,1252) UA:   190
		Interior loop (1047,1252) UA; (1048,1251) GU:  -100
		Interior loop (1048,1251) GU; (1049,1250) UA:  -140
		Interior loop (1049,1250) UA; (1050,1249) CG:  -240
		Interior loop (1050,1249) CG; (1051,1248) AU:  -210
		Interior loop (1051,1248) AU; (1052,1247) GC:  -210
		Interior loop (1052,1247) GC; (1054,1246) GC:    50
		Interior loop (1054,1246) GC; (1055,1245) GC:  -330
		Interior loop (1055,1245) GC; (1056,1244) UA:  -220
		Interior loop (1056,1244) UA; (1057,1243) UA:   -90
		Interior loop (1057,1243) UA; (1064,1239) GC:   430
		Interior loop (1064,1239) GC; (1065,1238) UA:  -220
		Interior loop (1065,1238) UA; (1066,1237) CG:  -240
		Interior loop (1066,1237) CG; (1068,1236) UA:   170
		Interior loop (1068,1236) UA; (1069,1235) CG:  -240
		Interior loop (1071,1169) CG; (1072,1168) CG:  -330
		Interior loop (1072,1168) CG; (1073,1167) GC:  -240
		Interior loop (1073,1167) GC; (1074,1166) AU:  -240
		Interior loop (1074,1166) AU; (1075,1165) AU:   -90
		Interior loop (1075,1165) AU; (1076,1164) AU:   -90
		Interior loop (1076,1164) AU; (1078,1160) GC:   300
		Interior loop (1078,1160) GC; (1079,1159) CG:  -340
		Interior loop (1079,1159) CG; (1080,1157) GU:   240
		Interior loop (1080,1157) GU; (1081,1156) CG:  -250
		Interior loop (1081,1156) CG; (1084,1153) GC:   -30
		Interior loop (1084,1153) GC; (1085,1152) AU:  -240
		Interior loop (1085,1152) AU; (1086,1151) CG:  -220
		Interior loop (1086,1151) CG; (1087,1150) GC:  -240
		Interior loop (1087,1150) GC; (1089,1146) AU:   300
		Interior loop (1089,1146) AU; (1090,1145) AU:   -90
		Interior loop (1090,1145) AU; (1091,1144) GC:  -210
		Interior loop (1091,1144) GC; (1092,1143) GU:  -150
		Interior loop (1092,1143) GU; (1093,1142) GU:   -50
		Interior loop (1093,1142) GU; (1094,1140) CG:   130
		Interior loop (1094,1140) CG; (1095,1139) CG:  -330
		Interior loop (1095,1139) CG; (1096,1131) UA:   510
		Interior loop (1096,1131) UA; (1097,1130) CG:  -240
		Interior loop (1097,1130) CG; (1098,1129) GU:  -140
		Interior loop (1098,1129) GU; (1099,1128) UA:  -140
		Interior loop (1099,1128) UA; (1100,1127) GC:  -210
		Interior loop (1100,1127) GC; (1101,1126) AU:  -240
		Interior loop (1101,1126) AU; (1102,1125) UG:  -140
		Interior loop (1102,1125) UG; (1103,1124) AU:  -100
		Interior loop (1103,1124) AU; (1105,1120) GU:   370
		Interior loop (1105,1120) GU; (1106,1119) CG:  -250
		Interior loop (1106,1119) CG; (1107,1118) CG:  -330
		Interior loop (1107,1118) CG; (1108,1117) UA:  -210
		Interior loop (1108,1117) UA; (1109,1116) AU:  -130
		Interior loop (1109,1116) AU; (1110,1115) UA:  -110
		Hairpin  loop (1110,1115) UA              :   470
		Interior loop (1176,1234) GU; (1177,1233) UA:  -140
		Interior loop (1177,1233) UA; (1178,1232) GC:  -210
		Interior loop (1178,1232) GC; (1180,1230) GC:     0
		Interior loop (1180,1230) GC; (1181,1229) CG:  -340
		Interior loop (1181,1229) CG; (1182,1228) GC:  -240
		Interior loop (1182,1228) GC; (1183,1227) GC:  -330
		Interior loop (1183,1227) GC; (1184,1226) AU:  -240
		Interior loop (1184,1226) AU; (1185,1224) AU:   290
		Interior loop (1185,1224) AU; (1186,1223) CG:  -220
		Interior loop (1186,1223) CG; (1190,1221) UA:   300
		Interior loop (1190,1221) UA; (1191,1220) AU:  -130
		Interior loop (1191,1220) AU; (1192,1219) UA:  -110
		Interior loop (1192,1219) UA; (1193,1218) UA:   -90
		Interior loop (1193,1218) UA; (1194,1217) UA:   -90
		Interior loop (1194,1217) UA; (1195,1216) GC:  -210
		Interior loop (1195,1216) GC; (1198,1211) UA:   330
		Interior loop (1198,1211) UA; (1199,1210) AU:  -130
		Interior loop (1199,1210) AU; (1200,1209) UA:  -110
		Interior loop (1200,1209) UA; (1201,1208) UA:   -90
		Interior loop (1201,1208) UA; (1202,1207) UA:   -90
		Hairpin  loop (1202,1207) UA              :   470
		Multi    loop (1069,1235) CG              :   320
		Interior loop (1258,1272) UA; (1259,1271) UA:   -90
		Interior loop (1259,1271) UA; (1260,1270) CG:  -240
		Interior loop (1260,1270) CG; (1261,1269) AU:  -210
		Interior loop (1261,1269) AU; (1262,1268) AU:   -90
		Interior loop (1262,1268) AU; (1263,1267) UA:  -110
		Hairpin  loop (1263,1267) UA              :   590
		Interior loop (1276,1302) GC; (1277,1301) GC:  -330
		Interior loop (1277,1301) GC; (1278,1300) AU:  -240
		Interior loop (1278,1300) AU; (1279,1299) AU:   -90
		Interior loop (1279,1299) AU; (1280,1298) GU:   -60
		Interior loop (1280,1298) GU; (1282,1296) GC:    60
		Interior loop (1282,1296) GC; (1283,1295) UA:  -220
		Interior loop (1283,1295) UA; (1285,1294) UA:   290
		Interior loop (1285,1294) UA; (1286,1293) GC:  -210
		Interior loop (1286,1293) GC; (1287,1292) AU:  -240
		Hairpin  loop (1287,1292) AU              :   510
		Multi    loop (956,1303) CG              :  -500
		Interior loop (1335,1626) UA; (1336,1625) GC:  -210
		Interior loop (1336,1625) GC; (1337,1624) CG:  -340
		Interior loop (1337,1624) CG; (1338,1623) CG:  -330
		Interior loop (1338,1623) CG; (1340,1621) UA:    80
		Interior loop (1340,1621) UA; (1341,1620) CG:  -240
		Interior loop (1341,1620) CG; (1342,1619) CG:  -330
		Interior loop (1348,1377) UA; (1349,1376) UA:   -90
		Interior loop (1349,1376) UA; (1350,1375) GU:  -100
		Interior loop (1350,1375) GU; (1351,1374) CG:  -250
		Interior loop (1351,1374) CG; (1352,1371) UA:   330
		Interior loop (1352,1371) UA; (1353,1370) CG:  -240
		Interior loop (1353,1370) CG; (1354,1369) AU:  -210
		Interior loop (1354,1369) AU; (1355,1368) CG:  -220
		Interior loop (1355,1368) CG; (1356,1367) CG:  -330
		Interior loop (1356,1367) CG; (1358,1366) AU:   170
		Interior loop (1358,1366) AU; (1359,1365) GC:  -210
		Hairpin  loop (1359,1365) GC              :   440
		Interior loop (1378,1616) AU; (1379,1615) AU:   -90
		Interior loop (1379,1615) AU; (1380,1614) GC:  -210
		Interior loop (1380,1614) GC; (1381,1613) AU:  -240
		Interior loop (1381,1613) AU; (1382,1612) UA:  -110
		Interior loop (1382,1612) UA; (1383,1611) GC:  -210
		Interior loop (1383,1611) GC; (1384,1610) CG:  -340
		Interior loop (1384,1610) CG; (1385,1609) UA:  -210
		Interior loop (1388,1571) AU; (1389,1570) GC:  -210
		Interior loop (1389,1570) GC; (1390,1569) AU:  -240
		Interior loop (1390,1569) AU; (1393,1566) AU:   170
		Interior loop (1393,1566) AU; (1394,1565) GC:  -210
		Interior loop (1394,1565) GC; (1395,1564) UA:  -220
		Interior loop (1395,1564) UA; (1396,1560) UA:   420
		Interior loop (1396,1560) UA; (1397,1559) GC:  -210
		Interior loop (1397,1559) GC; (1398,1557) GC:    50
		Interior loop (1398,1557) GC; (1399,1556) GC:  -330
		Interior loop (1399,1556) GC; (1400,1555) UG:  -250
		Interior loop (1400,1555) UG; (1401,1554) GC:  -140
		Interior loop (1401,1554) GC; (1404,1551) CG:   -10
		Interior loop (1404,1551) CG; (1405,1550) GC:  -240
		Interior loop (1405,1550) GC; (1406,1549) AU:  -240
		Interior loop (1406,1549) AU; (1407,1548) GC:  -210
		Interior loop (1407,1548) GC; (1408,1547) UA:  -220
		Interior loop (1409,1446) GC; (1410,1445) GC:  -330
		Interior loop (1410,1445) GC; (1411,1444) GU:  -150
		Interior loop (1411,1444) GU; (1412,1443) UA:  -140
		Interior loop (1412,1443) UA; (1413,1442) UG:  -130
		Interior loop (1413,1442) UG; (1416,1439) AU:   240
		Interior loop (1416,1439) AU; (1417,1438) UG:  -140
		Interior loop (1417,1438) UG; (1418,1437) CG:  -150
		Interior loop (1418,1437) CG; (1419,1436) GC:  -240
		Interior loop (1419,1436) GC; (1422,1435) CG:   280
		Interior loop (1422,1435) CG; (1423,1434) UA:  -210
		Interior loop (1423,1434) UA; (1424,1433) GC:  -210
		Hairpin  loop (1424,1433) GC              :   300
		Interior loop (1447,1495) UA; (1448,1494) UA:   -90
		Interior loop (1448,1494) UA; (1449,1493) GU:  -100
		Interior loop (1449,1493) GU; (1450,1492) AU:  -130
		Interior loop (1450,1492) AU; (1451,1491) GU:   -60
		Interior loop (1451,1491) GU; (1452,1490) AU:  -130
		Interior loop (1452,1490) AU; (1453,1489) GC:  -210
		Interior loop (1453,1489) GC; (1454,1488) UA:  -220
		Interior loop (1454,1488) UA; (1456,1486) UG:   190
		Interior loop (1456,1486) UG; (1457,1485) UA:   -60
		Interior loop (1457,1485) UA; (1458,1484) CG:  -240
		Interior loop (1458,1484) CG; (1459,1483) GU:  -140
		Interior loop (1459,1483) GU; (1464,1476) GC:   440
		Interior loop (1464,1476) GC; (1465,1475) AU:  -240
		Interior loop (1465,1475) AU; (1466,1474) AU:   -90
		Interior loop (1466,1474) AU; (1467,1473) GU:   -60
		Interior loop (1467,1473) GU; (1468,1472) AU:  -130
		Hairpin  loop (1468,1472) AU              :   590
		Interior loop (1497,1545) GC; (1498,1544) UG:  -250
		Interior loop (1498,1544) UG; (1499,1543) UA:   -60
		Interior loop (1499,1543) UA; (1500,1542) CG:  -240
		Interior loop (1500,1542) CG; (1501,1541) UA:  -210
		Interior loop (1501,1541) UA; (1502,1539) GC:   170
		Interior loop (1502,1539) GC; (1503,1538) CG:  -340
		Interior loop (1503,1538) CG; (1504,1537) UG:  -210
		Interior loop (1504,1537) UG; (1508,1536) UG:   420
		Interior loop (1508,1536) UG; (1509,1535) GC:  -140
		Interior loop (1509,1535) GC; (1510,1534) GC:  -330
		Interior loop (1510,1534) GC; (1511,1533) CG:  -340
		Interior loop (1511,1533) CG; (1512,1532) GC:  -240
		Interior loop (1512,1532) GC; (1513,1530) CG:    40
		Interior loop (1513,1530) CG; (1514,1529) GU:  -140
		Interior loop (1514,1529) GU; (1515,1528) GU:   -50
		Interior loop (1515,1528) GU; (1516,1527) UA:  -140
		Interior loop (1516,1527) UA; (1517,1526) AU:  -130
		Interior loop (1517,1526) AU; (1518,1525) UG:  -140
		Hairpin  loop (1518,1525) UG              :   530
		Multi    loop (1408,1547) UA              :   190
		Interior loop (1577,1602) UA; (1578,1601) GC:  -210
		Interior loop (1578,1601) GC; (1579,1600) AU:  -240
		Interior loop (1579,1600) AU; (1580,1599) CG:  -220
		Interior loop (1580,1599) CG; (1581,1598) UA:  -210
		Interior loop (1581,1598) UA; (1583,1597) GC:   170
		Interior loop (1583,1597) GC; (1584,1596) GC:  -330
		Interior loop (1584,1596) GC; (1585,1595) UA:  -220
		Interior loop (1585,1595) UA; (1587,1594) GC:   170
		Interior loop (1587,1594) GC; (1588,1593) AU:  -240
		Hairpin  loop (1588,1593) AU              :   510
		Multi    loop (1385,1609) UA              :   570
		Multi    loop (1342,1619) CG              :   420
		Multi    loop (871,1632) CG              :   110
		Multi    loop (759,1648) CG              :    70
		Multi    loop (714,1661) GC              :   100
		Interior loop (1685,2039) CG; (1686,2038) UA:  -210
		Interior loop (1686,2038) UA; (1687,2037) UA:   -90
		Interior loop (1687,2037) UA; (1688,2036) AU:  -130
		Interior loop (1688,2036) AU; (1689,2035) CG:  -220
		Interior loop (1689,2035) CG; (1690,2034) UG:  -210
		Interior loop (1690,2034) UG; (1691,2032) UA:   320
		Interior loop (1691,2032) UA; (1692,2031) CG:  -240
		Interior loop (1692,2031) CG; (1693,2030) UA:  -210
		Interior loop (1696,1771) CG; (1697,1770) AU:  -210
		Interior loop (1697,1770) AU; (1698,1769) AU:   -90
		Interior loop (1698,1769) AU; (1699,1768) CG:  -220
		Interior loop (1699,1768) CG; (1700,1767) GC:  -240
		Interior loop (1700,1767) GC; (1701,1766) AU:  -240
		Interior loop (1701,1766) AU; (1702,1765) UA:  -110
		Interior loop (1702,1765) UA; (1703,1764) CG:  -240
		Interior loop (1703,1764) CG; (1704,1763) GU:  -140
		Interior loop (1704,1763) GU; (1706,1762) AU:   250
		Interior loop (1706,1762) AU; (1707,1761) GC:  -210
		Interior loop (1707,1761) GC; (1708,1760) GC:  -330
		Interior loop (1711,1733) CG; (1712,1732) GU:  -140
		Interior loop (1712,1732) GU; (1713,1731) AU:  -130
		Interior loop (1713,1731) AU; (1714,1730) AU:   -90
		Interior loop (1714,1730) AU; (1715,1729) GU:   -60
		Interior loop (1715,1729) GU; (1716,1728) GU:   -50
		Interior loop (1716,1728) GU; (1717,1727) AU:  -130
		Interior loop (1717,1727) AU; (1718,1726) GC:  -210
		Interior loop (1718,1726) GC; (1719,1725) CG:  -340
		Hairpin  loop (1719,1725) CG              :   430
		Interior loop (1738,1753) AU; (1739,1752) CG:  -220
		Interior loop (1739,1752) CG; (1740,1751) AU:  -210
		Interior loop (1740,1751) AU; (1741,1750) UA:  -110
		Interior loop (1741,1750) UA; (1742,1749) GC:  -210
		Interior loop (1742,1749) GC; (1743,1748) GU:  -150
		Hairpin  loop (1743,1748) GU              :   460
		Multi    loop (1708,1760) GC              :   410
		Interior loop (1772,2029) GC; (1773,2028) GC:  -330
		Interior loop (1773,2028) GC; (1776,2025) CG:  -190
		Interior loop (1776,2025) CG; (1777,2024) CG:  -330
		Interior loop (1777,2024) CG; (1778,2023) GU:  -140
		Interior loop (1778,2023) GU; (1779,2022) GC:  -210
		Interior loop (1779,2022) GC; (1781,2020) GC:    90
		Interior loop (1781,2020) GC; (1782,2019) CG:  -340
		Interior loop (1782,2019) CG; (1783,2018) UA:  -210
		Interior loop (1783,2018) UA; (1784,2017) GC:  -210
		Interior loop (1784,2017) GC; (1785,2015) AU:   140
		Interior loop (1785,2015) AU; (1786,2014) AU:   -90
		Interior loop (1786,2014) AU; (1787,2013) UA:  -110
		Interior loop (1787,2013) UA; (1788,2012) GC:  -210
		Interior loop (1788,2012) GC; (1789,2011) AU:  -240
		Interior loop (1789,2011) AU; (1791,2009) GU:   190
		Interior loop (1791,2009) GU; (1792,2008) CG:  -250
		Interior loop (1792,2008) CG; (1793,2007) CG:  -330
		Interior loop (1793,2007) CG; (1808,2006) GC:   530
		Interior loop (1808,2006) GC; (1809,2005) CG:  -340
		Interior loop (1809,2005) CG; (1810,2004) GC:  -240
		Interior loop (1810,2004) GC; (1812,2002) GC:   -10
		Interior loop (1812,2002) GC; (1813,2001) AU:  -240
		Interior loop (1813,2001) AU; (1814,2000) CG:  -220
		Interior loop (1814,2000) CG; (1816,1999) CG:    50
		Interior loop (1816,1999) CG; (1817,1998) CG:  -330
		Interior loop (1817,1998) CG; (1818,1997) AU:  -210
		Interior loop (1818,1997) AU; (1819,1996) CG:  -220
		Interior loop (1819,1996) CG; (1820,1995) GC:  -240
		Interior loop (1826,1851) UA; (1827,1850) GC:  -210
		Interior loop (1827,1850) GC; (1828,1849) UG:  -250
		Interior loop (1828,1849) UG; (1830,1848) GC:   240
		Interior loop (1830,1848) GC; (1831,1847) CG:  -340
		Interior loop (1831,1847) CG; (1832,1846) AU:  -210
		Interior loop (1832,1846) AU; (1833,1845) AU:   -90
		Interior loop (1833,1845) AU; (1834,1844) UG:  -140
		Interior loop (1834,1844) UG; (1835,1843) GC:  -140
		Hairpin  loop (1835,1843) GC              :   350
		Interior loop (1854,1902) CG; (1855,1901) UA:  -210
		Interior loop (1855,1901) UA; (1856,1900) AU:  -130
		Interior loop (1856,1900) AU; (1857,1899) UA:  -110
		Interior loop (1857,1899) UA; (1858,1898) UA:   -90
		Interior loop (1858,1898) UA; (1859,1897) AU:  -130
		Interior loop (1859,1897) AU; (1860,1896) AU:   -90
		Interior loop (1860,1896) AU; (1861,1889) CG:   490
		Interior loop (1861,1889) CG; (1862,1888) UG:  -210
		Interior loop (1862,1888) UG; (1863,1887) GC:  -140
		Interior loop (1863,1887) GC; (1864,1886) GC:  -330
		Interior loop (1864,1886) GC; (1866,1885) GC:    50
		Interior loop (1866,1885) GC; (1867,1884) AU:  -240
		Interior loop (1867,1884) AU; (1868,1883) AU:   -90
		Interior loop (1868,1883) AU; (1869,1881) CG:   160
		Interior loop (1869,1881) CG; (1870,1880) UA:  -210
		Interior loop (1870,1880) UA; (1871,1879) AU:  -130
		Hairpin  loop (1871,1879) AU              :   580
		Interior loop (1904,1956) CG; (1905,1955) UG:  -210
		Interior loop (1905,1955) UG; (1906,1954) GC:  -140
		Interior loop (1906,1954) GC; (1907,1953) GC:  -330
		Interior loop (1907,1953) GC; (1908,1952) AU:  -240
		Interior loop (1908,1952) AU; (1910,1951) GU:   320
		Interior loop (1910,1951) GU; (1911,1950) GC:  -210
		Interior loop (1911,1950) GC; (1913,1949) GC:    50
		Interior loop (1913,1949) GC; (1914,1948) GC:  -330
		Interior loop (1914,1948) GC; (1915,1947) CG:  -340
		Interior loop (1915,1947) CG; (1917,1945) GC:  -140
		Interior loop (1917,1945) GC; (1918,1944) AU:  -240
		Interior loop (1918,1944) AU; (1923,1943) GC:   410
		Interior loop (1923,1943) GC; (1924,1942) UG:  -250
		Interior loop (1924,1942) UG; (1926,1941) GC:   240
		Interior loop (1926,1941) GC; (1927,1940) CG:  -340
		Interior loop (1927,1940) CG; (1928,1939) AU:  -210
		Interior loop (1928,1939) AU; (1929,1938) GC:  -210
		Interior loop (1929,1938) GC; (1930,1937) GU:  -150
		Interior loop (1930,1937) GU; (1931,1936) AU:  -130
		Hairpin  loop (1931,1936) AU              :   540
		Interior loop (1957,1990) CG; (1958,1989) UG:  -210
		Interior loop (1958,1989) UG; (1959,1988) GC:  -140
		Interior loop (1959,1988) GC; (1960,1987) GC:  -330
		Interior loop (1960,1987) GC; (1961,1986) CG:  -340
		Interior loop (1961,1986) CG; (1962,1985) UA:  -210
		Interior loop (1962,1985) UA; (1963,1981) GC:   370
		Interior loop (1963,1981) GC; (1964,1980) GU:  -150
		Interior loop (1964,1980) GU; (1965,1979) UA:  -140
		Interior loop (1965,1979) UA; (1966,1978) UA:   -90
		Interior loop (1966,1978) UA; (1967,1977) UA:   -90
		Interior loop (1967,1977) UA; (1968,1976) AU:  -130
		Interior loop (1968,1976) AU; (1969,1975) UA:  -110
		Interior loop (1969,1975) UA; (1970,1974) UG:  -130
		Hairpin  loop (1970,1974) UG              :   590
		Multi    loop (1820,1995) GC              :   -30
		Multi    loop (1693,2030) UA              :   330
		Interior loop (2043,2072) UG; (2044,2071) CG:  -150
		Interior loop (2044,2071) CG; (2045,2070) CG:  -330
		Interior loop (2045,2070) CG; (2046,2069) CG:  -330
		Interior loop (2046,2069) CG; (2047,2068) GC:  -240
		Interior loop (2047,2068) GC; (2048,2067) UA:  -220
		Interior loop (2048,2067) UA; (2052,2063) GC:   190
		Interior loop (2052,2063) GC; (2053,2062) UA:  -220
		Interior loop (2053,2062) UA; (2054,2061) AU:  -130
		Interior loop (2054,2061) AU; (2055,2060) GC:  -210
		Hairpin  loop (2055,2060) GC              :   400
		Interior loop (2073,2128) AU; (2074,2127) GC:  -210
		Interior loop (2074,2127) GC; (2075,2126) UA:  -220
		Interior loop (2075,2126) UA; (2077,2124) AU:   190
		Interior loop (2077,2124) AU; (2078,2123) GC:  -210
		Interior loop (2078,2123) GC; (2079,2122) GC:  -330
		Interior loop (2079,2122) GC; (2080,2121) CG:  -340
		Interior loop (2080,2121) CG; (2081,2120) AU:  -210
		Interior loop (2081,2120) AU; (2083,2118) CG:   120
		Interior loop (2083,2118) CG; (2084,2117) UA:  -210
		Interior loop (2084,2117) UA; (2085,2116) AU:  -130
		Interior loop (2085,2116) AU; (2086,2115) UA:  -110
		Interior loop (2086,2115) UA; (2094,2109) CG:   360
		Interior loop (2094,2109) CG; (2095,2108) GC:  -240
		Interior loop (2095,2108) GC; (2096,2107) AU:  -240
		Hairpin  loop (2096,2107) AU              :   620
		Interior loop (2129,2164) GC; (2130,2163) AU:  -240
		Interior loop (2130,2163) AU; (2132,2161) UA:   190
		Interior loop (2132,2161) UA; (2133,2160) AU:  -130
		Interior loop (2133,2160) AU; (2134,2159) AU:   -90
		Interior loop (2134,2159) AU; (2135,2158) GU:   -60
		Interior loop (2135,2158) GU; (2136,2157) CG:  -250
		Interior loop (2136,2157) CG; (2138,2156) UA:   170
		Interior loop (2138,2156) UA; (2139,2155) UA:   -90
		Interior loop (2139,2155) UA; (2140,2154) GC:  -210
		Interior loop (2140,2154) GC; (2141,2153) GC:  -330
		Interior loop (2141,2153) GC; (2142,2152) UA:  -220
		Hairpin  loop (2142,2152) UA              :   590
		Multi    loop (662,2172) GC              :  -440
		Interior loop (2201,2234) AU; (2202,2233) UA:  -110
		Interior loop (2202,2233) UA; (2203,2232) UG:  -130
		Interior loop (2203,2232) UG; (2205,2230) AU:   190
		Interior loop (2205,2230) AU; (2206,2229) AU:   -90
		Interior loop (2206,2229) AU; (2207,2228) AU:   -90
		Interior loop (2207,2228) AU; (2208,2227) AU:   -90
		Interior loop (2208,2227) AU; (2209,2226) GC:  -210
		Interior loop (2209,2226) GC; (2210,2225) GC:  -330
		Interior loop (2210,2225) GC; (2211,2224) AU:  -240
		Interior loop (2211,2224) AU; (2212,2223) UA:  -110
		Interior loop (2212,2223) UA; (2213,2222) CG:  -240
		Interior loop (2213,2222) CG; (2214,2221) UA:  -210
		Hairpin  loop (2214,2221) UA              :   490
		Multi    loop (  5,2236) UA              :  -220
		Interior loop (2260,2678) GC; (2261,2677) UA:  -220
		Interior loop (2261,2677) UA; (2262,2676) GU:  -100
		Interior loop (2262,2676) GU; (2264,2674) GC:   120
		Interior loop (2264,2674) GC; (2265,2673) UA:  -220
		Interior loop (2265,2673) UA; (2266,2672) UG:  -130
		Interior loop (2266,2672) UG; (2267,2670) UA:   320
		Interior loop (2267,2670) UA; (2268,2669) UA:   -90
		Interior loop (2268,2669) UA; (2269,2668) CG:  -240
		Interior loop (2269,2668) CG; (2270,2667) GC:  -240
		Interior loop (2270,2667) GC; (2271,2666) UG:  -250
		Interior loop (2271,2666) UG; (2272,2665) UA:   -60
		Interior loop (2272,2665) UA; (2273,2664) CG:  -240
		Interior loop (2273,2664) CG; (2274,2663) CG:  -330
		Interior loop (2274,2663) CG; (2275,2662) AU:  -210
		Interior loop (2278,2350) GU; (2279,2349) AU:  -130
		Interior loop (2279,2349) AU; (2280,2348) GC:  -210
		Interior loop (2280,2348) GC; (2281,2347) CG:  -340
		Interior loop (2281,2347) CG; (2282,2346) GU:  -140
		Interior loop (2282,2346) GU; (2284,2344) CG:   120
		Interior loop (2284,2344) CG; (2285,2343) AU:  -210
		Interior loop (2285,2343) AU; (2286,2342) GC:  -210
		Interior loop (2286,2342) GC; (2287,2341) AU:  -240
		Interior loop (2287,2341) AU; (2291,2335) CG:   420
		Interior loop (2291,2335) CG; (2292,2334) GC:  -240
		Interior loop (2292,2334) GC; (2293,2333) UG:  -250
		Interior loop (2293,2333) UG; (2294,2332) AU:  -100
		Interior loop (2294,2332) AU; (2295,2331) GC:  -210
		Interior loop (2295,2331) GC; (2296,2330) AU:  -240
		Interior loop (2296,2330) AU; (2297,2329) AU:   -90
		Interior loop (2297,2329) AU; (2298,2328) AU:   -90
		Interior loop (2298,2328) AU; (2299,2327) AU:   -90
		Interior loop (2299,2327) AU; (2300,2326) GU:   -60
		Interior loop (2300,2326) GU; (2301,2325) AU:  -130
		Interior loop (2301,2325) AU; (2306,2324) AU:   460
		Interior loop (2306,2324) AU; (2307,2323) GC:  -210
		Interior loop (2307,2323) GC; (2308,2322) GC:  -330
		Interior loop (2308,2322) GC; (2309,2321) AU:  -240
		Interior loop (2309,2321) AU; (2310,2320) UA:  -110
		Interior loop (2310,2320) UA; (2311,2319) CG:  -240
		Interior loop (2311,2319) CG; (2312,2318) UA:  -210
		Interior loop (2312,2318) UA; (2313,2317) UG:  -130
		Hairpin  loop (2313,2317) UG              :   590
		Interior loop (2351,2391) GC; (2352,2390) CG:  -340
		Interior loop (2352,2390) CG; (2353,2389) AU:  -210
		Interior loop (2353,2389) AU; (2354,2388) AU:   -90
		Interior loop (2354,2388) AU; (2355,2387) AU:   -90
		Interior loop (2355,2387) AU; (2356,2386) CG:  -220
		Interior loop (2356,2386) CG; (2357,2385) AU:  -210
		Interior loop (2357,2385) AU; (2358,2384) AU:   -90
		Interior loop (2358,2384) AU; (2363,2383) AU:   460
		Interior loop (2363,2383) AU; (2364,2382) CG:  -220
		Interior loop (2364,2382) CG; (2365,2381) CG:  -330
		Interior loop (2365,2381) CG; (2366,2380) AU:  -210
		Interior loop (2366,2380) AU; (2367,2379) CG:  -220
		Interior loop (2367,2379) CG; (2368,2378) CG:  -330
		Interior loop (2368,2378) CG; (2369,2377) GC:  -240
		Interior loop (2369,2377) GC; (2370,2376) CG:  -340
		Interior loop (2370,2376) CG; (2371,2375) UA:  -210
		Hairpin  loop (2371,2375) UA              :   590
		Interior loop (2392,2420) CG; (2393,2419) GC:  -240
		Interior loop (2393,2419) GC; (2394,2418) GC:  -330
		Interior loop (2394,2418) GC; (2395,2417) AU:  -240
		Interior loop (2395,2417) AU; (2398,2414) AU:   170
		Interior loop (2398,2414) AU; (2399,2413) AU:   -90
		Interior loop (2399,2413) AU; (2400,2412) GC:  -210
		Interior loop (2400,2412) GC; (2401,2411) AU:  -240
		Interior loop (2401,2411) AU; (2402,2410) GC:  -210
		Hairpin  loop (2402,2410) GC              :   490
		Interior loop (2423,2545) GC; (2424,2544) GC:  -330
		Interior loop (2424,2544) GC; (2425,2543) UA:  -220
		Interior loop (2425,2543) UA; (2426,2542) AU:  -130
		Interior loop (2426,2542) AU; (2427,2541) AU:   -90
		Interior loop (2427,2541) AU; (2428,2540) CG:  -220
		Interior loop (2428,2540) CG; (2430,2538) GC:    40
		Interior loop (2430,2538) GC; (2431,2537) GC:  -330
		Interior loop (2431,2537) GC; (2433,2535) UA:   120
		Interior loop (2433,2535) UA; (2434,2534) UA:   -90
		Interior loop (2434,2534) UA; (2436,2533) AU:   250
		Interior loop (2436,2533) AU; (2437,2532) GC:  -210
		Interior loop (2437,2532) GC; (2438,2531) CG:  -340
		Interior loop (2438,2531) CG; (2439,2530) AU:  -210
		Interior loop (2439,2530) AU; (2440,2529) GC:  -210
		Interior loop (2440,2529) GC; (2441,2528) AU:  -240
		Interior loop (2441,2528) AU; (2442,2527) GC:  -210
		Interior loop (2442,2527) GC; (2443,2526) CG:  -340
		Interior loop (2443,2526) CG; (2444,2525) GC:  -240
		Interior loop (2444,2525) GC; (2446,2524) AU:   140
		Interior loop (2446,2524) AU; (2447,2523) GC:  -210
		Interior loop (2460,2501) GC; (2461,2500) UA:  -220
		Interior loop (2461,2500) UA; (2462,2499) UA:   -90
		Interior loop (2462,2499) UA; (2463,2498) CG:  -240
		Interior loop (2463,2498) CG; (2464,2497) UA:  -210
		Interior loop (2464,2497) UA; (2465,2496) UA:   -90
		Interior loop (2465,2496) UA; (2468,2493) AU:   190
		Interior loop (2468,2493) AU; (2469,2492) GC:  -210
		Interior loop (2469,2492) GC; (2470,2491) UA:  -220
		Interior loop (2470,2491) UA; (2471,2490) GC:  -210
		Interior loop (2471,2490) GC; (2474,2486) GC:   260
		Interior loop (2474,2486) GC; (2475,2485) CG:  -340
		Interior loop (2475,2485) CG; (2476,2484) CG:  -330
		Hairpin  loop (2476,2484) CG              :   370
		Interior loop (2504,2519) UA; (2505,2518) GC:  -210
		Interior loop (2505,2518) GC; (2506,2517) UA:  -220
		Interior loop (2506,2517) UA; (2507,2516) AU:  -130
		Interior loop (2507,2516) AU; (2508,2515) GC:  -210
		Hairpin  loop (2508,2515) GC              :   470
		Multi    loop (2447,2523) GC              :   360
		Interior loop (2549,2626) GC; (2550,2625) GU:  -150
		Interior loop (2550,2625) GU; (2551,2624) CG:  -250
		Interior loop (2551,2624) CG; (2552,2623) UG:  -210
		Interior loop (2552,2623) UG; (2553,2622) GC:  -140
		Interior loop (2553,2622) GC; (2554,2621) CG:  -340
		Interior loop (2554,2621) CG; (2555,2620) UA:  -210
		Interior loop (2555,2620) UA; (2556,2619) GC:  -210
		Interior loop (2556,2619) GC; (2557,2618) CG:  -340
		Interior loop (2562,2572) GC; (2563,2571) GU:  -150
		Interior loop (2563,2571) GU; (2564,2570) CG:  -250
		Hairpin  loop (2564,2570) CG              :   340
		Interior loop (2575,2617) GC; (2576,2616) UG:  -250
		Interior loop (2576,2616) UG; (2577,2615) CG:  -150
		Interior loop (2577,2615) CG; (2578,2614) UA:  -210
		Interior loop (2578,2614) UA; (2579,2613) UA:   -90
		Interior loop (2579,2613) UA; (2580,2612) AU:  -130
		Interior loop (2580,2612) AU; (2581,2610) CG:   160
		Interior loop (2581,2610) CG; (2582,2609) CG:  -330
		Interior loop (2582,2609) CG; (2583,2608) GC:  -240
		Interior loop (2583,2608) GC; (2584,2607) GC:  -330
		Interior loop (2584,2607) GC; (2585,2601) GU:   450
		Interior loop (2585,2601) GU; (2586,2600) UA:  -140
		Interior loop (2586,2600) UA; (2587,2599) UG:  -130
		Interior loop (2587,2599) UG; (2588,2598) GC:  -140
		Hairpin  loop (2588,2598) GC              :   390
		Multi    loop (2557,2618) CG              :   250
		Interior loop (2627,2657) GC; (2628,2656) GC:  -330
		Interior loop (2628,2656) GC; (2629,2655) GC:  -330
		Interior loop (2629,2655) GC; (2630,2654) CG:  -340
		Interior loop (2630,2654) CG; (2631,2653) UA:  -210
		Interior loop (2631,2653) UA; (2632,2652) GC:  -210
		Interior loop (2632,2652) GC; (2634,2646) AU:   510
		Interior loop (2634,2646) AU; (2635,2645) CG:  -220
		Interior loop (2635,2645) CG; (2636,2644) GC:  -240
		Interior loop (2636,2644) GC; (2637,2643) GU:  -150
		Hairpin  loop (2637,2643) GU              :   520
		Multi    loop (2275,2662) AU              :  -430
		-936.50`, t)
}

func compareMFEOutputToViennaRNA(sequence, structure string, temperature float64, energyParamsSet energy_params.EnergyParamsSet, dangleModel DanglingEndsModel, expected_output string, t *testing.T) {
	output := captureOutput(func() {
		mfe, secondaryStructure, err := MinimumFreeEnergy(sequence, structure, temperature, energyParamsSet, dangleModel)
		if err != nil {
			panic(err)
		}
		logEnergyContributions(*secondaryStructure, sequence)
		// mfe = math.Floor(mfe*100) / 100
		log.Printf("%.2f", mfe)
	})

	if stripWhiteSpace(output) != stripWhiteSpace(expected_output) {
		t.Errorf("\n\nFailed to calcualte mfe for '%v'. \nExpected: \n'%v'\n\nGot: \n'%v'\n\n",
			sequence, expected_output, output)
	}
}

// Captures any logged output from running func `f` and returns it as a string
func captureOutput(f func()) string {
	var buf bytes.Buffer
	log.SetFlags(0)
	log.SetOutput(&buf)
	f()
	log.SetOutput(os.Stderr)
	return buf.String()
}

// removes all white space in a string
func stripWhiteSpace(s string) string {
	return strings.Join(strings.Fields(s), "")
}

/**
* Logs energy contributions in the same format as ViennaRNA does. This is used
* to compare output to ViennaRNA in test cases.
* Note: 1 is added to the indexes of the closing and enclosed base pairs as
* everyting in ViennaRNA is 1-indexed, but output from the `MinimumFreeEnergy` func
* is 0-indexed.
 */
func logEnergyContributions(secondaryStructure SecondaryStructure, sequence string) {
	log.Printf("External loop       : %v\n", int(secondaryStructure.ExteriorLoopsEnergy))
	doLogEnergyContributions(secondaryStructure.Structures, sequence)
}

func doLogEnergyContributions(structures []interface{}, sequence string) {
	for _, structure := range structures {
		switch structure := structure.(type) {
		case *MultiLoop:
			// process Stem
			logStemEnergyContributions(structure.Stem, sequence)
			// process multiloop substructures
			doLogEnergyContributions(structure.Substructures, sequence)
			// process multiloop energy
			var i, j int = structure.Stem.EnclosedFivePrimeIdx, structure.Stem.EnclosedThreePrimeIdx
			log.Printf("Multi   loop (%v,%v) %v%v              : %v\n",
				i+1, j+1,
				string(sequence[i]), string(sequence[j]),
				int(structure.Energy))
		case *Hairpin:
			// process Stem
			logStemEnergyContributions(structure.Stem, sequence)
			// process hairpin energy
			var i, j int = structure.Stem.EnclosedFivePrimeIdx, structure.Stem.EnclosedThreePrimeIdx
			log.Printf("Hairpin loop  (%v,%v) %v%v              : %v\n",
				i+1, j+1,
				string(sequence[i]), string(sequence[j]),
				int(structure.Energy))
		}
	}
}

func logStemEnergyContributions(stem Stem, sequence string) {
	for _, stemStructure := range stem.Structures {
		logStemStructureEnergyContribution(stemStructure, sequence)
	}
}

func logStemStructureEnergyContribution(stemStructure StemStructure, sequence string) {
	var i, j int = stemStructure.ClosingFivePrimeIdx, stemStructure.ClosingThreePrimeIdx
	var k, l int = stemStructure.EnclosedFivePrimeIdx, stemStructure.EnclosedThreePrimeIdx
	log.Printf("Interior loop (%v,%v) %v%v; (%v,%v) %v%v: %v\n",
		i+1, j+1,
		string(sequence[i]), string(sequence[j]),
		k+1, l+1,
		string(sequence[k]), string(sequence[l]),
		int(stemStructure.Energy))
}
