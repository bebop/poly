package poly

import (
	"container/heap"
	"errors"
	"fmt"
	"math"
	"sort"
	"strings"
	"unicode"
)

// vivek: Magic numbers
// Should add source and small description of what the number is about
var (
	external_unpaired float64 = -0.009729
	multi_unpaired    float64 = -0.198330
	external_paired   float64 = -0.000967
	multi_paired      float64 = -0.925388
	multi_base        float64 = -1.199055

	bulge_length              = [31]float64{0.0, -2.399548472, -3.2940667837, -4.2029218746, -5.0441693501, -5.4807172844, -6.0506360645, -5.8503526421, -5.0964765063, -5.7009810518, -6.4210758616, -6.9347480537, -7.2962207216, -7.5576661608, -7.7170588501, -7.80330553291, -7.83437644287, -7.84534866319, -7.81533646036, -7.76774522247, -7.81070694312, -7.82862593974, -7.90663145496, -7.97762471926, -8.03530424822, -8.08164219503, -8.11723639959, -8.14398574353, -8.16217532325, -8.17269833057, -8.17785195742}
	internal_length           = [31]float64{0.0, 0.0, -0.429061443, -0.7822725931, -1.1786523466, -1.4897722641, -1.7449668113, -1.79645798028, -1.83964800435, -1.83766251486, -2.01381382847, -2.27778244917, -2.62384380687, -2.91650411476, -2.95274661783, -3.07274199393, -3.11628971319, -3.19838264454, -3.20549587058, -3.18194762206, -3.15127788635, -3.21746029729, -3.34906953559, -3.48986908699, -3.55587200561, -3.63366405305, -3.6845060657, -3.72590482171, -3.72262823831, -3.71670365547, -3.70982791746}
	internal_explicit         = [21]float64{0.0, 0.0, 0.0, 0.0, 0.0, -0.1754591076, 0.03083787104, -0.171565435, -0.2294680983, 0.0, -0.1304072693, -0.07730329553, 0.2782767264, 0.0, 0.0, -0.02898949617, 0.3112350694, 0.0, 0.0, 0.0, -0.3226348245}
	internal_symmetric_length = [16]float64{0.0, -0.5467082599, -0.9321784246, -1.1910250647, -1.4251087392, -1.2800509627, -1.9363442142, -2.2384530511, -2.26877580377, -2.62057020957, -2.83648346017, -2.95931050557, -3.11453136507, -3.1999425725, -3.24586367049, -3.26818601285}
	internal_asymmetry        = [29]float64{0.0, -2.105646719, -2.6576607621, -3.2347315291, -3.8483983138, -4.1541139979, -4.269619198, -4.4801804211, -4.7947547341, -5.1096509022, -5.19983279712, -5.41983547652, -5.56048380082, -5.77672492672, -5.94927807022, -6.10516925682, -6.20925512312, -6.2789319654, -6.31999174034, -6.3356979835, -6.32187797711, -6.28055809148, -6.24461623198, -6.21639436916, -6.20002851042, -6.17452794867, -6.14104762074, -6.10132837662, -6.10387349055}
	hairpin_length            = [31]float64{-5.993180158, -9.10128592, -8.6843882853, -6.4789692193, -4.5522195273, -5.1395440602, -5.222301238, -4.6439122536, -5.3660005908, -5.5385880532, -5.8410970399, -5.8707286338, -6.7976282286, -6.82920576838, -6.93145297848, -6.74131224388, -6.83412134214, -6.66507650134, -6.74680216605, -7.09139606915, -7.20054636315, -7.49089873245, -7.83027009915, -8.02180651085, -8.07199860464, -8.11074481388, -8.06323010636, -7.9957868871, -7.89856812984, -7.73125495654, -7.49826123164}
	terminal_mismatch         = [625]float64{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.184546064, -0.1181844187, -0.4461469607, -0.6175254495, 0.0, 0.004788458708, 0.08319395146, -0.2249479995, -0.3981327204, 0.0, 0.5191110288, -0.3524119307, -0.4056429433, -0.7733932162, 0.0, -0.01574403519, 0.268570042, -0.0934388741, 0.3373711531, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.08386423535, -0.2520716816, -0.6711841881, -0.3816350028, 0.0, 0.1117852189, -0.1704393624, -0.2179987732, -0.459267635, 0.0, 0.8520640313, -0.9332488517, -0.3289551692, -0.7778822056, 0.0, -0.2422339958, -0.03780509247, -0.4322334143, -0.2419976114, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.1703136025, -0.09154056357, -0.2522413002, -0.8520314799, 0.0, 0.04763224188, -0.2428654283, -0.2079275061, -0.1874270053, 0.0, 0.6540033983, -0.7823988605, 0.1995898255, -0.4432169392, 0.0, -0.1736921762, 0.288494362, -0.01638238057, 0.6757988971, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.4871607613, 0.1105031953, 0.363373916, -0.6193199348, 0.0, 0.3451056056, 0.0314944976, -0.3799172956, -0.03222973182, 0.0, 0.4948638637, -0.2821952552, -0.2702227211, -0.06658395291, 0.0, -0.4306154451, -0.09497863465, -0.3130794485, -0.2283242981, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0115363879, -0.3923408221, 0.05661063599, -0.1251485388, 0.0, -0.06545074758, -0.3167200568, 0.002258383981, -0.422217724, 0.0, 0.5458416646, -0.2085887954, -0.1971766062, -0.4722410132, 0.0, -0.1779642496, 0.1643454344, -0.5005617032, 0.1333867679, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1218741278, 0.1990260141, 0.04681893928, 0.3256264491, 0.0, 0.1186812326, -0.1851065102, -0.04311512683, -0.6150608139, 0.0, 0.754933218, -0.3150708483, 0.1569582926, -0.514970007, 0.0, -0.2926246029, 0.1373068149, -0.05422333363, 0.03086776921, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
	helix_closing             = [25]float64{0.0, 0.0, 0.0, -0.9770893163, 0.0, 0.0, 0.0, -0.4574650937, 0.0, 0.0, 0.0, -0.8265995623, 0.0, -1.051678928, 0.0, -0.9246140521, 0.0, -0.3698708172, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
	base_pair                 = [25]float64{0.0, 0.0, 0.0, 0.59791199, 0.0, 0.0, 0.0, 1.544290641, 0.0, 0.0, 0.0, 1.544290641, 0.0, -0.01304754992, 0.0, 0.59791199, 0.0, -0.01304754992, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
	dangle_left               = [125]float64{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.1251037681, 0.0441606708, -0.02541879082, 0.00785098466, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.07224381372, 0.05279281874, 0.1009554299, -0.1515059013, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.1829535099, 0.03393000394, 0.1335339061, -0.1604274506, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.06517511341, -0.04250882422, 0.02875971806, -0.04359727428, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.03373847659, -0.005070324324, -0.1186861149, -0.01162357727, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.08047139148, 0.001608000669, 0.1016272216, -0.09200842832, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
	dangle_right              = [125]float64{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.03232578201, -0.09096819493, -0.0740750973, -0.01621157379, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2133964379, -0.06234810991, -0.07008531041, -0.2141912285, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01581957549, 0.005644320058, -0.00943297687, -0.2597793095, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.04480271781, -0.07321213002, 0.01270494867, -0.05717033985, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.1631918513, 0.06769304994, -0.08789074414, -0.05525570007, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.04105458185, -0.008136642572, -0.03808592022, -0.08629373429, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
	helix_stacking            = [625]float64{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1482005248, 0.0, 0.0, 0.0, 0.4343497127, 0.0, 0.0, 0.0, 0.7079642577, 0.0, -0.1010777582, 0.0, 0.243256656, 0.0, 0.1623654243, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4878707793, 0.0, 0.0, 0.0, 0.8481320247, 0.0, 0.0, 0.0, 0.4784248478, 0.0, -0.1811268205, 0.0, 0.7079642577, 0.0, 0.4849351028, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5551785831, 0.0, 0.0, 0.0, 0.5008324248, 0.0, 0.0, 0.0, 0.8481320247, 0.0, 0.2165962476, 0.0, 0.4343497127, 0.0, 0.4864603589, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.04665365028, 0.0, 0.0, 0.0, 0.4864603589, 0.0, 0.0, 0.0, 0.4849351028, 0.0, 0.1833447295, 0.0, 0.1623654243, 0.0, -0.2858970755, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3897593783, 0.0, 0.0, 0.0, 0.5551785831, 0.0, 0.0, 0.0, 0.4878707793, 0.0, -0.1157333764, 0.0, 0.1482005248, 0.0, -0.04665365028, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.1157333764, 0.0, 0.0, 0.0, 0.2165962476, 0.0, 0.0, 0.0, -0.1811268205, 0.0, 0.120296538, 0.0, -0.1010777582, 0.0, 0.1833447295, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
	bulge_0x1_nucleotides     = [5]float64{-0.1216861662, -0.07111241127, 0.008947026647, -0.002685763742, 0.0}
	internal_1x1_nucleotides  = [25]float64{0.2944404686, 0.08641360967, -0.3664197228, -0.2053107048, 0.0, 0.08641360967, -0.1582543624, 0.4175273724, 0.1368762582, 0.0, -0.3664197228, 0.4175273724, -0.1193514754, -0.4188101413, 0.0, -0.2053107048, 0.1368762582, -0.4188101413, 0.147140653, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}

	_allowed_pairs  [NOTON][NOTON]bool
	_helix_stacking [NOTON][NOTON][NOTON][NOTON]bool
	cache_single    [SINGLE_MAX_LEN + 1][SINGLE_MAX_LEN + 1]float64
	scores          []Pair
	bestC           []*State
	bestH           [](map[int]*State)
	bestP           [](map[int]*State)
	bestM           [](map[int]*State)
	bestM2          [](map[int]*State)
	bestMulti       [](map[int]*State)
	sorted_bestM    [][]Pair

	beam_size int = 100 // vivek: recommended beam size based on paper
)

const (
	NOTON                 int = 5 // NUM_OF_TYPE_OF_NUCS
	NOTOND                int = 25
	NOTONT                int = 125
	SINGLE_MAX_LEN        int = 30 // NOTE: *must* <= sizeof(char), otherwise modify State::TraceInfo accordingly
	SINGLE_MIN_LEN        int = 0
	MIN_CUBE_PRUNING_SIZE int = 20

	HAIRPIN_MAX_LEN   int = 30
	EXPLICIT_MAX_LEN  int = 4
	SYMMETRIC_MAX_LEN int = 15
	ASYMMETRY_MAX_LEN int = 28

	INTERNAL_MAX_LEN int     = SINGLE_MAX_LEN
	VALUE_MIN        float64 = -math.MaxFloat64
)

// vivek: Func `initalize` in the original source
func init() {
	_allowed_pairs[GET_ACGU_NUM('A')][GET_ACGU_NUM('U')] = true
	_allowed_pairs[GET_ACGU_NUM('U')][GET_ACGU_NUM('A')] = true
	_allowed_pairs[GET_ACGU_NUM('C')][GET_ACGU_NUM('G')] = true
	_allowed_pairs[GET_ACGU_NUM('G')][GET_ACGU_NUM('C')] = true
	_allowed_pairs[GET_ACGU_NUM('G')][GET_ACGU_NUM('U')] = true
	_allowed_pairs[GET_ACGU_NUM('U')][GET_ACGU_NUM('G')] = true

	SET_HELIX_STACKING('A', 'U', 'A', 'U', true)
	SET_HELIX_STACKING('A', 'U', 'C', 'G', true)
	SET_HELIX_STACKING('A', 'U', 'G', 'C', true)
	SET_HELIX_STACKING('A', 'U', 'G', 'U', true)
	SET_HELIX_STACKING('A', 'U', 'U', 'A', true)
	SET_HELIX_STACKING('A', 'U', 'U', 'G', true)
	SET_HELIX_STACKING('C', 'G', 'A', 'U', true)
	SET_HELIX_STACKING('C', 'G', 'C', 'G', true)
	SET_HELIX_STACKING('C', 'G', 'G', 'C', true)
	SET_HELIX_STACKING('C', 'G', 'G', 'U', true)
	SET_HELIX_STACKING('C', 'G', 'U', 'G', true)
	SET_HELIX_STACKING('G', 'C', 'A', 'U', true)
	SET_HELIX_STACKING('G', 'C', 'C', 'G', true)
	SET_HELIX_STACKING('G', 'C', 'G', 'U', true)
	SET_HELIX_STACKING('G', 'C', 'U', 'G', true)
	SET_HELIX_STACKING('G', 'U', 'A', 'U', true)
	SET_HELIX_STACKING('G', 'U', 'G', 'U', true)
	SET_HELIX_STACKING('G', 'U', 'U', 'G', true)
	SET_HELIX_STACKING('U', 'A', 'A', 'U', true)
	SET_HELIX_STACKING('U', 'A', 'G', 'U', true)
	SET_HELIX_STACKING('U', 'G', 'G', 'U', true)

	InitCachesingle()
}

// vivek: Func `initalize_cachesingle` in the original source
func InitCachesingle() {
	for i := 0; i < SINGLE_MAX_LEN+1; i++ {
		for j := 0; j < SINGLE_MAX_LEN+1; j++ {
			cache_single[i][j] = 0
		}
	}

	for l1 := SINGLE_MIN_LEN; l1 <= SINGLE_MAX_LEN; l1++ {
		for l2 := SINGLE_MIN_LEN; l2 <= SINGLE_MAX_LEN; l2++ {
			if l1 == 0 && l2 == 0 {
				continue
			} else if l1 == 0 {
				cache_single[l1][l2] += bulge_length[l2]
			} else if l2 == 0 {
				cache_single[l1][l2] += bulge_length[l1]
			} else {
				// internal
				cache_single[l1][l2] += internal_length[Min(l1+l2, INTERNAL_MAX_LEN)]

				// internal explicit
				if l1 <= EXPLICIT_MAX_LEN && l2 <= EXPLICIT_MAX_LEN {
					if l1 <= l2 {
						cache_single[l1][l2] += internal_explicit[l1*EXPLICIT_MAX_LEN+l2]
					} else {
						cache_single[l1][l2] += internal_explicit[l2*EXPLICIT_MAX_LEN+l1]
					}
				}

				if l1 == l2 {
					// internal symmetry
					cache_single[l1][l2] += internal_symmetric_length[Min(l1, SYMMETRIC_MAX_LEN)]
				} else {
					// internal asymmetry
					var diff int = l1 - l2
					if diff < 0 {
						diff = -diff
					}
					cache_single[l1][l2] += internal_asymmetry[Min(diff, ASYMMETRY_MAX_LEN)]
				}
			}

		}
	}
}

func SET_HELIX_STACKING(x rune, y rune, z rune, w rune, val bool) {
	_helix_stacking[GET_ACGU_NUM(x)][GET_ACGU_NUM(y)][GET_ACGU_NUM(z)][GET_ACGU_NUM(w)] = val
}

func GET_ACGU_NUM(base_pair rune) int {
	switch base_pair {
	case 'A':
		return 0
	case 'C':
		return 1
	case 'G':
		return 2
	case 'U':
		return 3
	default:
		return 4
	}
}

func Min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func Max(a, b int) int {
	if a > b {
		return a
	}
	return b
}

// vivek: What is manner?
const (
	MANNER_NONE                  = iota // 0: empty
	MANNER_H                            // 1: hairpin candidate
	MANNER_HAIRPIN                      // 2: hairpin
	MANNER_SINGLE                       // 3: single
	MANNER_HELIX                        // 4: helix
	MANNER_MULTI                        // 5: multi = ..M2. [30 restriction on the left and jump on the right]
	MANNER_MULTI_eq_MULTI_plus_U        // 6: multi = multi + U
	MANNER_P_eq_MULTI                   // 7: P = (multi)
	MANNER_M2_eq_M_plus_P               // 8: M2 = M + P
	MANNER_M_eq_M2                      // 9: M = M2
	MANNER_M_eq_M_plus_U                // 10: M = M + U
	MANNER_M_eq_P                       // 11: M = P
	/* MANNER_C_eq_U, */
	/* MANNER_C_eq_P, */
	MANNER_C_eq_C_plus_U // 12: C = C + U
	MANNER_C_eq_C_plus_P // 13: C = C + P
)

type State struct {
	score  float64
	manner int
	trace  struct {
		split    int
		paddings struct {
			l1 rune
			l2 int
		}
	}
}

func newState() *State {
	return &State{score: VALUE_MIN, manner: MANNER_NONE}
}

func update_if_better(s *State, newscore float64, manner int) {
	if s.score < newscore {
		s.Set(newscore, manner)
	}
}

// vivek: update_if_better is an overloaded function in the original source
// I couldn't think of a better way to translate it to Go so I naively made
// three functions. Same for the function set.
func update_if_better2(s *State, newscore float64, manner int, l1 rune, l2 int) {
	if s.score < newscore || s.manner == MANNER_NONE {
		s.Set2(newscore, manner, l1, l2)
	}
}

func update_if_better3(s *State, newscore float64, manner int, split int) {
	if s.score < newscore || s.manner == MANNER_NONE {
		// ++ nos_set_update;
		s.Set3(newscore, manner, split)
	}
}

func (s *State) Set(score float64, manner int) {
	s.score = score
	s.manner = manner
}

func (s *State) Set2(score float64, manner int, l1 rune, l2 int) {
	s.score = score
	s.manner = manner
	s.trace.paddings.l1 = l1
	s.trace.paddings.l2 = l2
}

func (s *State) Set3(score float64, manner int, split int) {
	s.score = score
	s.manner = manner
	s.trace.split = split
}

func LinearFold(sequence string) (string, float64) {
	// convert to uppercase
	sequence = strings.ToUpper(sequence)

	// convert T to U
	sequence = strings.Replace(sequence, "T", "U", -1)

	// lhuang: moved inside loop, fixing an obscure but crucial bug in initialization
	// BeamCKYParser parser(beam_size);
	// BeamCKYParser::DecoderResult result = parser.parse(seq, NULL);
	return Parse(sequence)

	// #ifdef lv
	//         double printscore = (result.score / -100.0);
	// #else
	// var printscore float64 = result.score
	// #endif
	// fmt.Printf("%s (%.2f)\n", result.structure, printscore)

}

type Pair struct {
	first, second interface{}
}

// An PairHeap is a max-heap of pairs.
type PairHeap []Pair

// Methods required by heap.Interface
func (h PairHeap) Len() int           { return len(h) }
func (h PairHeap) Less(i, j int) bool { return h[i].first.(float64) < h[j].first.(float64) }
func (h PairHeap) Swap(i, j int)      { h[i], h[j] = h[j], h[i] }

func (h *PairHeap) Push(x interface{}) {
	*h = append(*h, x.(Pair))
}

func (h *PairHeap) Pop() interface{} {
	old := *h
	n := len(old)
	x := old[n-1]
	*h = old[0 : n-1]
	return x
}

// vivek: returns folded sequence and score
func Parse(sequence string) (string, float64) {

	// number of states
	var nos_H, nos_P, nos_M2, nos_M, nos_C, nos_Multi uint64 = 0, 0, 0, 0, 0, 0

	seq_length := len(sequence)

	nucs := make([]int, seq_length)

	bestC = make([]*State, seq_length)
	for i := 0; i < seq_length; i++ {
		bestC[i] = newState()
	}
	bestH = make([](map[int]*State), seq_length)
	bestP = make([](map[int]*State), seq_length)
	bestM = make([](map[int]*State), seq_length)
	bestM2 = make([](map[int]*State), seq_length)
	bestMulti = make([](map[int]*State), seq_length)
	for i := 0; i < seq_length; i++ {
		bestH[i] = make(map[int]*State)
		bestP[i] = make(map[int]*State)
		bestM[i] = make(map[int]*State)
		bestM2[i] = make(map[int]*State)
		bestMulti[i] = make(map[int]*State)
	}

	sorted_bestM = make([][]Pair, seq_length)

	// vector to store the scores at each beam temporarily for beam pruning
	scores = make([]Pair, seq_length)

	// Convert from ACGU to integers
	for i := 0; i < seq_length; i++ {
		nucs[i] = GET_ACGU_NUM((rune)(sequence[i]))
	}

	// vivek: What is next_pair?
	next_pair := make([][]int, NOTON)

	// Iterate through ACGU
	for nuci := 0; nuci < NOTON; nuci++ {
		next_pair[nuci] = make([]int, seq_length)
		for j := 0; j < seq_length; j++ {
			next_pair[nuci][j] = -1
		}

		next := -1
		for j := seq_length - 1; j >= 0; j-- {
			next_pair[nuci][j] = next
			if _allowed_pairs[nuci][nucs[j]] {
				next = j
			}
		}
	}

	if seq_length > 0 {
		bestC[0].Set(score_external_unpaired(0, 0), MANNER_C_eq_C_plus_U)
	}

	if seq_length > 1 {
		bestC[1].Set(score_external_unpaired(0, 1), MANNER_C_eq_C_plus_U)
	}

	nos_C++

	// from left to right
	for j := 0; j < seq_length; j++ {
		nucj := nucs[j]
		nucj1 := 0
		if (j + 1) < seq_length {
			nucj1 = nucs[j+1]
		} else {
			nucj1 = -1
		}

		var beamstepH *map[int]*State = &bestH[j]
		var beamstepMulti *map[int]*State = &bestMulti[j]
		var beamstepP *map[int]*State = &bestP[j]
		var beamstepM2 *map[int]*State = &bestM2[j]
		var beamstepM *map[int]*State = &bestM[j]
		var beamstepC *State = bestC[j]

		// beam of H
		{
			if beam_size > 0 && len(*beamstepH) > beam_size {
				BeamPrune(beamstepH)
			}

			{
				// for nucj put H(j, j_next) into H[j_next]
				jnext := next_pair[nucj][j]

				/*
					vivek: The original repo has a no_sharp_turn flag that's enabled by
					default. It wasn't explained in the paper nor in the source code.
					I emailed Liang Huang about this and this was the explanation he gave
					me:
					"—sharpturn means it allows hairpins with less than 3 unpaired
					nucleotides, e.g.,

					((()))
					(((.)))
					(((..)))

					These three structures are not allowed by default, unless —sharpturn
					is turned on."
				*/
				// if (no_sharp_turn) {
				for jnext-j < 4 && jnext != -1 {
					jnext = next_pair[nucj][jnext]
				}
				// }

				if jnext != -1 {
					nucjnext := nucs[jnext]
					var nucjnext_1 int
					if (jnext - 1) > -1 {
						nucjnext_1 = nucs[jnext-1]
					} else {
						nucjnext_1 = -1
					}
					var newscore float64
					newscore = score_hairpin(j, jnext, nucj, nucj1, nucjnext_1, nucjnext)

					// this candidate must be the best one at [j, jnext]
					// so no need to check the score
					if bestH[jnext][j] == nil {
						bestH[jnext][j] = newState()
					}

					update_if_better(bestH[jnext][j], newscore, MANNER_H)
					nos_H++
				}
			}

			{
				// for every state h in H[j]
				//   1. extend h(i, j) to h(i, jnext)
				//   2. generate p(i, j)
				for i, state := range *beamstepH {
					nuci := nucs[i]
					jnext := next_pair[nuci][j]

					// 2. generate p(i, j)
					// lisiz, change the order because of the constriants
					{
						if (*beamstepP)[i] == nil {
							(*beamstepP)[i] = newState()
						}
						update_if_better((*beamstepP)[i], state.score, MANNER_HAIRPIN)
						nos_P++
					}

					if jnext != -1 {
						var nuci1 int
						if (i + 1) < seq_length {
							nuci1 = nucs[i+1]
						} else {
							nuci1 = -1
						}

						nucjnext := nucs[jnext]

						var nucjnext_1 int
						if (jnext - 1) > -1 {
							nucjnext_1 = nucs[jnext-1]
						} else {
							nucjnext_1 = -1
						}

						// 1. extend h(i, j) to h(i, jnext)
						var newscore float64

						newscore = score_hairpin(i, jnext, nuci, nuci1, nucjnext_1, nucjnext)
						// this candidate must be the best one at [i, jnext]
						// so no need to check the score

						if bestH[jnext][i] == nil {
							bestH[jnext][i] = newState()
						}
						update_if_better(bestH[jnext][i], newscore, MANNER_H)
						nos_H++
					}
				}
			}
		}

		if j == 0 {
			continue
		}

		// beam of Multi
		{
			if beam_size > 0 && len(*beamstepMulti) > beam_size {
				BeamPrune(beamstepMulti)
			}

			// for every state in Multi[j]
			//   1. extend (i, j) to (i, jnext)
			//   2. generate P (i, j)
			for i, state := range *beamstepMulti {

				nuci := nucs[i]
				nuci1 := nucs[i+1]
				jnext := next_pair[nuci][j]

				// 2. generate P (i, j)
				// lisiz, change the order because of the constraits
				{
					var newscore float64
					newscore = state.score + score_multi(i, j, nuci, nuci1, nucs[j-1], nucj, seq_length)
					if (*beamstepP)[i] == nil {
						(*beamstepP)[i] = newState()
					}
					update_if_better((*beamstepP)[i], newscore, MANNER_P_eq_MULTI)
					nos_P++
				}

				// 1. extend (i, j) to (i, jnext)
				{
					new_l1 := state.trace.paddings.l1
					new_l2 := state.trace.paddings.l2 + jnext - j
					// if (jnext != -1 && new_l1 + new_l2 <= SINGLE_MAX_LEN) {
					if jnext != -1 {
						// 1. extend (i, j) to (i, jnext)
						newscore := state.score + score_multi_unpaired(j, jnext-1)
						// this candidate must be the best one at [i, jnext]
						// so no need to check the score
						if bestMulti[jnext][i] == nil {
							bestMulti[jnext][i] = newState()
						}
						update_if_better2(bestMulti[jnext][i], newscore, MANNER_MULTI_eq_MULTI_plus_U,
							new_l1,
							new_l2)
						nos_Multi++
					}
				}
			}
		}

		// beam of P
		{
			if beam_size > 0 && len(*beamstepP) > beam_size {
				BeamPrune(beamstepP)
			}

			// for every state in P[j]
			//   1. generate new helix/bulge
			//   2. M = P
			//   3. M2 = M + P
			//   4. C = C + P
			var use_cube_pruning bool = beam_size > MIN_CUBE_PRUNING_SIZE && len(*beamstepP) > MIN_CUBE_PRUNING_SIZE

			for i, state := range *beamstepP {
				nuci := nucs[i]
				var nuci_1 int
				if i-1 > -1 {
					nuci_1 = nucs[i-1]
				} else {
					nuci_1 = -1
				}

				// 2. M = P
				if i > 0 && j < seq_length-1 {
					newscore := score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length) + state.score
					if (*beamstepM)[i] == nil {
						(*beamstepM)[i] = newState()
					}
					update_if_better((*beamstepM)[i], newscore, MANNER_M_eq_P)
					nos_M++
				}

				// 3. M2 = M + P
				if !use_cube_pruning {
					k := i - 1
					if k > 0 && len(bestM[k]) != 0 {
						M1_score := score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length) + state.score
						// candidate list
						// bestM2_iter := (*beamstepM2)[i]
						for newi, state := range bestM[k] {
							// eq. to first convert P to M1, then M2/M = M + M1
							newscore := M1_score + state.score
							if (*beamstepM2)[newi] == nil {
								(*beamstepM2)[newi] = newState()
							}
							update_if_better3((*beamstepM2)[newi], newscore, MANNER_M2_eq_M_plus_P, k)
							//update_if_better(bestM[j][newi], newscore, MANNER_M_eq_M_plus_P, k);
							nos_M2++
							//++nos_M;
						}
					}
				}

				// 4. C = C + P
				{
					k := i - 1
					if k >= 0 {
						var prefix_C State = *bestC[k]
						if prefix_C.manner != MANNER_NONE {
							nuck := nuci_1
							nuck1 := nuci

							newscore := score_external_paired(k+1, j, nuck, nuck1,
								nucj, nucj1, seq_length) + prefix_C.score + state.score

							if beamstepC == nil {
								beamstepC = newState()
							}
							update_if_better3(beamstepC, newscore, MANNER_C_eq_C_plus_P, k)
							nos_C++
						}
					} else {
						newscore := score_external_paired(0, j, -1, nucs[0],
							nucj, nucj1, seq_length) + state.score
						if beamstepC == nil {
							beamstepC = newState()
						}
						update_if_better3(beamstepC, newscore, MANNER_C_eq_C_plus_P, -1)
						nos_C++
					}
				}
				//printf(" C = C + P at %d\n", j); fflush(stdout);

				// 1. generate new helix / single_branch
				// new state is of shape p..i..j..q
				if i > 0 && j < seq_length-1 {
					var precomputed float64 = score_junction_B(j, i, nucj, nucj1, nuci_1, nuci)
					for p := i - 1; p >= Max(i-SINGLE_MAX_LEN, 0); p-- {
						nucp := nucs[p]
						nucp1 := nucs[p+1]
						q := next_pair[nucp][j]

						for q != -1 && ((i-p)+(q-j)-2 <= SINGLE_MAX_LEN) {
							nucq := nucs[q]
							nucq_1 := nucs[q-1]

							if p == i-1 && q == j+1 {
								// helix
								var newscore float64 = score_helix(nucp, nucp1, nucq_1, nucq) + state.score
								if bestP[q][p] == nil {
									bestP[q][p] = newState()
								}
								update_if_better(bestP[q][p], newscore, MANNER_HELIX)
								nos_P++
							} else {
								// single branch

								var newscore float64 = score_junction_B(p, q, nucp, nucp1, nucq_1, nucq) +
									precomputed +
									score_single_without_junctionB(p, q, i, j,
										nuci_1, nuci, nucj, nucj1) +
									state.score

								if bestP[q][p] == nil {
									bestP[q][p] = newState()
								}
								update_if_better2(bestP[q][p], newscore, MANNER_SINGLE,
									rune(i-p), q-j)
								nos_P++
							}
							q = next_pair[nucp][q]
						}
					}
				}
				//printf(" helix / single at %d\n", j); fflush(stdout);
			}

			if use_cube_pruning {
				// 3. M2 = M + P with cube pruning
				var valid_Ps []int
				var M1_scores []float64

				for i, state := range *beamstepP {
					nuci := nucs[i]
					var nuci_1 int
					if i-1 > -1 {
						nuci_1 = nucs[i-1]
					} else {
						nuci_1 = -1
					}

					k := i - 1

					// group candidate Ps
					if k > 0 && len(bestM[k]) != 0 {
						if len(bestM[k]) != len(sorted_bestM[k]) {
							panic("bestM size != sorted_bestM size")
						}

						M1_score := score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length) + state.score

						// bestM2_iter = beamstepM2[i]

						valid_Ps = append(valid_Ps, i)
						M1_scores = append(M1_scores, M1_score)

					}
				}

				// build max heap
				// heap is of form (heuristic score, (index of i in valid_Ps, index of M in bestM[i-1]))
				// vector<pair<value_type, pair<int, int>>> heap;
				// var heap []Pair
				h := &PairHeap{}
				heap.Init(h)
				for p := 0; p < len(valid_Ps); p++ {
					i := valid_Ps[p]
					k := i - 1
					heap.Push(h, Pair{M1_scores[p] + sorted_bestM[k][0].first.(float64), Pair{p, 0}})
					// heap = append(heap, Pair{M1_scores[p] + sorted_bestM[k][0].first.(float64), Pair{p, 0}})
					// HEAP STUFF
					// ±±±±±±±±push_heap(heap.begin(), heap.end())
				}

				// start cube pruning
				// stop after beam size M2 states being filled
				var filled int = 0
				// exit when filled >= beam and current score < prev score
				var prev_score float64 = VALUE_MIN
				var current_score float64 = VALUE_MIN
				for (filled < beam_size || current_score == prev_score) && len(*h) != 0 {

					// HEAP STUFF
					top := (*h)[0]
					// auto & top = heap.front()
					prev_score = current_score
					current_score = top.first.(float64)
					index_P := top.second.(Pair).first.(int)
					index_M := top.second.(Pair).second.(int)
					i := valid_Ps[top.second.(Pair).first.(int)]
					k := i - 1
					newi := sorted_bestM[k][index_M].second.(int)
					var newscore float64 = M1_scores[index_P] + bestM[k][newi].score
					// HEAP STUFF
					// pop the greatest element off the heap
					// pop_heap(heap.begin(), heap.end())
					// heap.pop_back()
					heap.Pop(h)

					if (*beamstepM2)[newi].manner == MANNER_NONE {
						filled++
						if (*beamstepM2)[newi] == nil {
							(*beamstepM2)[newi] = newState()
						}
						update_if_better3((*beamstepM2)[newi], newscore, MANNER_M2_eq_M_plus_P, k)
						nos_M2++
					} else {
						if !((*beamstepM2)[newi].score > newscore-1e-8) {
							panic("beamstepM2[newi].score <= newscore - 1e-8")
						}
					}

					index_M++
					for index_M < len(sorted_bestM[k]) {
						// candidate_score is a heuristic score
						var candidate_score float64 = M1_scores[index_P] + sorted_bestM[k][index_M].first.(float64)
						var candidate_newi int = sorted_bestM[k][index_M].second.(int)

						var last State
						for _, v := range *beamstepM2 {
							last = *v
						}

						if *(*beamstepM2)[candidate_newi] == last {
							// HEAP STUFF
							heap.Push(h, Pair{candidate_score, Pair{index_P, index_M}})
							break
						} else {
							// based on the property of cube pruning, the new score must be worse
							// than the state already inserted
							// so we keep iterate through the candidate list to find the next
							// candidate
							index_M++
							if !((*beamstepM2)[candidate_newi].score >
								M1_scores[index_P]+bestM[k][candidate_newi].score-1e-8) {
								panic("beamstepM2[candidate_newi].score <= M1_scores[index_P] + bestM[k][candidate_newi].score - 1e-8")
							}
						}
					}
				}
			}
		}
		//printf("P at %d\n", j); fflush(stdout);

		// beam of M2
		{
			if beam_size > 0 && len(*beamstepM2) > beam_size {
				BeamPrune(beamstepM2)
			}

			// for every state in M2[j]
			//   1. multi-loop  (by extending M2 on the left)
			//   2. M = M2
			for i, state := range *beamstepM2 {
				// 2. M = M2
				{
					if (*beamstepM)[i] == nil {
						(*beamstepM)[i] = newState()
					}
					update_if_better((*beamstepM)[i], state.score, MANNER_M_eq_M2)
					nos_M++
				}

				// 1. multi-loop
				{
					for p := i - 1; p >= Max(i-SINGLE_MAX_LEN, 0); p-- {
						nucp := nucs[p]
						q := next_pair[nucp][j]

						if q != -1 && ((i - p - 1) <= SINGLE_MAX_LEN) {
							// the current shape is p..i M2 j ..q
							var newscore float64 = score_multi_unpaired(p+1, i-1) +
								score_multi_unpaired(j+1, q-1) + state.score

							if bestMulti[q][p] == nil {
								bestMulti[q][p] = newState()
							}
							update_if_better2(bestMulti[q][p], newscore, MANNER_MULTI,
								rune(i-p),
								q-j)
							nos_Multi++
							//q = next_pair[nucp][q];
						}
					}
				}
			}
		}
		//printf("M2 at %d\n", j); fflush(stdout);

		// beam of M
		{
			var threshold float64 = VALUE_MIN
			if beam_size > 0 && len(*beamstepM) > beam_size {
				threshold = BeamPrune(beamstepM)
			}

			sortM(threshold, beamstepM, sorted_bestM[j])

			// for every state in M[j]
			//   1. M = M + unpaired
			for i, state := range *beamstepM {
				if j < seq_length-1 {
					var newscore float64 = score_multi_unpaired(j+1, j+1) + state.score
					if bestM[j+1][i] == nil {
						bestM[j+1][i] = newState()
					}
					update_if_better(bestM[j+1][i], newscore, MANNER_M_eq_M_plus_U)
					nos_M++
				}
			}
		}
		//printf("M at %d\n", j); fflush(stdout);

		// beam of C
		{
			// C = C + U
			if j < seq_length-1 {
				var newscore float64 = score_external_unpaired(j+1, j+1) + beamstepC.score
				if bestC[j+1] == nil {
					bestC[j+1] = newState()
				}
				update_if_better(bestC[j+1], newscore, MANNER_C_eq_C_plus_U)
				nos_C++
			}
		}
		//printf("C at %d\n", j); fflush(stdout);
	} // end of for-loo j

	var viterbi *State = bestC[seq_length-1]

	// char result[seq_length + 1];
	result := get_parentheses(sequence)

	// gettimeofday(&parse_endtime, NULL);
	// double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec - parse_starttime.tv_usec) / 1000000.0;

	// var nos_tot uint64 = nos_H + nos_P + nos_M2 + nos_Multi + nos_M + nos_C
	// if (is_verbose)
	// {
	// 		printf("Parse Time: %f len: %d score %f #states %lu H %lu P %lu M2 %lu Multi %lu M %lu C %lu\n",
	// 						parse_elapsed_time, seq_length, double(viterbi.score), nos_tot,
	// 						nos_H, nos_P, nos_M2, nos_Multi, nos_M, nos_C);
	// }

	//fmt.Printf("result: %v\nscore: %v\n", result, viterbi.score)
	return result, viterbi.score
	// return {string(result), viterbi.score, nos_tot};
}

// vivek: func beam_prune in original source
func BeamPrune(beamstep *map[int]*State) float64 {
	scores = make([]Pair, 0)
	for i, cand := range *beamstep {
		k := i - 1
		var newscore float64
		if (k >= 0) && (bestC[k].score == VALUE_MIN) {
			newscore = VALUE_MIN
		} else if k >= 0 {
			newscore = bestC[k].score + cand.score
		} else {
			newscore = cand.score
		}
		scores = append(scores, Pair{newscore, i})
	}

	if len(scores) <= beam_size {
		return VALUE_MIN
	}

	var threshold float64 = quickselect(0, len(scores)-1, len(scores)-beam_size)

	for _, pair := range scores {
		if pair.first.(float64) < threshold {
			delete(*beamstep, pair.second.(int))
		}
	}

	return threshold
}

// ------------- nucs based scores -------------

func score_hairpin(i, j, nuci, nuci1, nucj_1, nucj int) float64 {
	return hairpin_length[Min(j-i-1, HAIRPIN_MAX_LEN)] +
		score_junction_B(i, j, nuci, nuci1, nucj_1, nucj)
}

func score_junction_B(i, j, nuci, nuci1, nucj_1, nucj int) float64 {
	return helix_closing_score(nuci, nucj) + terminal_mismatch_score(nuci, nuci1, nucj_1, nucj)
}

// parameters: nucs[i], nucs[j]
func helix_closing_score(nuci, nucj int) float64 {
	return helix_closing[nuci*NOTON+nucj]
}

// parameters: nucs[i], nucs[i+1], nucs[j-1], nucs[j]
func terminal_mismatch_score(nuci, nuci1, nucj_1, nucj int) float64 {
	return terminal_mismatch[nuci*NOTONT+nucj*NOTOND+nuci1*NOTON+nucj_1]
}

// in-place quick-select
func quickselect(lower int, upper int, k int) float64 {

	if lower == upper {
		return scores[lower].first.(float64)
	}
	var split int = quickselect_partition(lower, upper)
	var length int = split - lower + 1

	if length == k {
		return scores[split].first.(float64)
	} else if k < length {
		return quickselect(lower, split-1, k)
	} else {
		return quickselect(split+1, upper, k-length)
	}
}

func quickselect_partition(lower int, upper int) int {
	var pivot float64 = scores[upper].first.(float64)
	for lower < upper {
		for scores[lower].first.(float64) < pivot {
			lower++
		}
		for scores[upper].first.(float64) > pivot {
			upper--
		}
		if scores[lower].first == scores[upper].first {
			lower++
		} else if lower < upper {
			scores[lower], scores[upper] = scores[upper], scores[lower]
		}
	}
	return upper
}

func score_M1(i, j, k, nuci_1, nuci, nuck, nuck1, len int) float64 {
	return score_junction_A(k, i, nuck, nuck1, nuci_1, nuci, len) +
		score_multi_unpaired(k+1, j) + base_pair_score(nuci, nuck) + multi_paired
}

func score_junction_A(i, j, nuci, nuci1, nucj_1, nucj, len int) float64 {
	var result float64 = helix_closing_score(nuci, nucj)
	if i < len-1 {
		result += dangle_left_score(nuci, nuci1, nucj)
	}
	if j > 0 {
		result += dangle_right_score(nuci, nucj_1, nucj)
	}
	return result
}

func score_multi_unpaired(i, j int) float64 {
	return (float64(j) - float64(i) + 1) * multi_unpaired
}

// parameters: nucs[i], nucs[j]
func base_pair_score(nuci, nucj int) float64 {
	return base_pair[nucj*NOTON+nuci]
}

// parameters: nucs[i], nucs[i+1], nucs[j]
func dangle_left_score(nuci, nuci1, nucj int) float64 {
	return dangle_left[nuci*NOTOND+nucj*NOTON+nuci1]
}

// parameters: nucs[i], nucs[j-1], nucs[j]
func dangle_right_score(nuci, nucj_1, nucj int) float64 {
	return dangle_right[nuci*NOTOND+nucj*NOTON+nucj_1]
}

func score_external_paired(i, j, nuci_1, nuci, nucj, nucj1, len int) float64 {
	return score_junction_A(j, i, nucj, nucj1, nuci_1, nuci, len) +
		external_paired + base_pair_score(nuci, nucj)
}

func score_helix(nuci, nuci1, nucj_1, nucj int) float64 {
	return helix_stacking_score(nuci, nuci1, nucj_1, nucj) + base_pair_score(nuci1, nucj_1)
}

// parameters: nucs[i], nucs[i+1], nucs[j-1], nucs[j]
func helix_stacking_score(nuci, nuci1, nucj_1, nucj int) float64 {
	return helix_stacking[nuci*NOTONT+nucj*NOTOND+nuci1*NOTON+nucj_1]
}

func sortM(threshold float64, beamstep *map[int]*State, sorted_stepM []Pair) {
	sorted_stepM = make([]Pair, 0)
	if threshold == VALUE_MIN {
		// no beam pruning before, so scores vector not usable
		for i, cand := range *beamstep {
			k := i - 1
			var newscore float64
			if k >= 0 {
				newscore = bestC[k].score + cand.score
			} else {
				newscore = cand.score
			}
			sorted_stepM = append(sorted_stepM, Pair{newscore, i})
		}
	} else {
		for _, pair := range scores {
			if pair.first.(float64) >= threshold {
				sorted_stepM = append(sorted_stepM, pair)
			}
		}
	}

	sort.Slice(sorted_stepM, func(i, j int) bool {
		return sorted_stepM[i].first.(float64) > sorted_stepM[j].first.(float64)
	})
}

func score_external_unpaired(i, j int) float64 {
	return (float64(j) - float64(i) + 1) * external_unpaired
}

type Tuple struct {
	first, second, third interface{}
}

func get_parentheses(seq string) string {
	// is_verbose := true

	seq_length := len(seq)
	result := make([]rune, seq_length)
	for i := 0; i < seq_length; i++ {
		result[i] = '.'
	}

	var stk []Tuple
	// stack<tuple<int, int, State>> stk;
	stk = append(stk, Tuple{0, seq_length - 1, bestC[seq_length-1]})

	// verbose stuff
	// var multi_todo []Pair
	// var mbp map[int]int // multi bp

	// var total_energy float64 = .0
	// var external_energy float64 = .0

	for len(stk) != 0 {
		last_elem_idx := len(stk) - 1
		top := stk[last_elem_idx]
		i, j, state := top.first.(int), top.second.(int), top.third.(*State)

		// pop off stack
		// stk[last_elem_idx] = nil
		stk = stk[:last_elem_idx]

		var k, p, q int

		switch state.manner {
		case MANNER_H:
			// this state should not be traced
		case MANNER_HAIRPIN:
			{
				result[i] = '('
				result[j] = ')'
				// if (is_verbose) {
				// 	var tetra_hex_tri int = -1
				// 	if (j - i - 1 == 4) {
				// 		// 6:tetra
				// 		tetra_hex_tri = if_tetraloops[i]
				// 	} else if (j - i - 1 == 6) {
				// 		// 8:hexa
				// 		tetra_hex_tri = if_hexaloops[i]
				// 	} else if (j - i - 1 == 3) {
				// 		// 5:tri
				// 		tetra_hex_tri = if_triloops[i]
				// 		var nuci, nucj int = nucs[i], nucs[j]

				// 		var nuci1 int
				// 		if (i + 1) < seq_length {
				// 			nuci1 = nucs[i + 1]
				// 		} else {
				// 			nuci1 = -1
				// 		}

				// 		var nucj_1 int
				// 		if (j - 1) > -1 {
				// 			nucj_1 = nucs[j - 1]
				// 		} else {
				// 			nucj_1 = -1
				// 		}

				// 		var newscore float64 = -v_score_hairpin(i, j, nuci, nuci1, nucj_1, nucj, tetra_hex_tri);
				// 		fmt.Printf("Hairpin loop ( %d, %d) %c%c : %.2f\n", i + 1, j + 1, seq[i], seq[j], newscore / -100.0);
				// 		total_energy += newscore;
				// 	}
				// }
			}
		case MANNER_SINGLE:
			{
				result[i] = '('
				result[j] = ')'
				p = i + int(state.trace.paddings.l1)
				q = j - state.trace.paddings.l2
				stk = append(stk, Tuple{p, q, bestP[q][p]})
				// if (is_verbose)
				// {
				// 		int nuci = nucs[i], nuci1 = nucs[i + 1], nucj_1 = nucs[j - 1], nucj = nucs[j];
				// 		int nucp_1 = nucs[p - 1], nucp = nucs[p], nucq = nucs[q], nucq1 = nucs[q + 1];

				// 		value_type newscore = -v_score_single(i, j, p, q, nuci, nuci1, nucj_1, nucj,
				// 																					nucp_1, nucp, nucq, nucq1);
				// 		printf("Interior loop ( %d, %d) %c%c; ( %d, %d) %c%c : %.2f\n", i + 1, j + 1, seq[i], seq[j], p + 1, q + 1, seq[p], seq[q], newscore / -100.0);
				// 		total_energy += newscore;
				// }
			}
		case MANNER_HELIX:
			{
				result[i] = '('
				result[j] = ')'
				stk = append(stk, Tuple{i + 1, j - 1, bestP[j-1][i+1]})
				// if (is_verbose)
				// {
				// 		p = i + 1;
				// 		q = j - 1;
				// 		int nuci = nucs[i], nuci1 = nucs[i + 1], nucj_1 = nucs[j - 1], nucj = nucs[j];
				// 		int nucp_1 = nucs[p - 1], nucp = nucs[p], nucq = nucs[q], nucq1 = nucs[q + 1];

				// 		value_type newscore = -v_score_single(i, j, p, q, nuci, nuci1, nucj_1, nucj,
				// 																					nucp_1, nucp, nucq, nucq1);
				// 		printf("Interior loop ( %d, %d) %c%c; ( %d, %d) %c%c : %.2f\n", i + 1, j + 1, seq[i], seq[j], p + 1, q + 1, seq[p], seq[q], newscore / -100.0);
				// 		total_energy += newscore;
				// }
			}
		case MANNER_MULTI:
			p = i + int(state.trace.paddings.l1)
			q = j - state.trace.paddings.l2
			stk = append(stk, Tuple{p, q, bestM2[q][p]})
			// stk.push(make_tuple(p, q, bestM2[q][p]));
		case MANNER_MULTI_eq_MULTI_plus_U:
			p = i + int(state.trace.paddings.l1)
			q = j - state.trace.paddings.l2
			stk = append(stk, Tuple{p, q, bestM2[q][p]})
		case MANNER_P_eq_MULTI:
			result[i] = '('
			result[j] = ')'
			stk = append(stk, Tuple{i, j, bestMulti[j][i]})
			// if (is_verbose) {
			// 		multi_todo.push_back(make_pair(i, j));
			// }
		case MANNER_M2_eq_M_plus_P:
			k = state.trace.split
			stk = append(stk, Tuple{i, k, bestM[k][i]})
			stk = append(stk, Tuple{k + 1, j, bestP[j][k+1]})
			// if (is_verbose) {
			// 	mbp[k + 1] = j
			// }
			// break
		case MANNER_M_eq_M2:
			stk = append(stk, Tuple{i, j, bestM2[j][i]})
			// stk.push(make_tuple(i, j, bestM2[j][i]));
			// break
		case MANNER_M_eq_M_plus_U:
			stk = append(stk, Tuple{i, j - 1, bestM[j-1][i]})
			// break
		case MANNER_M_eq_P:
			stk = append(stk, Tuple{i, j, bestP[j][i]})
			// if (is_verbose) {
			// 	mbp[i] = j;
			// }
		case MANNER_C_eq_C_plus_U:
			k = j - 1
			if k != -1 {
				stk = append(stk, Tuple{0, k, bestC[k]})
			}
		// if (is_verbose)
		// 		external_energy += -v_score_external_unpaired(0, 0); // zero at this moment
		case MANNER_C_eq_C_plus_P:
			{
				k = state.trace.split
				if k != -1 {
					stk = append(stk, Tuple{0, k, bestC[k]})
					stk = append(stk, Tuple{k + 1, j, bestP[j][k+1]})
				} else {
					stk = append(stk, Tuple{i, j, bestP[j][i]})
				}
				// if (is_verbose)
				// {
				// 		int nuck = k > -1 ? nucs[k] : -1;
				// 		int nuck1 = nucs[k + 1], nucj = nucs[j];
				// 		int nucj1 = (j + 1) < seq_length ? nucs[j + 1] : -1;
				// 		external_energy += -v_score_external_paired(k + 1, j, nuck, nuck1,
				// 																								nucj, nucj1, seq_length);
				// }
			}
		default: // MANNER_NONE or other cases
			fmt.Printf("wrong manner at %d, %d: manner %d\n", i, j, state.manner)
			// fflush(stdout);
			// assert(false);
			panic("wrong manner")
		}
	}

	// 	if (is_verbose)
	// 	{
	// 			for (auto item : multi_todo)
	// 			{
	// 					int i = item.first;
	// 					int j = item.second;
	// 					int nuci = nucs[i], nuci1 = nucs[i + 1], nucj_1 = nucs[j - 1], nucj = nucs[j];
	// 					value_type multi_energy = -v_score_multi(i, j, nuci, nuci1, nucj_1, nucj, seq_length);
	// 					int num_unpaired = 0;
	// 					for (int k = i + 1; k < j; ++k)
	// 					{
	// 							if (result[k] == '.')
	// 									num_unpaired += 1;
	// 							else if (result[k] == '(')
	// 							{
	// 									int p = k, q = mbp[k];
	// 									int nucp_1 = nucs[p - 1], nucp = nucs[p], nucq = nucs[q], nucq1 = nucs[q + 1];

	// 									multi_energy += -v_score_M1(p, q, q, nucp_1, nucp, nucq, nucq1, seq_length);
	// 									k = q;
	// 							}
	// 					}
	// 					multi_energy += -v_score_multi_unpaired(1, num_unpaired);

	// 					printf("Multi loop ( %d, %d) %c%c : %.2f\n", i + 1, j + 1, seq[i], seq[j], multi_energy / -100.0);
	// 					total_energy += multi_energy;
	// 			}

	// 			printf("External loop : %.2f\n", external_energy / -100.0);
	// 			total_energy += external_energy;

	// #ifndef lv
	// 			printf("Energy(kcal/mol): %.2f\n", total_energy / -100.0);
	// #endif
	// 	}

	return string(result)
}

func NUM_TO_NUC(x int) int {
	switch x {
	case -1:
		return -1
	case 4:
		return 0
	default:
		return x + 1
	}
}

func NUM_TO_PAIR(x, y int) int {
	switch x {
	case 0:
		if y == 3 {
			return 5
		} else {
			return 0
		}
	case 1:
		if y == 2 {
			return 1
		} else {
			return 0
		}
	case 2:
		if y == 1 {
			return 2
		} else {
			if y == 3 {
				return 3
			} else {
				return 0
			}
		}
	case 3:
		if y == 2 {
			return 4
		} else if y == 0 {
			return 6
		} else {
			return 0
		}
	default:
		return 0
	}
}

func score_multi(i, j, nuci, nuci1, nucj_1, nucj, len int) float64 {
	return score_junction_A(i, j, nuci, nuci1, nucj_1, nucj, len) +
		multi_paired + multi_base
}

func score_single_without_junctionB(i, j, p, q, nucp_1, nucp, nucq, nucq1 int) float64 {
	var l1, l2 int = p - i - 1, j - q - 1
	return cache_single[l1][l2] + base_pair_score(nucp, nucq) +
		score_single_nuc(i, j, p, q, nucp_1, nucq1)
}

func score_single_nuc(i, j, p, q, nucp_1, nucq1 int) float64 {
	var l1, l2 int = p - i - 1, j - q - 1
	if l1 == 0 && l2 == 1 {
		return bulge_nuc_score(nucq1)
	} else if l1 == 1 && l2 == 0 {
		return bulge_nuc_score(nucp_1)
	} else if l1 == 1 && l2 == 1 {
		return internal_nuc_score(nucp_1, nucq1)
	}
	return 0
}

// parameter: nucs[i]
func bulge_nuc_score(nuci int) float64 {
	return bulge_0x1_nucleotides[nuci]
}

// parameters: nucs[i], nucs[j]
func internal_nuc_score(nuci, nucj int) float64 {
	return internal_1x1_nucleotides[nuci*NOTON+nucj]
}

// process_record

// record.sequence = rec_sequence;
// record.rest = rec_rest;
// record.options = opt; // default options
// var rec_rest []string // constians input lines
// vrna_file_fasta_read_record(&rec_id,
// 	&rec_sequence,
// 	&rec_rest,
// 	input_stream,
// 	read_opt);

type vrna_fold_compound_t struct {
	length   uint64 /**<  @brief  The length of the sequence (or sequence alignment) */
	cutpoint int    /*  @brief  The position of the (cofold) cutpoint within the provided sequence.
	* If there is no cutpoint, this field will be set to -1
	 */

	strand_number []uint64 /**<  @brief  The strand number a particular nucleotide is associated with */
	strand_order  []uint64 /**<  @brief  The strand order, i.e. permutation of current concatenated sequence */
	strand_start  []uint64 /**<  @brief  The start position of a particular strand within the current concatenated sequence */
	strand_end    []uint64 /**<  @brief  The end (last) position of a particular strand within the current concatenated sequence */

	strands     uint64
	nucleotides []vrna_seq_t
	alignment   *vrna_msa_t

	hc *vrna_hc_t /**<  @brief  The hard constraints data structure used for structure prediction */

	matrices *vrna_mx_mfe_t /**<  @brief  The MFE DP matrices */
	// exp_matrices *vrna_mx_pf_t  /**<  @brief  The PF DP matrices  */

	params *vrna_param_t /**<  @brief  The precomputed free energy contributions for each type of loop */
	// exp_params *vrna_exp_param_t    /**<  @brief  The precomputed free energy contributions as Boltzmann factors  */

	iindx []int /**<  @brief  DP matrix accessor  */
	jindx []int /**<  @brief  DP matrix accessor  */

	/**
	 *  @}
	 *
	 *  @name Secondary Structure Decomposition (grammar) related data fields
	 *  @{
	 */

	/* data structure to adjust additional structural domains, such as G-quadruplexes */
	// domains_struc *vrna_sd_t             /**<  @brief  Additional structured domains */

	/* data structure to adjust additional contributions to unpaired stretches, e.g. due to protein binding */
	// domains_up     *vrna_ud_                /**<  @brief  Additional unstructured domains */

	/**
	 *  @}
	 */

	/**
	 *  @name Data fields available for single/hybrid structure prediction
	 *  @{
	 */
	sequence string /**<  @brief  The input sequence string
	 *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_SINGLE @endverbatim
	 */
	sequence_encoding []int /**<  @brief  Numerical encoding of the input sequence
	 *    @see    vrna_sequence_encode()
	 *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_SINGLE @endverbatim
	 */
	sequence_encoding2 []int
	ptype              string /**<  @brief  Pair type array
	 *
	 *    Contains the numerical encoding of the pair type for each pair (i,j) used
	 *    in MFE, Partition function and Evaluation computations.
	 *    @note This array is always indexed via jindx, in contrast to previously
	 *    different indexing between mfe and pf variants!
	 *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_SINGLE @endverbatim
	 *    @see    vrna_idx_col_wise(), vrna_ptypes()
	 */
	ptype_pf_compat string /**<  @brief  ptype array indexed via iindx
	 *    @deprecated  This attribute will vanish in the future!
	 *    It's meant for backward compatibility only!
	 *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_SINGLE @endverbatim
	 */
	sc *vrna_sc_t /**<  @brief  The soft constraints for usage in structure prediction and evaluation
	 *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_SINGLE @endverbatim
	 */

	/**
	 *  @}
	 */

	// /**
	//  *  @name Data fields for consensus structure prediction
	//  *  @{
	//  */
	//     sequences []string        /**<  @brief  The aligned sequences
	//                                        *    @note   The end of the alignment is indicated by a NULL pointer in the second dimension
	//                                        *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
	//                                        */
	//     n_seq uint64              /**<  @brief  The number of sequences in the alignment
	//                                        *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
	//                                        */
	//     cons_seq string          /**<  @brief  The consensus sequence of the aligned sequences
	//                                        *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
	//                                        */
	//     S_cons []int           /**<  @brief  Numerical encoding of the consensus sequence
	//                                        *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
	//                                        */
	//     S   [][]int            /**<  @brief  Numerical encoding of the sequences in the alignment
	//                                        *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
	//                                        */
	//     S5  [][]int             /**<  @brief    S5[s][i] holds next base 5' of i in sequence s
	//                                        *    @warning  Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
	//                                        */
	//     S3  [][]int             /**<  @brief    Sl[s][i] holds next base 3' of i in sequence s
	//                                        *    @warning  Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
	//                                        */
	// 		Ss []string
	// 		a2s [][]uint64
	//     pscore []int              /**<  @brief  Precomputed array of pair types expressed as pairing scores
	//                                          *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
	//                                          */
	//     pscore_local [][]int       /**<  @brief  Precomputed array of pair types expressed as pairing scores
	//                                          *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
	//                                          */
	//     pscore_pf_compat []int    /**<  @brief  Precomputed array of pair types expressed as pairing scores indexed via iindx
	//                                          *    @deprecated  This attribute will vanish in the future!
	//                                          *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
	//                                          */
	//     scs []*vrna_sc_t                /**<  @brief  A set of soft constraints (for each sequence in the alignment)
	//                                          *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
	//                                          */
	// 		oldAliEn int

	/**
	 *  @}
	 */

	/**
	 *  @name Additional data fields for Distance Class Partitioning
	 *
	 *  These data fields are typically populated with meaningful data only if used in the context of Distance Class Partitioning
	 *  @{
	 */
	maxD1         uint64 /**<  @brief  Maximum allowed base pair distance to first reference */
	maxD2         uint64 /**<  @brief  Maximum allowed base pair distance to second reference */
	reference_pt1 []int  /**<  @brief  A pairtable of the first reference structure */
	reference_pt2 []int  /**<  @brief  A pairtable of the second reference structure */

	referenceBPs1 []uint64 /**<  @brief  Matrix containing number of basepairs of reference structure1 in interval [i,j] */
	referenceBPs2 []uint64 /**<  @brief  Matrix containing number of basepairs of reference structure2 in interval [i,j] */
	bpdist        []uint64 /**<  @brief  Matrix containing base pair distance of reference structure 1 and 2 on interval [i,j] */

	mm1 []uint64 /**<  @brief  Maximum matching matrix, reference struct 1 disallowed */
	mm2 []uint64 /**<  @brief  Maximum matching matrix, reference struct 2 disallowed */

	/**
	 *  @}
	 */

	/**
	 *  @name Additional data fields for local folding
	 *
	 *  These data fields are typically populated with meaningful data only if used in the context of local folding
	 *  @{
	 */
	window_size int      /**<  @brief  window size for local folding sliding window approach */
	ptype_local []string /**<  @brief  Pair type array (for local folding) */

	/**
	 *  @}
	 */
}

func default_vrna_fold_compound_t() *vrna_fold_compound_t {
	return &vrna_fold_compound_t{
		length:        0,
		strands:       0,
		cutpoint:      -1,
		strand_number: nil,
		strand_order:  nil,
		strand_start:  nil,
		strand_end:    nil,
		nucleotides:   nil,
		alignment:     nil,

		hc:       nil,
		matrices: nil,
		// exp_matrices: nil,
		params: nil,
		// exp_params:   nil,
		iindx: nil,
		jindx: nil,

		// stat_cb:      nil,
		// auxdata:      nil,
		// free_auxdata: nil,

		// domains_struc: nil,
		// domains_up:    nil,
		// aux_grammar:   nil,

		sequence:           "",
		sequence_encoding:  nil,
		sequence_encoding2: nil,
		ptype:              "",
		ptype_pf_compat:    "",
		sc:                 nil,

		// axD1:          0,
		// axD2:          0,
		reference_pt1: nil,
		reference_pt2: nil,
		referenceBPs1: nil,
		referenceBPs2: nil,
		bpdist:        nil,
		mm1:           nil,
		mm2:           nil,

		window_size: -1,
		ptype_local: nil,
	}
}

type vrna_seq_t struct {
	sequence  string /**< @brief The string representation of the sequence */
	encoding  []int  /**< @brief The integer representation of the sequence */
	encoding5 []int
	encoding3 []int
	length    uint64 /**< @brief The length of the sequence */
}

type vrna_msa_t struct {
	n_seq        uint64
	sequences    []vrna_seq_t
	gapfree_seq  []string
	gapfree_size []int    /* for MAF alignment coordinates */
	genome_size  []uint64 /* for MAF alignment coordinates */
	start        []uint64 /* for MAF alignment coordinates */
	orientation  string   /* for MAF alignment coordinates */
	a2s          [][]uint64
}

type vrna_hc_t struct {
	n uint64

	state rune

	mx []rune

	matrix_local [][]rune

	up_ext []int /**<  @brief  A linear array that holds the number of allowed
	 *            unpaired nucleotides in an exterior loop
	 */
	up_hp []int /**<  @brief  A linear array that holds the number of allowed
	 *            unpaired nucleotides in a hairpin loop
	 */
	up_int []int /**<  @brief  A linear array that holds the number of allowed
	 *            unpaired nucleotides in an interior loop
	 */
	up_ml []int /**<  @brief  A linear array that holds the number of allowed
	 *            unpaired nucleotides in a multi branched loop
	 */
}

/**
 *  @brief  Minimum Free Energy (MFE) Dynamic Programming (DP) matrices data structure required within the #vrna_fold_compound_t
 */
type vrna_mx_mfe_t struct {
	/** @name Common fields for MFE matrices
	 *  @{
	 */

	length uint64 /**<  @brief  Length of the sequence, therefore an indicator of the size of the DP matrices */
	/**
	 *  @}
	 */

	/** @name Default DP matrices
	 *  @note These data fields are available if
	 *        @code vrna_mx_mfe_t.type == VRNA_MX_DEFAULT @endcode
	 * @{
	 */
	c   []int /**<  @brief  Energy array, given that i-j pair */
	f5  []int /**<  @brief  Energy of 5' end */
	f3  []int /**<  @brief  Energy of 3' end */
	fc  []int /**<  @brief  Energy from i to cutpoint (and vice versa if i>cut) */
	fML []int /**<  @brief  Multi-loop auxiliary energy array */
	fM1 []int /**<  @brief  Second ML array, only for unique multibrnach loop decomposition */
	fM2 []int /**<  @brief  Energy for a multibranch loop region with exactly two stems, extending to 3' end */
	ggg []int /**<  @brief  Energies of g-quadruplexes */
	Fc  int   /**<  @brief  Minimum Free Energy of entire circular RNA */
	FcH int
	FcI int
	FcM int
	/**
	 * @}
	 */

	/** @name Local Folding DP matrices using window approach
	 *  @note These data fields are available if
	 *        @code vrna_mx_mfe_t.type == VRNA_MX_WINDOW @endcode
	 * @{
	 */
	c_local   [][]int /**<  @brief  Energy array, given that i-j pair */
	f3_local  []int   /**<  @brief  Energy of 5' end */
	fML_local [][]int /**<  @brief  Multi-loop auxiliary energy array */
	ggg_local [][]int /**<  @brief  Energies of g-quadruplexes */
	/**
	 * @}
	 */

	/** @name Distance Class DP matrices
	 *  @note These data fields are available if
	 *        @code vrna_mx_mfe_t.type == VRNA_MX_2DFOLD @endcode
	 * @{
	 */
	E_F5     [][][]int
	l_min_F5 [][]int
	l_max_F5 [][]int
	k_min_F5 []int
	k_max_F5 []int

	E_F3     [][][]int
	l_min_F3 [][]int
	l_max_F3 [][]int
	k_min_F3 []int
	k_max_F3 []int

	E_C     [][][]int
	l_min_C [][]int
	l_max_C [][]int
	k_min_C []int
	k_max_C []int

	E_M     [][][]int
	l_min_M [][]int
	l_max_M [][]int
	k_min_M []int
	k_max_M []int

	E_M1     [][][]int
	l_min_M1 [][]int
	l_max_M1 [][]int
	k_min_M1 []int
	k_max_M1 []int

	E_M2     [][][]int
	l_min_M2 [][]int
	l_max_M2 [][]int
	k_min_M2 []int
	k_max_M2 []int

	E_Fc     [][]int
	l_min_Fc []int
	l_max_Fc []int
	k_min_Fc int
	k_max_Fc int

	E_FcH     [][]int
	l_min_FcH []int
	l_max_FcH []int
	k_min_FcH int
	k_max_FcH int

	E_FcI     [][]int
	l_min_FcI []int
	l_max_FcI []int
	k_min_FcI int
	k_max_FcI int

	E_FcM     [][]int
	l_min_FcM []int
	l_max_FcM []int
	k_min_FcM int
	k_max_FcM int

	/* auxilary arrays for remaining set of coarse graining (k,l) > (k_max, l_max) */
	E_F5_rem []int
	E_F3_rem []int
	E_C_rem  []int
	E_M_rem  []int
	E_M1_rem []int
	E_M2_rem []int

	E_Fc_rem  int
	E_FcH_rem int
	E_FcI_rem int
	E_FcM_rem int

	// #ifdef COUNT_STATES
	//   unsigned long ***N_F5;
	//   unsigned long ***N_C;
	//   unsigned long ***N_M;
	//   unsigned long ***N_M1;
	// #endif

	/**
	 * @}
	 */
}

/**
 *  @brief  The soft constraints data structure
 *
 *  @ingroup soft_constraints
 */
type vrna_sc_t struct {
	n uint64

	state rune

	energy_up     [][]int     /**<  @brief Energy contribution for stretches of unpaired nucleotides */
	exp_energy_up [][]float64 /**<  @brief Boltzmann Factors of the energy contributions for unpaired sequence stretches */

	up_storage []int                    /**<  @brief  Storage container for energy contributions per unpaired nucleotide */
	bp_storage [][]vrna_sc_bp_storage_t /**<  @brief  Storage container for energy contributions per base pair */

	energy_bp     []int     /**<  @brief Energy contribution for base pairs */
	exp_energy_bp []float64 /**<  @brief Boltzmann Factors of the energy contribution for base pairs */

	energy_bp_local     [][]int     /**<  @brief Energy contribution for base pairs (sliding window approach) */
	exp_energy_bp_local [][]float64 /**<  @brief Boltzmann Factors of the energy contribution for base pairs (sliding window approach) */

	energy_stack     []int     /**<  @brief Pseudo Energy contribution per base pair involved in a stack */
	exp_energy_stack []float64 /**<  @brief Boltzmann weighted pseudo energy contribution per nucleotide involved in a stack */

	// /* generic soft contraints below */
	// vrna_callback_sc_energy     *f;     /**<  @brief  A function pointer used for pseudo
	//                                      *            energy contribution in MFE calculations
	//                                      *    @see    vrna_sc_add_f()
	//                                      */

	// vrna_callback_sc_backtrack  *bt;    /**<  @brief  A function pointer used to obtain backtraced
	//                                      *            base pairs in loop regions that were altered
	//                                      *            by soft constrained pseudo energy contributions
	//                                      *    @see    vrna_sc_add_bt()
	//                                      */

	// vrna_callback_sc_exp_energy *exp_f; /**<  @brief  A function pointer used for pseudo energy
	//                                      *            contribution boltzmann factors in PF
	//                                      *            calculations
	//                                      *    @see    vrna_sc_add_exp_f()
	//                                      */

	// void                        *data;  /**<  @brief  A pointer to the data object provided for
	//                                      *            for pseudo energy contribution functions of the
	//                                      *            generic soft constraints feature
	//                                      */
	// vrna_callback_free_auxdata  *free_data;
}

/**
 *  @brief  A base pair constraint
 */

// /**
//  *  @brief  Partition function (PF) Dynamic Programming (DP) matrices data structure required within the #vrna_fold_compound_t
//  */
// type vrna_mx_pf_t struct {
//   /** @name Common fields for DP matrices
//    *  @{
//    */

//   length uint64
//   scale []float64
//   expMLbase []float64

//   /**
//    *  @}
//    */

//   /** @name Default PF matrices
//    *  @note These data fields are available if
//    *        @code vrna_mx_pf_t.type == VRNA_MX_DEFAULT @endcode
//    *  @{
//    */
//   q []float64
//   qb []float64
//   qm []float64
//   qm1 []float64
//   probs []float64
//   q1k []float64
//   qln []float64
//   G []float64

//   qo float64
//   qm2 []float64
//   qho float64
//   qio float64
//   qmo float64

//   /**
//    *  @}
//    */

// /** @name Local Folding DP matrices using window approach
//    *  @note These data fields are available if
//    *        @code vrna_mx_mfe_t.type == VRNA_MX_WINDOW @endcode
//    * @{
//    */
// 	 q_local [][]float64
// 	 qb_local [][]float64
// 	 qm_local [][]float64
// 	 pR [][]float64
// 	 qm2_local [][]float64
// 	 QI5 [][]float64
// 	 q2l [][]float64
// 	 qmb [][]float64
// 	 G_local [][]float64
// 	 /**
// 		*  @}
// 		*/

// 	 /** @name Distance Class DP matrices
// 		*  @note These data fields are available if
// 		*        @code vrna_mx_pf_t.type == VRNA_MX_2DFOLD @endcode
// 		*  @{
// 		*/
// 	 Q [][][]float64
// 	 l_min_Q [][]int
// 	 l_max_Q [][]int
// 	 int *k_min_Q;
// 	 int *k_max_Q;

// 	 Q_B [][][]float64
// 	 l_min_Q_B [][]int
// 	 l_max_Q_B [][]int
// 	 int *k_min_Q_B;
// 	 int *k_max_Q_B;

// 	 Q_M [][][]float64
// 	 l_min_Q_M [][]int
// 	 l_max_Q_M [][]int
// 	 int *k_min_Q_M;
// 	 int *k_max_Q_M;

// 	 Q_M1 [][][]float64
// 	 l_min_Q_M1 [][]int
// 	 l_max_Q_M1 [][]int
// 	 int *k_min_Q_M1;
// 	 int *k_max_Q_M1;

// 	 Q_M2 [][][]float64
// 	 l_min_Q_M2 [][]int
// 	 l_max_Q_M2 [][]int
// 	 int *k_min_Q_M2;
// 	 int *k_max_Q_M2;

// 	 Q_c [][]float64
// 	 int *l_min_Q_c;
// 	 int *l_max_Q_c;
// 	 int k_min_Q_c;
// 	 int k_max_Q_c;

// 	 Q_cH [][]float64
// 	 int *l_min_Q_cH;
// 	 int *l_max_Q_cH;
// 	 int k_min_Q_cH;
// 	 int k_max_Q_cH;

// 	 Q_cI [][]float64
// 	 int *l_min_Q_cI;
// 	 int *l_max_Q_cI;
// 	 int k_min_Q_cI;
// 	 int k_max_Q_cI;

// 	 Q_cM [][]float64
// 	 int *l_min_Q_cM;
// 	 int *l_max_Q_cM;
// 	 int k_min_Q_cM;
// 	 int k_max_Q_cM;

// 	 /* auxilary arrays for remaining set of coarse graining (k,l) > (k_max, l_max) */
// 	 FLT_OR_DBL *Q_rem;
// 	 FLT_OR_DBL *Q_B_rem;
// 	 FLT_OR_DBL *Q_M_rem;
// 	 FLT_OR_DBL *Q_M1_rem;
// 	 FLT_OR_DBL *Q_M2_rem;

// 	 FLT_OR_DBL Q_c_rem;
// 	 FLT_OR_DBL Q_cH_rem;
// 	 FLT_OR_DBL Q_cI_rem;
// 	 FLT_OR_DBL Q_cM_rem;
// 	 /**
// 		*  @}
// 		*/

//  #ifndef VRNA_DISABLE_C11_FEATURES
// 	 /* C11 support for unnamed unions/structs */
//  };
//  };
//  #endif
//  };

func CalculateMfe(seq, structure string) (float64, error) {
	if len(seq) != len(structure) {
		return 0, errors.New("Length of sequence (%v) != Length of structure (%v)", len(seq), len(structure))
	} else if len(seq) == 0 {
		return 0, errors.New("Lengths of sequence and structure cannot be 0")
	}

	// fc is everything with NIL
	// md is the default model details
	// add_params(fc, &md, options);
	/*
	 * ALWAYS provide regular energy parameters
	 * remove previous parameters if present and they differ from current model
	 */

	fc := default_vrna_fold_compound_t()
	fc.params = vrna_params()
	fc.sequence = seq
	fc.length = len(seq)
	// fc := &vrna_fold_compound_t{params: vrna_params(), sequence: seq, structure: structure}

	// vrna_params_prepare(fc, options);
	sanitize_bp_span(fc)
	set_fold_compound(fc, sequence)

	// tmp is structure

	energy := vrna_eval_structure_cstr(fc, structure)
	return energy, nil
}

func vrna_params() *vrna_param_t {
	return get_scaled_params()
}

func RESCALE_dG_int(dG, dH, dT int) int {
	return (dH - (dH-dG)*dT)
}

func RESCALE_dG_float64(dG, dH, dT float64) float64 {
	return (dH - (dH-dG)*dT)
}

var (
	DEFAULT_TEMP float64 = 37.0
	NB_PAIRS     int     = 7  /** The number of distinguishable base pairs */
	MAXLOOP      int     = 30 /** The maximum loop length */

)

type vrna_md_t struct {
	temperature float64 /**<  @brief  The temperature used to scale the thermodynamic parameters */
	betaScale   float64 /**<  @brief  A scaling factor for the thermodynamic temperature of the Boltzmann factors */
	pf_smooth   int     /**<  @brief  A flat specifying whether energies in Boltzmann factors need to be smoothed */
	dangles     int     /**<  @brief  Specifies the dangle model used in any energy evaluation (0,1,2 or 3)
	 *
	 *    If set to 0 no stabilizing energies are assigned to bases adjacent to
	 *    helices in free ends and multiloops (so called dangling ends). Normally
	 *    (dangles = 1) dangling end energies are assigned only to unpaired
	 *    bases and a base cannot participate simultaneously in two dangling ends. In
	 *    the partition function algorithm vrna_pf() these checks are neglected.
	 *    To provide comparability between free energy minimization and partition function
	 *    algorithms, the default setting is 2.
	 *    This treatment of dangling ends gives more favorable energies to helices
	 *    directly adjacent to one another, which can be beneficial since such
	 *    helices often do engage in stabilizing interactions through co-axial
	 *    stacking.\n
	 *    If set to 3 co-axial stacking is explicitly included for
	 *    adjacent helices in multiloops. The option affects only mfe folding
	 *    and energy evaluation (vrna_mfe() and vrna_eval_structure()), as
	 *    well as suboptimal folding (vrna_subopt()) via re-evaluation of energies.
	 *    Co-axial stacking with one intervening mismatch is not considered so far.
	 *    @note   Some function do not implement all dangle model but only a subset of
	 *            (0,1,2,3). In particular, partition function algorithms can only handle
	 *            0 and 2. Read the documentation of the particular recurrences or
	 *            energy evaluation function for information about the provided dangle
	 *            model.
	 */
	special_hp     int      /**<  @brief  Include special hairpin contributions for tri, tetra and hexaloops */
	noLP           int      /**<  @brief  Only consider canonical structures, i.e. no 'lonely' base pairs */
	noGU           int      /**<  @brief  Do not allow GU pairs */
	noGUclosure    int      /**<  @brief  Do not allow loops to be closed by GU pair */
	logML          int      /**<  @brief  Use logarithmic scaling for multiloops */
	circ           int      /**<  @brief  Assume RNA to be circular instead of linear */
	gquad          int      /**<  @brief  Include G-quadruplexes in structure prediction */
	uniq_ML        int      /**<  @brief  Flag to ensure unique multi-branch loop decomposition during folding */
	energy_set     int      /**<  @brief  Specifies the energy set that defines set of compatible base pairs */
	backtrack      int      /**<  @brief  Specifies whether or not secondary structures should be backtraced */
	backtrack_type rune     /**<  @brief  Specifies in which matrix to backtrack */
	compute_bpp    int      /**<  @brief  Specifies whether or not backward recursions for base pair probability (bpp) computation will be performed */
	nonstandards   [64]rune /**<  @brief  contains allowed non standard bases */
	max_bp_span    int      /**<  @brief  maximum allowed base pair span */

	min_loop_size int /**<  @brief  Minimum size of hairpin loops
	 *    @note The default value for this field is #TURN, however, it may
	 *    be 0 in cofolding context.
	 */
	window_size int                             /**<  @brief  Size of the sliding window for locally optimal structure prediction */
	oldAliEn    int                             /**<  @brief  Use old alifold energy model */
	ribo        int                             /**<  @brief  Use ribosum scoring table in alifold energy model */
	cv_fact     float64                         /**<  @brief  Co-variance scaling factor for consensus structure prediction */
	nc_fact     float64                         /**<  @brief  Scaling factor to weight co-variance contributions of non-canonical pairs */
	sfact       float64                         /**<  @brief  Scaling factor for partition function scaling */
	rtype       [8]int                          /**<  @brief  Reverse base pair type array */
	alias       [MAXALPHA + 1]int               /**<  @brief  alias of an integer nucleotide representation */
	pair        [MAXALPHA + 1][MAXALPHA + 1]int /**<  @brief  Integer representation of a base pair */
}

var (
	// Default temperature for structure prediction and free energy evaluation in $^\circ C$
	VRNA_MODEL_DEFAULT_TEMPERATURE float64 = 37.0
	VRNA_MODEL_DEFAULT_PF_SMOOTH   int     = 1
	// Default dangling end model
	VRNA_MODEL_DEFAULT_DANGLES int = 2
	// Default model behavior for lookup of special tri-, tetra-, and hexa-loops
	VRNA_MODEL_DEFAULT_SPECIAL_HP int = 1
	// Default model behavior for so-called 'lonely pairs'
	VRNA_MODEL_DEFAULT_NO_LP int = 1
	// Default model behavior for G-U base pairs
	VRNA_MODEL_DEFAULT_NO_GU int = 0
	// Default model behavior for G-U base pairs closing a loop
	VRNA_MODEL_DEFAULT_NO_GU_CLOSURE int = 0
	// Default model behavior on how to evaluate the energy contribution of multi-branch loops
	VRNA_MODEL_DEFAULT_LOG_ML int = 0
	// Default model behavior to treat a molecule as a circular RNA (DNA)
	VRNA_MODEL_DEFAULT_CIRC int = 0
	// Default model behavior regarding the treatment of G-Quadruplexes
	VRNA_MODEL_DEFAULT_GQUAD int = 0
	// Default behavior of the model regarding unique multi-branch loop decomposition
	VRNA_MODEL_DEFAULT_UNIQ_ML int = 0
	// Default model behavior on which energy set to use
	VRNA_MODEL_DEFAULT_ENERGY_SET int = 0
	// Default model behavior with regards to backtracking of structures
	VRNA_MODEL_DEFAULT_BACKTRACK int = 1
	// Default model behavior on what type of backtracking to perform
	VRNA_MODEL_DEFAULT_BACKTRACK_TYPE rune = 'F'
	// Default model behavior with regards to computing base pair probabilities
	VRNA_MODEL_DEFAULT_COMPUTE_BPP int = 1
	// Default model behavior for the allowed maximum base pair span
	VRNA_MODEL_DEFAULT_MAX_BP_SPAN int = -1
	// The minimum loop length
	TURN int = 3
	// Default model behavior for the sliding window approach
	VRNA_MODEL_DEFAULT_WINDOW_SIZE int = -1
	// Default model behavior for consensus structure energy evaluation
	VRNA_MODEL_DEFAULT_ALI_OLD_EN int = 0
	// Default model behavior for consensus structure co-variance contribution assessment
	VRNA_MODEL_DEFAULT_ALI_RIBO int = 0
	// Default model behavior for weighting the co-variance score in consensus structure prediction
	VRNA_MODEL_DEFAULT_ALI_CV_FACT float64 = 1.0
	// Default model behavior for weighting the nucleotide conservation? in consensus structure prediction
	VRNA_MODEL_DEFAULT_ALI_NC_FACT float64 = 1.0
)

func DefaultModelDetails() *vrna_md_t {
	return &vrna_md_t{
		temperature:    VRNA_MODEL_DEFAULT_TEMPERATURE,
		betaScale:      1.,
		pf_smooth:      VRNA_MODEL_DEFAULT_PF_SMOOTH,
		dangles:        VRNA_MODEL_DEFAULT_DANGLES,
		special_hp:     VRNA_MODEL_DEFAULT_SPECIAL_HP,
		noLP:           VRNA_MODEL_DEFAULT_NO_LP,
		noGU:           VRNA_MODEL_DEFAULT_NO_GU,
		noGUclosure:    VRNA_MODEL_DEFAULT_NO_GU_CLOSURE,
		logML:          VRNA_MODEL_DEFAULT_LOG_ML,
		circ:           VRNA_MODEL_DEFAULT_CIRC,
		gquad:          VRNA_MODEL_DEFAULT_GQUAD,
		uniq_ML:        VRNA_MODEL_DEFAULT_UNIQ_ML,
		energy_set:     VRNA_MODEL_DEFAULT_ENERGY_SET,
		backtrack:      VRNA_MODEL_DEFAULT_BACKTRACK,
		backtrack_type: VRNA_MODEL_DEFAULT_BACKTRACK_TYPE,
		compute_bpp:    VRNA_MODEL_DEFAULT_COMPUTE_BPP,
		nonstandards:   []rune{0},
		max_bp_span:    VRNA_MODEL_DEFAULT_MAX_BP_SPAN,
		min_loop_size:  TURN,
		window_size:    VRNA_MODEL_DEFAULT_WINDOW_SIZE,
		oldAliEn:       VRNA_MODEL_DEFAULT_ALI_OLD_EN,
		ribo:           VRNA_MODEL_DEFAULT_ALI_RIBO,
		cv_fact:        VRNA_MODEL_DEFAULT_ALI_CV_FACT,
		nc_fact:        VRNA_MODEL_DEFAULT_ALI_NC_FACT,
		sfact:          1.07,
		rtype:          []int{0, 2, 1, 4, 3, 6, 5, 7},
		alias:          []int{0, 1, 2, 3, 4, 3, 2, 0},
		pair: [][]int{
			{0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 5, 0, 0, 5},
			{0, 0, 0, 1, 0, 0, 0, 0},
			{0, 0, 2, 0, 3, 0, 0, 0},
			{0, 6, 0, 4, 0, 0, 0, 6},
			{0, 0, 0, 0, 0, 0, 2, 0},
			{0, 0, 0, 0, 0, 1, 0, 0},
			{0, 6, 0, 0, 5, 0, 0, 0},
		},
	}
}

type vrna_param_t struct {
	stack                                     [NBPAIRS + 1][NBPAIRS + 1]int
	hairpin                                   [31]int
	bulge                                     [MAXLOOP + 1]int
	internal_loop                             [MAXLOOP + 1]int
	mismatchExt                               [NBPAIRS + 1][5][5]int
	mismatchI                                 [NBPAIRS + 1][5][5]int
	mismatch1nI                               [NBPAIRS + 1][5][5]int
	mismatch23I                               [NBPAIRS + 1][5][5]int
	mismatchH                                 [NBPAIRS + 1][5][5]int
	mismatchM                                 [NBPAIRS + 1][5][5]int
	dangle5                                   [NBPAIRS + 1][5]int
	dangle3                                   [NBPAIRS + 1][5]int
	int11                                     [NBPAIRS + 1][NBPAIRS + 1][5][5]int
	int21                                     [NBPAIRS + 1][NBPAIRS + 1][5][5][5]int
	int22                                     [NBPAIRS + 1][NBPAIRS + 1][5][5][5][5]int
	ninio                                     [5]int
	lxc                                       float64
	MLbase, MLclosing, TerminalAU, DuplexInit int
	MLintern                                  [NBPAIRS + 1]int
	Tetraloop_E                               [200]int
	Tetraloops                                [1401]rune
	Triloop_E                                 [40]int
	Triloops                                  [241]rune
	Hexaloop_E                                [40]int
	Hexaloops                                 [1801]rune
	TripleC, MultipleCA, MultipleCB           int
	gquad                                     [VRNA_GQUAD_MAX_STACK_SIZE + 1][3*VRNA_GQUAD_MAX_LINKER_LENGTH + 1]int
	gquadLayerMismatch, gquadLayerMismatchMax int
	temperature                               float64    /**<  @brief  Temperature used for loop contribution scaling */
	model_details                             *vrna_md_t /**<  @brief  Model details to be used in the recursions */
}

// 0 deg Celsius in Kelvin
var K0 float64 = 273.15

// temperature of param measurements
var Tmeasure float64 = 37 + K0

var (
	lxc37        float64 = 107.856
	ML_intern37  int     = -90
	ML_interndH  int     = -220
	ML_closing37 int     = 930
	ML_closingdH int     = 3000
	ML_BASE37    int     = 0
	ML_BASEdH    int     = 0
	MAX_NINIO    int     = 300
	ninio37      int     = 60
	niniodH      int     = 320
	TerminalAU37 int     = 50
	TerminalAUdH int     = 370
	DuplexInit37 int     = 410
	DuplexInitdH int     = 360
	TripleC37    int     = 100
	TripleCdH    int     = 1860
	MultipleCA37 int     = 30
	MultipleCAdH int     = 340
	MultipleCB37 int     = 160
	MultipleCBdH int     = 760

	GQuadAlpha37          int = -1800
	GQuadAlphadH          int = -11934
	GQuadBeta37           int = 1200
	GQuadBetadH           int = 0
	GQuadLayerMismatch37  int = 300
	GQuadLayerMismatchH   int = 0
	GQuadLayerMismatchMax int = 1

	VRNA_GQUAD_MIN_STACK_SIZE    int = 2
	VRNA_GQUAD_MAX_STACK_SIZE    int = 7
	VRNA_GQUAD_MAX_LINKER_LENGTH int = 15
	VRNA_GQUAD_MIN_LINKER_LENGTH int = 1
	VRNA_GQUAD_MIN_BOX_SIZE      int = (4 * VRNA_GQUAD_MIN_STACK_SIZE) +
		(3 * VRNA_GQUAD_MIN_LINKER_LENGTH)
	VRNA_GQUAD_MAX_BOX_SIZE int = (4 * VRNA_GQUAD_MAX_STACK_SIZE) + (3 * VRNA_GQUAD_MAX_LINKER_LENGTH)
)

var (
	/** Infinity as used in minimization routines */
	INF             int = 10000000 /* (INT_MAX/10) */
	hairpin37           = [31]int{INF, INF, INF, 540, 560, 570, 540, 600, 550, 640, 650, 660, 670, 680, 690, 690, 700, 710, 710, 720, 720, 730, 730, 740, 740, 750, 750, 750, 760, 760, 770}
	hairpindH           = [31]int{INF, INF, INF, 130, 480, 360, -290, 130, -290, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500}
	bulge37             = [31]int{INF, 380, 280, 320, 360, 400, 440, 460, 470, 480, 490, 500, 510, 520, 530, 540, 540, 550, 550, 560, 570, 570, 580, 580, 580, 590, 590, 600, 600, 600, 610}
	bulgedH             = [31]int{INF, 1060, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710, 710}
	internal_loop37     = [31]int{INF, INF, 100, 100, 110, 200, 200, 210, 230, 240, 250, 260, 270, 280, 290, 290, 300, 310, 310, 320, 330, 330, 340, 340, 350, 350, 350, 360, 360, 370, 370}
	internal_loopdH     = [31]int{INF, INF, -720, -720, -720, -680, -130, -130, -130, -130, -130, -130, -130, -130, -130, -130, -130, -130, -130, -130, -130, -130, -130, -130, -130, -130, -130, -130, -130, -130, -130}

	stack37 = [NBPAIRS + 1][NBPAIRS + 1]int{
		{INF, INF, INF, INF, INF, INF, INF, INF},
		{INF, -240, -330, -210, -140, -210, -210, -140},
		{INF, -330, -340, -250, -150, -220, -240, -150},
		{INF, -210, -250, 130, -50, -140, -130, 130},
		{INF, -140, -150, -50, 30, -60, -100, 30},
		{INF, -210, -220, -140, -60, -110, -90, -60},
		{INF, -210, -240, -130, -100, -90, -130, -90},
		{INF, -140, -150, 130, 30, -60, -90, 130},
	}
	stackdH = [NBPAIRS + 1][NBPAIRS + 1]int{
		{INF, INF, INF, INF, INF, INF, INF, INF},
		{INF, -1060, -1340, -1210, -560, -1050, -1040, -560},
		{INF, -1340, -1490, -1260, -830, -1140, -1240, -830},
		{INF, -1210, -1260, -1460, -1350, -880, -1280, -880},
		{INF, -560, -830, -1350, -930, -320, -700, -320},
		{INF, -1050, -1140, -880, -320, -940, -680, -320},
		{INF, -1040, -1240, -1280, -700, -680, -770, -680},
		{INF, -560, -830, -880, -320, -320, -680, -320},
	}

	// There are 40 Tetraloops (each 7 in length) = 280 chars (what's the extra char for? why aren't there 40 Tetraloops here?)
	Tetraloops = [281]rune{
		"CAACGG ",
		"CCAAGG ",
		"CCACGG ",
		"CCCAGG ",
		"CCGAGG ",
		"CCGCGG ",
		"CCUAGG ",
		"CCUCGG ",
		"CUAAGG ",
		"CUACGG ",
		"CUCAGG ",
		"CUCCGG ",
		"CUGCGG ",
		"CUUAGG ",
		"CUUCGG ",
		"CUUUGG "}
	Tetraloop37 = [40]int{550, 330, 370, 340, 350, 360, 370, 250, 360, 280, 370, 270, 280, 350, 370, 370}
	TetraloopdH = [40]int{690, -1030, -330, -890, -660, -750, -350, -1390, -760, -1070, -660, -1290, -1070, -620, -1530, -680}

	Triloops = [241]rune{
		"CAACG ",
		"GUUAC ",
	}
	Triloop37 = [40]int{680, 690}
	TriloopdH = [40]int{2370, 1080}

	Hexaloops = [361]rune{"ACAGUACU ",
		"ACAGUGAU ",
		"ACAGUGCU ",
		"ACAGUGUU ",
	}
	Hexaloop37 = [40]int{280, 360, 290, 180}
	HexaloopdH = [40]int{-1680, -1140, -1280, -1540}
)

func get_scaled_params() *vrna_param_t {
	var model_details *vrna_md_t = DefaultModelDetails()
	// float64
	var i, j, k, l int
	var tempf float64 = (model_details.temperature + K0) / Tmeasure

	var params *vrna_param_t = &vrna_param_t{
		model_details:         model_details,
		temperature:           model_details.temperature,
		lxc:                   lxc37 * tempf,
		TripleC:               RESCALE_dG_int(TripleC37, TripleCdH, int(tempf)),
		MultipleCA:            RESCALE_dG_int(MultipleCA37, MultipleCAdH, int(tempf)),
		MultipleCB:            RESCALE_dG_int(TerminalAU37, TerminalAUdH, int(tempf)),
		TerminalAU:            RESCALE_dG_int(TerminalAU37, TerminalAUdH, int(tempf)),
		DuplexInit:            RESCALE_dG_int(DuplexInit37, DuplexInitdH, int(tempf)),
		MLbase:                RESCALE_dG_int(ML_BASE37, ML_BASEdH, int(tempf)),
		MLclosing:             RESCALE_dG_int(ML_closing37, ML_closingdH, int(tempf)),
		gquadLayerMismatch:    RESCALE_dG_int(GQuadLayerMismatch37, GQuadLayerMismatchH, int(tempf)),
		gquadLayerMismatchMax: GQuadLayerMismatchMax,
	}

	params.ninio[2] = RESCALE_dG_int(ninio37, niniodH, int(tempf))

	for i = VRNA_GQUAD_MIN_STACK_SIZE; i <= VRNA_GQUAD_MAX_STACK_SIZE; i++ {
		for j = 3 * VRNA_GQUAD_MIN_LINKER_LENGTH; j <= 3*VRNA_GQUAD_MAX_LINKER_LENGTH; j++ {
			var GQuadAlpha_T float64 = float64(RESCALE_dG_int(GQuadAlpha37, GQuadAlphadH, int(tempf)))
			var GQuadBeta_T float64 = float64(RESCALE_dG_int(GQuadBeta37, GQuadBetadH, int(tempf)))
			params.gquad[i][j] = int(GQuadAlpha_T)*(i-1) + int(float64(GQuadBeta_T)*math.Log(float64(j)-2.0))
		}
	}

	for i = 0; i < 31; i++ {
		params.hairpin[i] = RESCALE_dG_int(hairpin37[i], hairpindH[i], int(tempf))
	}

	for i = 0; i <= Min(30, MAXLOOP); i++ {
		params.bulge[i] = RESCALE_dG_int(bulge37[i], bulgedH[i], int(tempf))
		params.internal_loop[i] = RESCALE_dG_int(internal_loop37[i], internal_loopdH[i], int(tempf))
	}

	for ; i <= MAXLOOP; i++ {
		params.bulge[i] = params.bulge[30] +
			int(params.lxc*math.Log(float64(i)/30.0))
		params.internal_loop[i] = params.internal_loop[30] +
			int(params.lxc*math.Log(float64(i)/30.0))
	}

	// This could be 281 or only the amount of chars
	for i = 0; (i * 7) < len(Tetraloops); i++ {
		params.Tetraloop_E[i] = RESCALE_dG_int(Tetraloop37[i], TetraloopdH[i], int(tempf))
	}

	for i = 0; (i * 5) < len(Triloops); i++ {
		params.Triloop_E[i] = RESCALE_dG_int(Triloop37[i], TriloopdH[i], int(tempf))
	}

	for i = 0; (i * 9) < len(Hexaloops); i++ {
		params.Hexaloop_E[i] = RESCALE_dG_int(Hexaloop37[i], HexaloopdH[i], int(tempf))
	}

	for i = 0; i <= NBPAIRS; i++ {
		params.MLintern[i] = RESCALE_dG_int(ML_intern37, ML_interndH, int(tempf))
	}

	/* stacks    G(T) = H - [H - G(T0)]*T/T0 */
	for i = 0; i <= NBPAIRS; i++ {
		for j = 0; j <= NBPAIRS; j++ {
			params.stack[i][j] = RESCALE_dG_int(stack37[i][j],
				stackdH[i][j],
				int(tempf))
		}
	}

	/* mismatches */
	for i = 0; i <= NBPAIRS; i++ {
		for j = 0; j < 5; j++ {
			for k = 0; k < 5; k++ {
				var mm int
				params.mismatchI[i][j][k] = RESCALE_dG_int(mismatchI37[i][j][k],
					mismatchIdH[i][j][k],
					int(tempf))
				params.mismatchH[i][j][k] = RESCALE_dG_int(mismatchH37[i][j][k],
					mismatchHdH[i][j][k],
					int(tempf))
				params.mismatch1nI[i][j][k] = RESCALE_dG_int(mismatch1nI37[i][j][k],
					mismatch1nIdH[i][j][k],
					int(tempf))
				params.mismatch23I[i][j][k] = RESCALE_dG_int(mismatch23I37[i][j][k],
					mismatch23IdH[i][j][k],
					int(tempf))
				if model_details.dangles > 0 {
					mm = RESCALE_dG_int(mismatchM37[i][j][k],
						mismatchMdH[i][j][k],
						int(tempf))
					if mm > 0 {
						params.mismatchM[i][j][k] = 0
					} else {
						params.mismatchM[i][j][k] = mm
					}

					mm = RESCALE_dG_int(mismatchExt37[i][j][k],
						mismatchExtdH[i][j][k],
						int(tempf))
					if mm > 0 {
						params.mismatchExt[i][j][k] = 0
					} else {
						params.mismatchExt[i][j][k] = mm
					}
				} else {
					params.mismatchExt[i][j][k] = 0
					params.mismatchM[i][j][k] = 0
				}
			}
		}
	}

	/* dangles */
	for i = 0; i <= NBPAIRS; i++ {
		for j = 0; j < 5; j++ {
			var dd int
			dd = RESCALE_dG_int(dangle5_37[i][j],
				dangle5_dH[i][j],
				int(tempf))
			if dd > 0 {
				params.dangle5[i][j] = 0
			} else {
				params.dangle5[i][j] = dd
			}

			dd = RESCALE_dG_int(dangle3_37[i][j],
				dangle3_dH[i][j],
				int(tempf))
			if dd > 0 {
				params.dangle3[i][j] = 0
			} else {
				params.dangle3[i][j] = dd
			}
		}
	}

	/* interior 1x1 loops */
	for i = 0; i <= NBPAIRS; i++ {
		for j = 0; j <= NBPAIRS; j++ {
			for k = 0; k < 5; k++ {
				for l = 0; l < 5; l++ {
					params.int11[i][j][k][l] = RESCALE_dG_int(int11_37[i][j][k][l],
						int11_dH[i][j][k][l],
						int(tempf))
				}
			}
		}
	}

	/* interior 2x1 loops */
	for i = 0; i <= NBPAIRS; i++ {
		for j = 0; j <= NBPAIRS; j++ {
			for k = 0; k < 5; k++ {
				for l = 0; l < 5; l++ {
					var m int
					for m = 0; m < 5; m++ {
						params.int21[i][j][k][l][m] = RESCALE_dG_int(int21_37[i][j][k][l][m],
							int21_dH[i][j][k][l][m],
							int(tempf))
					}
				}
			}
		}
	}

	/* interior 2x2 loops */
	for i = 0; i <= NBPAIRS; i++ {
		for j = 0; j <= NBPAIRS; j++ {
			for k = 0; k < 5; k++ {
				for l = 0; l < 5; l++ {
					var m, n int
					for m = 0; m < 5; m++ {
						for n = 0; n < 5; n++ {
							params.int22[i][j][k][l][m][n] = RESCALE_dG_int(int22_37[i][j][k][l][m][n],
								int22_dH[i][j][k][l][m][n],
								int(tempf))
						}
					}
				}
			}
		}
	}

	params.Tetraloops = Tetraloops
	params.Triloops = Triloops
	params.Hexaloops = Hexaloops

	return params
}

func vrna_eval_structure_cstr(fc *vrna_fold_compound_t, structure string) (float64, error) {
	var pt []int
	pt, err = vrna_ptable_from_string(structure)
	if err != nil {
		return nil, err
	}

	var en float64
	en, err = wrap_eval_structure(fc, structure, pt)
	if err != nil {
		return nil, err
	}

	return en, nil
}

func vrna_ptable_from_string(str string) ([]int, error) {
	var pairs [3]rune
	var pt []int
	var i, n uint64

	n = len(str)

	if n > SHRT_MAX {
		err := fmt.Sprintf("vrna_ptable_from_string: Structure too long to be converted to pair table (n=%v, max=%v)", n, SHRT_MAX)
		return nil, errors.New(err)
	}

	pt = make([]int, n+2)
	pt[0] = int(n)

	pt, err = extract_pairs(pt, str, []rune("()"))
	if err != nil {
		return nil, err
	}

	return pt, nil
}

/* requires that pt[0] already contains the length of the string! */
func extract_pairs(pt []int, structure string, pair []rune) ([]int, error) {
	var ptr []rune = []rune(structure)
	var open, close rune
	var stack []int
	var i, j, n uint64
	var hx, ptr_idx int

	n = uint64(pt[0])
	stack = make([]int, n+1)

	open = pair[0]
	close = pair[1]

	for hx, i, ptr_idx := 0, 1, 0; i <= n && ptr_idx < len(ptr); ptr_idx++ {
		if ptr[ptr_idx] == open {
			hx++
			stack[hx] = i
		} else if ptr[ptr_idx] == close {
			hx--
			j = stack[hx]

			if hx < 0 {
				// vrna_message_warning("%s\nunbalanced brackets '%2s' found while extracting base pairs",
				//                      );
				return nil, errors.New(fmt.Sprintf("%v\nunbalanced brackets '%v' found while extracting base pairs", structure,
					pair))
				// free(stack);
				// return 0;
			}

			pt[i] = j
			pt[j] = i
		}
		i++
	}

	// free(stack);

	if hx != 0 {
		return nil, errors.New(fmt.Sprintf("%v\nunbalanced brackets '%v' found while extracting base pairs", structure,
			pair))
		// return 0;
	}

	return pt, nil
	// return 1; /* success */
}

func wrap_eval_structure(fc *vrna_fold_compound_t, structure string, pt []int) (float64, error) {
	var res, gq, L int
	var l [3]int
	var energy float64

	energy = float64(INF) / 100.0

	gq = fc.params.model_details.gquad
	fc.params.model_details.gquad = 0

	// can add support for circular strands with this
	// if (vc.params.model_details.circ)
	//   res = eval_circ_pt(vc, pt, output_stream, verbosity);
	// else
	res = eval_pt(vc, pt)
	fc.params.model_details.gquad = gq

	// if gq == 1 && (parse_gquad(structure, &L, l) > 0) {
	// 	if (verbosity > 0)
	// 		vrna_cstr_print_eval_sd_corr(output_stream);

	// 	res += en_corr_of_loop_gquad(vc, 1, vc.length, structure, pt, output_stream, verbosity);
	// }
	energy = float64(res) / 100.0
	return energy
}

func eval_pt(fc *vrna_fold_compound_t, pt []int) int {
	var sn []uint64
	var i, length, energy int

	length = fc.length
	sn = fc.strand_number

	// vrna_sc_prepare(fc)
	// if fc.params.model_details.backtrack_type == 'M' {
	// 	energy = energy_of_ml_pt(vc, 0, pt)
	// } else {

	// }
	energy = energy_of_extLoop_pt(fc, 0, pt)

	for i = 1; i <= length; i++ {
		if pt[i] == 0 {
			continue
		}

		energy += stack_energy(vc, i, pt)
		i = pt[i]
	}
	for i = 1; sn[i] != sn[length]; i++ {
		if sn[i] != sn[pt[i]] {
			energy += fc.params.DuplexInit
			break
		}
	}

	return energy
}

// func vrna_sc_prepare(fc *vrna_fold_compound_t) {
// 	// prepare_sc_up_mfe(fc)
// 	// prepare_sc_bp_mfe(fc)
// }

/* populate sc.energy_up arrays for usage in MFE computations */
// func prepare_sc_up_mfe(fc *vrna_fold_compound_t) {
//   var  i, n uint64
// 	var sc *vrna_sc_t
//   n = fc.length
// 	/* prepare sc for unpaired nucleotides only if we actually have some to apply */
// 	// if (sc.state & STATE_DIRTY_UP_MFE) {
// 	/*  allocate memory such that we can access the soft constraint
// 	*  energies of a subsequence of length j starting at position i
// 	*  via sc.energy_up[i][j]
// 	*/
// 	sc = fc.sc
// 	sc_energy_up = make([][]int, n + 2)
// 	// (int **)vrna_realloc(sc.energy_up, sizeof(int *) * (n + 2));

// 	for i := 0; i < n; i++ {
// 		sc_energy_up[i] = make([]int, n - i + 2)
// 	}
// 	// sc.energy_up[i] = (int *)vrna_realloc(sc.energy_up[i], sizeof(int) * (n - i + 2));

// 	sc_energy_up[0] = make([]int, 1)
// 	sc_energy_up[n + 1] = make([]int, 1)
// 	// sc.energy_up[0]     = (int *)vrna_realloc(sc.energy_up[0], sizeof(int));
// 	// sc.energy_up[n + 1] = (int *)vrna_realloc(sc.energy_up[n + 1], sizeof(int));

// 	/* now add soft constraints as stored in container for unpaired sc */
// 	for i := 1; i <= n; i++ {
// 		populate_sc_up_mfe(fc, i, (n - i + 1))
// 	}

// 	sc_energy_up[0][0]     = 0;
// 	sc_energy_up[n + 1][0] = 0;

// 		// sc.state &= ~STATE_DIRTY_UP_MFE;
// 	// }
// }

/* pupulate sc.energy_up array at position i from sc.up_storage data */
func populate_sc_up_mfe(i, n uint64) {
	sc_energy_up[i][0] = 0
	for j := 1; j <= n; j++ {
		sc_energy_up[i][j] = sc_energy_up[i][j-1] + sc_up_storage[i+j-1]
	}
}

/* populate sc.energy_bp arrays for usage in MFE computations */
// func prepare_sc_bp_mfe(sequence string) {
//   var i, n uint64
//   n = len(sequence)

// 	/* prepare sc for base paired positions only if we actually have some to apply */
// 	// if (sc.bp_storage) {
// 	// if (sc.state & STATE_DIRTY_BP_MFE) {

// 	sc_energy_bp = make([]int, (((n + 1) * (n + 2)) / 2))
// 		// (int *)vrna_realloc(sc.energy_bp, sizeof(int) * (((n + 1) * (n + 2)) / 2));

// 	for i = 1; i < n; i++ {
// 		populate_sc_bp_mfe(fc, i, n)
// 	}

// 		// sc.state &= ~STATE_DIRTY_BP_MFE;
// 	// }
// 	// } else {
// 	// 	free_sc_bp(sc);
// 	// }
// }

var jindx []int /**<  @brief  DP matrix accessor  */

type vrna_sc_bp_storage_t struct {
	interval_start, interval_end uint64
	e                            int
}

var sc_bp_storage [][]vrna_sc_bp_storage_t /**<  @brief  Storage container for energy contributions per base pair */

func populate_sc_bp_mfe(sequence string, i, maxdist uint64) {
	var j, k, turn, n uint64
	var e int
	var idx []int

	n = len(sequence)
	turn = 3 // could be 0?
	idx = jindx

	if sc_bp_storage[i] {
		// for (k = turn + 1; k < maxdist; k++) {
		// 	j = i + k;

		// 	if (j > n)
		// 		break;

		// 	e = get_stored_bp_contributions(sc.bp_storage[i], j);

		// 	switch (sc.type) {
		// 		case VRNA_SC_DEFAULT:
		// 			sc.energy_bp[idx[j] + i] = e;
		// 			break;

		// 		case VRNA_SC_WINDOW:
		// 			sc.energy_bp_local[i][j - i] = e;
		// 			break;
		// 		}
		// }
	} else {
		for k = turn + 1; k < maxdist; k++ {
			j = i + k
			if j > n {
				break
			}

			sc_energy_bp[idx[j]+i] = 0
		}
	}
}

var (
	mismatchI37 = [NBPAIRS + 1][5][5]int{
		{
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
		},
		{
			{0, 0, 0, 0, 0},
			{0, 0, 0, -80, 0},
			{0, 0, 0, 0, 0},
			{0, -100, 0, -100, 0},
			{0, 0, 0, 0, -60},
		},
		{
			{0, 0, 0, 0, 0},
			{0, 0, 0, -80, 0},
			{0, 0, 0, 0, 0},
			{0, -100, 0, -100, 0},
			{0, 0, 0, 0, -60},
		},
		{
			{70, 70, 70, 70, 70},
			{70, 70, 70, -10, 70},
			{70, 70, 70, 70, 70},
			{70, -30, 70, -30, 70},
			{70, 70, 70, 70, 10},
		},
		{
			{70, 70, 70, 70, 70},
			{70, 70, 70, -10, 70},
			{70, 70, 70, 70, 70},
			{70, -30, 70, -30, 70},
			{70, 70, 70, 70, 10},
		},
		{
			{70, 70, 70, 70, 70},
			{70, 70, 70, -10, 70},
			{70, 70, 70, 70, 70},
			{70, -30, 70, -30, 70},
			{70, 70, 70, 70, 10},
		},
		{
			{70, 70, 70, 70, 70},
			{70, 70, 70, -10, 70},
			{70, 70, 70, 70, 70},
			{70, -30, 70, -30, 70},
			{70, 70, 70, 70, 10},
		},
		{
			{70, 70, 70, 70, 70},
			{70, 70, 70, -10, 70},
			{70, 70, 70, 70, 70},
			{70, -30, 70, -30, 70},
			{70, 70, 70, 70, 10},
		},
	}
	mismatchIdH = [NBPAIRS + 1][5][5]int{
		{
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
		},
		{
			{280, 0, 0, 280, 0},
			{0, 0, 0, -340, 0},
			{0, 0, 0, 0, 0},
			{280, -760, 0, 280, 0},
			{0, 0, 0, 0, -580},
		},
		{
			{280, 0, 0, 280, 0},
			{0, 0, 0, -340, 0},
			{0, 0, 0, 0, 0},
			{280, -760, 0, 280, 0},
			{0, 0, 0, 0, -580},
		},
		{
			{790, 500, 500, 790, 500},
			{500, 500, 500, 170, 500},
			{500, 500, 500, 500, 500},
			{790, -260, 500, 790, 500},
			{500, 500, 500, 500, -80},
		},
		{
			{790, 500, 500, 790, 500},
			{500, 500, 500, 170, 500},
			{500, 500, 500, 500, 500},
			{790, -260, 500, 790, 500},
			{500, 500, 500, 500, -80},
		},
		{
			{790, 500, 500, 790, 500},
			{500, 500, 500, 170, 500},
			{500, 500, 500, 500, 500},
			{790, -260, 500, 790, 500},
			{500, 500, 500, 500, -80},
		},
		{
			{790, 500, 500, 790, 500},
			{500, 500, 500, 170, 500},
			{500, 500, 500, 500, 500},
			{790, -260, 500, 790, 500},
			{500, 500, 500, 500, -80},
		},
		{
			{790, 500, 500, 790, 500},
			{500, 500, 500, 170, 500},
			{500, 500, 500, 500, 500},
			{790, -260, 500, 790, 500},
			{500, 500, 500, 500, -80},
		},
	}

	mismatchH37 = [NBPAIRS + 1][5][5]int{
		{
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
		},
		{
			{-80, -100, -110, -100, -80},
			{-140, -150, -150, -140, -150},
			{-80, -100, -110, -100, -80},
			{-150, -230, -150, -240, -150},
			{-100, -100, -140, -100, -210},
		},
		{
			{-50, -110, -70, -110, -50},
			{-110, -110, -150, -130, -150},
			{-50, -110, -70, -110, -50},
			{-150, -250, -150, -220, -150},
			{-100, -110, -100, -110, -160},
		},
		{
			{20, 20, -20, -10, -20},
			{20, 20, -50, -30, -50},
			{-10, -10, -20, -10, -20},
			{-50, -100, -50, -110, -50},
			{-10, -10, -30, -10, -100},
		},
		{
			{0, -20, -10, -20, 0},
			{-30, -50, -30, -60, -30},
			{0, -20, -10, -20, 0},
			{-30, -90, -30, -110, -30},
			{-10, -20, -10, -20, -90},
		},
		{
			{-10, -10, -20, -10, -20},
			{-30, -30, -50, -30, -50},
			{-10, -10, -20, -10, -20},
			{-50, -120, -50, -110, -50},
			{-10, -10, -30, -10, -120},
		},
		{
			{0, -20, -10, -20, 0},
			{-30, -50, -30, -50, -30},
			{0, -20, -10, -20, 0},
			{-30, -150, -30, -150, -30},
			{-10, -20, -10, -20, -90},
		},
		{
			{20, 20, -10, -10, 0},
			{20, 20, -30, -30, -30},
			{0, -10, -10, -10, 0},
			{-30, -90, -30, -110, -30},
			{-10, -10, -10, -10, -90},
		},
	}
	mismatchHdH = [NBPAIRS + 1][5][5]int{
		{
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
		},
		{
			{560, -570, 560, -560, -270},
			{-560, -910, -560, -560, -560},
			{-270, -570, -340, -570, -270},
			{560, -1400, 560, -920, -560},
			{-530, -570, -530, -570, -1440},
		},
		{
			{50, -520, 50, -560, -400},
			{-400, -520, -400, -560, -400},
			{50, -720, 50, -720, -420},
			{-400, -1290, -400, -620, -400},
			{-30, -720, -30, -720, -1080},
		},
		{
			{970, 140, 970, 140, 570},
			{570, 30, 570, 20, 570},
			{970, 140, 970, 140, 340},
			{570, -270, 570, 20, 570},
			{830, 140, 830, 140, -50},
		},
		{
			{230, 100, 230, 220, 190},
			{-110, -110, -260, -520, -260},
			{190, -60, -140, -60, 190},
			{220, 100, -260, 220, -260},
			{230, -60, 230, -60, -70},
		},
		{
			{970, 140, 970, 140, 570},
			{570, -20, 570, 20, 570},
			{970, 140, 970, 140, 340},
			{570, -520, 570, 20, 570},
			{830, 140, 830, 140, -380},
		},
		{
			{230, -30, 230, -60, 190},
			{-30, -30, -260, -520, -260},
			{190, -60, -140, -60, 190},
			{-260, -590, -260, -520, -260},
			{230, -60, 230, -60, -70},
		},
		{
			{970, 140, 970, 220, 570},
			{570, 30, 570, 20, 570},
			{970, 140, 970, 140, 340},
			{570, 100, 570, 220, 570},
			{830, 140, 830, 140, -50},
		},
	}

	/* mismatch_multi */
	mismatchM37 = [NBPAIRS + 1][5][5]int{
		{ /* NP.. */
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
		},
		{ /* CG.. */
			{-50, -110, -50, -140, -70},
			{-110, -110, -110, -160, -110},
			{-70, -150, -70, -150, -100},
			{-110, -130, -110, -140, -110},
			{-50, -150, -50, -150, -70},
		},
		{ /* GC.. */
			{-80, -140, -80, -140, -100},
			{-100, -150, -100, -140, -100},
			{-110, -150, -110, -150, -140},
			{-100, -140, -100, -160, -100},
			{-80, -150, -80, -150, -120},
		},
		{ /* GU.. */
			{-50, -80, -50, -50, -50},
			{-50, -100, -70, -50, -70},
			{-60, -80, -60, -80, -60},
			{-70, -110, -70, -80, -70},
			{-50, -80, -50, -80, -50},
		},
		{ /* UG.. */
			{-30, -30, -60, -60, -60},
			{-30, -30, -60, -60, -60},
			{-70, -100, -70, -100, -80},
			{-60, -80, -60, -80, -60},
			{-60, -100, -70, -100, -60},
		},
		{ /* AU.. */
			{-50, -80, -50, -80, -50},
			{-70, -100, -70, -110, -70},
			{-60, -80, -60, -80, -60},
			{-70, -110, -70, -120, -70},
			{-50, -80, -50, -80, -50},
		},
		{ /* UA.. */
			{-60, -80, -60, -80, -60},
			{-60, -80, -60, -80, -60},
			{-70, -100, -70, -100, -80},
			{-60, -80, -60, -80, -60},
			{-70, -100, -70, -100, -80},
		},
		{ /* NN.. */
			{-30, -30, -50, -50, -50},
			{-30, -30, -60, -50, -60},
			{-60, -80, -60, -80, -60},
			{-60, -80, -60, -80, -60},
			{-50, -80, -50, -80, -50},
		},
	}

	/* mismatch_multi_enthalpies */
	mismatchMdH = [NBPAIRS + 1][5][5]int{
		{ /* NP.. */
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
		},
		{ /* CG.. */
			{50, -400, 50, -400, -30},
			{-520, -520, -720, -710, -720},
			{50, -400, 50, -400, -30},
			{-560, -560, -720, -620, -720},
			{-400, -400, -420, -400, -500},
		},
		{ /* GC.. */
			{-270, -560, -270, -560, -530},
			{-570, -910, -570, -820, -570},
			{-340, -560, -340, -560, -530},
			{-560, -560, -570, -920, -570},
			{-270, -560, -270, -560, -860},
		},
		{ /* GU.. */
			{310, -480, -180, 310, 140},
			{310, -480, -430, 310, -430},
			{-140, -630, -510, -630, -140},
			{-150, -890, -430, -150, -430},
			{140, -630, -180, -630, 140},
		},
		{ /* UG.. */
			{600, 200, 600, 200, 460},
			{-60, -340, -230, -60, -230},
			{600, 200, 600, 200, 460},
			{-230, -350, -230, -350, -230},
			{200, 200, -30, 200, 160},
		},
		{ /* AU.. */
			{140, -400, -180, -380, 140},
			{-380, -400, -430, -380, -430},
			{-140, -630, -510, -630, -140},
			{-430, -890, -430, -890, -430},
			{140, -630, -180, -630, 140},
		},
		{ /* UA.. */
			{600, 200, 600, 200, 460},
			{-230, -390, -230, -310, -230},
			{600, 200, 600, 200, 460},
			{-230, -350, -230, -350, -230},
			{200, 200, -30, 200, -170},
		},
		{ /* NN.. */
			{600, 200, 600, 310, 460},
			{310, -340, -230, 310, -230},
			{600, 200, 600, 200, 460},
			{-150, -350, -230, -150, -230},
			{200, 200, -30, 200, 160},
		},
	}

	mismatch1nI37 = [NBPAIRS + 1][5][5]int{
		{
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
		},
		{
			{0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0},
		},
		{
			{0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0},
		},
		{
			{70, 70, 70, 70, 70},
			{70, 70, 70, 70, 70},
			{70, 70, 70, 70, 70},
			{70, 70, 70, 70, 70},
			{70, 70, 70, 70, 70},
		},
		{
			{70, 70, 70, 70, 70},
			{70, 70, 70, 70, 70},
			{70, 70, 70, 70, 70},
			{70, 70, 70, 70, 70},
			{70, 70, 70, 70, 70},
		},
		{
			{70, 70, 70, 70, 70},
			{70, 70, 70, 70, 70},
			{70, 70, 70, 70, 70},
			{70, 70, 70, 70, 70},
			{70, 70, 70, 70, 70},
		},
		{
			{70, 70, 70, 70, 70},
			{70, 70, 70, 70, 70},
			{70, 70, 70, 70, 70},
			{70, 70, 70, 70, 70},
			{70, 70, 70, 70, 70},
		},
		{
			{70, 70, 70, 70, 70},
			{70, 70, 70, 70, 70},
			{70, 70, 70, 70, 70},
			{70, 70, 70, 70, 70},
			{70, 70, 70, 70, 70},
		},
	}
	mismatch1nIdH = [NBPAIRS + 1][5][5]int{
		{
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
		},
		{
			{0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0},
		},
		{
			{0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0},
		},
		{
			{500, 500, 500, 500, 500},
			{500, 500, 500, 500, 500},
			{500, 500, 500, 500, 500},
			{500, 500, 500, 500, 500},
			{500, 500, 500, 500, 500},
		},
		{
			{500, 500, 500, 500, 500},
			{500, 500, 500, 500, 500},
			{500, 500, 500, 500, 500},
			{500, 500, 500, 500, 500},
			{500, 500, 500, 500, 500},
		},
		{
			{500, 500, 500, 500, 500},
			{500, 500, 500, 500, 500},
			{500, 500, 500, 500, 500},
			{500, 500, 500, 500, 500},
			{500, 500, 500, 500, 500},
		},
		{
			{500, 500, 500, 500, 500},
			{500, 500, 500, 500, 500},
			{500, 500, 500, 500, 500},
			{500, 500, 500, 500, 500},
			{500, 500, 500, 500, 500},
		},
		{
			{500, 500, 500, 500, 500},
			{500, 500, 500, 500, 500},
			{500, 500, 500, 500, 500},
			{500, 500, 500, 500, 500},
			{500, 500, 500, 500, 500},
		},
	}

	mismatch23I37 = [NBPAIRS + 1][5][5]int{
		{
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
		},
		{
			{0, 0, 0, 0, 0},
			{0, 0, 0, -50, 0},
			{0, 0, 0, 0, 0},
			{0, -110, 0, -70, 0},
			{0, 0, 0, 0, -30},
		},
		{
			{0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0},
			{0, -120, 0, -70, 0},
			{0, 0, 0, 0, -30},
		},
		{
			{70, 70, 70, 70, 70},
			{70, 70, 70, 70, 70},
			{70, 70, 70, 70, 70},
			{70, -40, 70, 0, 70},
			{70, 70, 70, 70, 40},
		},
		{
			{70, 70, 70, 70, 70},
			{70, 70, 70, 20, 70},
			{70, 70, 70, 70, 70},
			{70, -40, 70, 0, 70},
			{70, 70, 70, 70, 40},
		},
		{
			{70, 70, 70, 70, 70},
			{70, 70, 70, 70, 70},
			{70, 70, 70, 70, 70},
			{70, -40, 70, 0, 70},
			{70, 70, 70, 70, 40},
		},
		{
			{70, 70, 70, 70, 70},
			{70, 70, 70, 20, 70},
			{70, 70, 70, 70, 70},
			{70, -40, 70, 0, 70},
			{70, 70, 70, 70, 40},
		},
		{
			{70, 70, 70, 70, 70},
			{70, 70, 70, 70, 70},
			{70, 70, 70, 70, 70},
			{70, -40, 70, 0, 70},
			{70, 70, 70, 70, 40},
		},
	}
	mismatch23IdH = [NBPAIRS + 1][5][5]int{
		{
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
		},
		{
			{0, 0, 0, 0, 0},
			{0, 0, 0, -570, 0},
			{0, 0, 0, 0, 0},
			{0, -860, 0, -900, 0},
			{0, 0, 0, 0, -640},
		},
		{
			{0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0},
			{0, 0, 0, 0, 0},
			{0, -1090, 0, -900, 0},
			{0, 0, 0, 0, -640},
		},
		{
			{500, 500, 500, 500, 500},
			{500, 500, 500, 500, 500},
			{500, 500, 500, 500, 500},
			{500, -580, 500, -400, 500},
			{500, 500, 500, 500, -140},
		},
		{
			{500, 500, 500, 500, 500},
			{500, 500, 500, -60, 500},
			{500, 500, 500, 500, 500},
			{500, -360, 500, -400, 500},
			{500, 500, 500, 500, -140},
		},
		{
			{500, 500, 500, 500, 500},
			{500, 500, 500, 500, 500},
			{500, 500, 500, 500, 500},
			{500, -580, 500, -400, 500},
			{500, 500, 500, 500, -140},
		},
		{
			{500, 500, 500, 500, 500},
			{500, 500, 500, -60, 500},
			{500, 500, 500, 500, 500},
			{500, -360, 500, -400, 500},
			{500, 500, 500, 500, -140},
		},
		{
			{500, 500, 500, 500, 500},
			{500, 500, 500, 500, 500},
			{500, 500, 500, 500, 500},
			{500, -360, 500, -400, 500},
			{500, 500, 500, 500, -140},
		},
	}

	/* mismatch_exterior */
	mismatchExt37 = [NBPAIRS + 1][5][5]int{
		{ /* NP.. */
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
		},
		{ /* CG.. */
			{-50, -110, -50, -140, -70},
			{-110, -110, -110, -160, -110},
			{-70, -150, -70, -150, -100},
			{-110, -130, -110, -140, -110},
			{-50, -150, -50, -150, -70},
		},
		{ /* GC.. */
			{-80, -140, -80, -140, -100},
			{-100, -150, -100, -140, -100},
			{-110, -150, -110, -150, -140},
			{-100, -140, -100, -160, -100},
			{-80, -150, -80, -150, -120},
		},
		{ /* GU.. */
			{-50, -80, -50, -50, -50},
			{-50, -100, -70, -50, -70},
			{-60, -80, -60, -80, -60},
			{-70, -110, -70, -80, -70},
			{-50, -80, -50, -80, -50},
		},
		{ /* UG.. */
			{-30, -30, -60, -60, -60},
			{-30, -30, -60, -60, -60},
			{-70, -100, -70, -100, -80},
			{-60, -80, -60, -80, -60},
			{-60, -100, -70, -100, -60},
		},
		{ /* AU.. */
			{-50, -80, -50, -80, -50},
			{-70, -100, -70, -110, -70},
			{-60, -80, -60, -80, -60},
			{-70, -110, -70, -120, -70},
			{-50, -80, -50, -80, -50},
		},
		{ /* UA.. */
			{-60, -80, -60, -80, -60},
			{-60, -80, -60, -80, -60},
			{-70, -100, -70, -100, -80},
			{-60, -80, -60, -80, -60},
			{-70, -100, -70, -100, -80},
		},
		{ /* NN.. */
			{-30, -30, -50, -50, -50},
			{-30, -30, -60, -50, -60},
			{-60, -80, -60, -80, -60},
			{-60, -80, -60, -80, -60},
			{-50, -80, -50, -80, -50},
		},
	}

	/* mismatch_exterior_enthalpies */
	mismatchExtdH = [NBPAIRS + 1][5][5]int{
		{ /* NP.. */
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
			{INF, INF, INF, INF, INF},
		},
		{ /* CG.. */
			{50, -400, 50, -400, -30},
			{-520, -520, -720, -710, -720},
			{50, -400, 50, -400, -30},
			{-560, -560, -720, -620, -720},
			{-400, -400, -420, -400, -500},
		},
		{ /* GC.. */
			{-270, -560, -270, -560, -530},
			{-570, -910, -570, -820, -570},
			{-340, -560, -340, -560, -530},
			{-560, -560, -570, -920, -570},
			{-270, -560, -270, -560, -860},
		},
		{ /* GU.. */
			{310, -480, -180, 310, 140},
			{310, -480, -430, 310, -430},
			{-140, -630, -510, -630, -140},
			{-150, -890, -430, -150, -430},
			{140, -630, -180, -630, 140},
		},
		{ /* UG.. */
			{600, 200, 600, 200, 460},
			{-60, -340, -230, -60, -230},
			{600, 200, 600, 200, 460},
			{-230, -350, -230, -350, -230},
			{200, 200, -30, 200, 160},
		},
		{ /* AU.. */
			{140, -400, -180, -380, 140},
			{-380, -400, -430, -380, -430},
			{-140, -630, -510, -630, -140},
			{-430, -890, -430, -890, -430},
			{140, -630, -180, -630, 140},
		},
		{ /* UA.. */
			{600, 200, 600, 200, 460},
			{-230, -390, -230, -310, -230},
			{600, 200, 600, 200, 460},
			{-230, -350, -230, -350, -230},
			{200, 200, -30, 200, -170},
		},
		{ /* NN.. */
			{600, 200, 600, 310, 460},
			{310, -340, -230, 310, -230},
			{600, 200, 600, 200, 460},
			{-150, -350, -230, -150, -230},
			{200, 200, -30, 200, 160},
		},
	}

	/* dangle5 */
	dangle5_37 = [NBPAIRS + 1][5]int{
		/*           N      A      C      G      U */
		/* NP */ {INF, INF, INF, INF, INF},
		/* CG */ {-10, -50, -30, -20, -10},
		/* GC */ {-0, -20, -30, -0, -0},
		/* GU */ {-20, -30, -30, -40, -20},
		/* UG */ {-10, -30, -10, -20, -20},
		/* AU */ {-20, -30, -30, -40, -20},
		/* UA */ {-10, -30, -10, -20, -20},
		/* NN */ {-0, -20, -10, -0, -0},
	}

	/* dangle3 */
	dangle3_37 = [NBPAIRS + 1][5]int{
		/*           N      A      C      G      U */
		/* NP */ {INF, INF, INF, INF, INF},
		/* CG */ {-40, -110, -40, -130, -60},
		/* GC */ {-80, -170, -80, -170, -120},
		/* GU */ {-10, -70, -10, -70, -10},
		/* UG */ {-50, -80, -50, -80, -60},
		/* AU */ {-10, -70, -10, -70, -10},
		/* UA */ {-50, -80, -50, -80, -60},
		/* NN */ {-10, -70, -10, -70, -10},
	}

	/* dangle5_enthalpies */
	dangle5_dH = [NBPAIRS + 1][5]int{
		/*           N      A      C      G      U */
		/* NP */ {INF, INF, INF, INF, INF},
		/* CG */ {330, -240, 330, 80, -140},
		/* GC */ {70, -160, 70, -460, -40},
		/* GU */ {310, 160, 220, 70, 310},
		/* UG */ {690, -50, 690, 60, 60},
		/* AU */ {310, 160, 220, 70, 310},
		/* UA */ {690, -50, 690, 60, 60},
		/* NN */ {690, 160, 690, 80, 310},
	}

	/* dangle3_enthalpies */
	dangle3_dH = [NBPAIRS + 1][5]int{
		/*           N      A      C      G      U */
		/* NP */ {INF, INF, INF, INF, INF},
		/* CG */ {-280, -740, -280, -640, -360},
		/* GC */ {-410, -900, -410, -860, -750},
		/* GU */ {-70, -570, -70, -580, -220},
		/* UG */ {-90, -490, -90, -550, -230},
		/* AU */ {-70, -570, -70, -580, -220},
		/* UA */ {-90, -490, -90, -550, -230},
		/* NN */ {-70, -490, -70, -550, -220},
	}
)

// func vrna_params_prepare(fc *vrna_fold_compound_t) {
// // VRNA_OPTION_MFE | VRNA_OPTION_EVAL_ONLY

//     var md_p *vrna_md_t

//     /*
//      *  every vrna_fold_compound_t must have a vrna_paramt_t structure attached
//      *  to it that holds the current model details. So we just use this here as
//      *  the reference model
//      */
//     md_p = &(fc.params.model_details)

//     if (options & VRNA_OPTION_PF)
//     {
//       /* remove previous parameters if present and they differ from reference model */
//       if (fc.exp_params)
//       {
//         if (memcmp(md_p, &(fc.exp_params.model_details), sizeof(vrna_md_t)) != 0)
//         {
//           free(fc.exp_params);
//           fc.exp_params = NULL;
//         }
//       }

//       if (!fc.exp_params)
//         fc.exp_params = (fc.type == VRNA_FC_TYPE_SINGLE) ? vrna_exp_params(md_p) : vrna_exp_params_comparative(fc.n_seq, md_p);
//     }

// }

func sanitize_bp_span(fc *vrna_fold_compound_t) {
	var md *vrna_md_t
	md = fc.params.model_details

	/* non-local fold mode */
	md.window_size = fc.length

	if md.max_bp_span <= 0 || md.max_bp_span > md.window_size {
		md.max_bp_span = md.window_size
	}
}

func set_fold_compound(fc *vrna_fold_compound_t) {
	// VRNA_OPTION_MFE | VRNA_OPTION_EVAL_ONLY
	var sequence string
	var sequences []string
	var length, s int
	var md_p *vrna_md_t

	md_p = fc.params.model_details

	sequence = fc.sequence

	fc.sequence = nil

	/* split input sequences at default delimiter '&' */
	// sequences = vrna_strsplit(sequence, NULL);

	vrna_sequence_add(fc, sequence)

	if fc.strands > 1 {
		fc.cutpoint = fc.nucleotides[0].length + 1

		if md_p.min_loop_size == TURN {
			md_p.min_loop_size = 0 /* is it safe to set this here? */
		}
	}

	// if !(options & VRNA_OPTION_EVAL_ONLY) {
	// 	fc.ptype = (aux & WITH_PTYPE) ? vrna_ptypes(fc.sequence_encoding2, md_p) : NULL;
	// 	/* backward compatibility ptypes */
	// 	fc.ptype_pf_compat =
	// 			(aux & WITH_PTYPE_COMPAT) ? get_ptypes(fc.sequence_encoding2, md_p, 1) : NULL;
	// }

	vrna_sequence_prepare(fc)

	fc.iindx = vrna_idx_row_wise(fc.length)
	fc.jindx = vrna_idx_col_wise(fc.length)

}

// vivek: big potential point of failure fiven the cryptic memcpy-ing
func vrna_sequence_add(fc *vrna_fold_compound_t, sequence string) {
	// VRNA_SEQUENCE_RNA
	var add_length uint64

	add_length = len(sequence)

	/* add the sequence to the nucleotides container */
	// fc.nucleotides = (vrna_seq_t *)vrna_realloc(vc.nucleotides,
	// 																							sizeof(vrna_seq_t) *
	// 																							(vc.strands + 1));
	fc.nucleotides = make([]vrna_seq_t, 1)
	set_sequence(&(fc.nucleotides[fc.strands]),
		sequence,
		&(fc.params.model_details))

	/* increase strands counter */
	fc.strands++

	/* add new sequence to initial order of all strands */
	// fc.sequence = (char *)vrna_realloc(vc.sequence,
	// 																		sizeof(char) *
	// 																		(vc.length + add_length + 1));
	// vivek: should we add an ampersand between sequences?
	fc.sequence = vs.sequence + sequence

	/* add encoding for new strand */
	fc.sequence_encoding = make([]int, fc.length+add_length+2)

	// fc.sequence_encoding = fc.sequence_encoding + '\0' + fc.nucleotides[fc.strands - 1].encoding

	/* restore circular encoding */
	fc.sequence_encoding = fc.sequence_encoding + fc.sequence_encoding[1]
	fc.sequence_encoding[0] = fc.sequence_encoding[fc.length+add_length]

	/* add encoding2 (simple encoding) for new strand */
	fc.sequence_encoding2 = make([]int, fc.length+add_length+2)

	enc := vrna_seq_encode_simple(fc.nucleotides[fc.strands-1].sequence,
		&(fc.params.model_details))
	fc.sequence_encoding2[fc.length+1:] = enc[1 : 1+add_length]
	// memcpy(vc.sequence_encoding2 + vc.length + 1,
	// 				enc + 1,
	// 				add_length * sizeof(short))

	fc.sequence_encoding2[fc.length+add_length+1 : 0] = fc.sequence_encoding2[1:1]
	fc.sequence_encoding2[0] = int(fc.length + add_length)

	/* finally, increase length property of the fold compound */
	fc.length = fc.length + add_length
}

func set_sequence(obj *vrna_seq_t, sequence string, md *vrna_md_t) {
	obj.str = strings.ToUpper(sequence)
	obj.length = len(sequence)

	obj.encoding = vrna_seq_encode(obj.sequence, md)
	obj.encoding3 = make([]int, obj.length+1)
	obj.encoding5 = make([]int, obj.length+1)

	if md.circ {
		// for (size_t i = obj.length; i > 0; i--) {
		//   if (obj.encoding[i] == 0) /* no nucleotide, i.e. gap */
		//     continue;

		//   obj.encoding5[1] = obj.encoding[i];
		//   break;
		// }
		// for (size_t i = 1; i <= obj.length; i++) {
		//   if (obj.encoding[i] == 0) /* no nucleotide, i.e. gap */
		//     continue;

		//   obj.encoding3[obj.length] = obj.encoding[i];
		//   break;
		// }
	} else {
		obj.encoding5[1] = 0
		obj.encoding3[obj.length] = 0
	}

	var i, p uint64
	for i, p = 1, 0; i < obj.length; i++ {
		if obj.encoding[i] == 0 {
			obj.encoding5[i+1] = obj.encoding5[i]
		} else {
			obj.encoding5[i+1] = obj.encoding[i]
		}
	}

	for i := obj.length; i > 1; i-- {
		if obj.encoding[i] == 0 {
			obj.encoding3[i-1] = obj.encoding3[i]
		} else {
			obj.encoding3[i-1] = obj.encoding[i]
		}
	}
}

func vrna_seq_encode(sequence string, md *vrna_md_t) []int {
	var i, l uint64
	var S []int

	if sequence && md {
		S = vrna_seq_encode_simple(sequence, md)

		l = uint64(len(sequence))

		for i = 1; i <= l; i++ {
			S[i] = md.alias[S[i]]
		}

		S[l+1] = S[1]
		S[0] = S[l]
	}

	return S
}

// returns [len(sequence) + 2]int
func vrna_seq_encode_simple(sequence string,
	md *vrna_md_t) []int {
	var i, l uint64
	var S []int

	if sequence && md {
		l = uint64(len(sequence))
		S = make([]int, l+2)

		for i = 1; i <= l; i++ { /* make numerical encoding of sequence */
			S[i] = vrna_nucleotide_encode(sequence[i-1], md)
		}

		S[l+1] = S[1]
		S[0] = int(l)
	}

	return S
}

// vivek: I don't know if this variable name is a joke or not...
var Law_and_Order []rune = []rune("_ACGUTXKI")

func vrna_nucleotide_encode(c rune,
	md *vrna_md_t) int {
	/* return numerical representation of nucleotide used e.g. in vrna_md_t.pair[][] */
	var code int = -1

	c = unicode.ToUpper(c)

	if md {
		if md.energy_set > 0 {
			code = int(c-'A') + 1
		} else {
			pos := strings.Index(Law_and_Order, string(c))
			if pos == -1 {
				code = 0
			} else {
				code = pos
			}

			if code > 5 {
				code = 0
			}

			if code > 4 {
				code-- /* make T and U equivalent */
			}
		}
	}

	return code
}

func vrna_sequence_prepare(fc *vrna_fold_compound_t) {
	var cnt, i uint64

	fc.strand_order = nil
	fc.strand_start = nil
	fc.strand_end = nil
	fc.strand_number = make([]uint64, fc.length+2)

	/* 1. store initial strand order */
	fc.strand_order = make([]uint64, fc.strands+1)
	for cnt = 0; cnt < fc.strands; cnt++ {
		fc.strand_order[cnt] = cnt
	}

	/* 2. mark start and end positions of sequences */
	fc.strand_start = make([]uint64, fc.strands+1)
	fc.strand_end = make([]uint64, fc.strands+1)

	fc.strand_start[0] = 1
	fc.strand_end[0] = fc.strand_start[0] + fc.nucleotides[0].length - 1

	for cnt = 1; cnt < fc.strands; cnt++ {
		fc.strand_start[cnt] = fc.strand_end[cnt-1] + 1
		fc.strand_end[cnt] = fc.strand_start[cnt] + fc.nucleotides[cnt].length - 1
		for i = fc.strand_start[cnt]; i <= fc.strand_end[cnt]; i++ {
			fc.strand_number[i] = cnt
		}
	}

	/* this sets pos. n + 1 as well */
	fc.strand_number[fc.length+1] = fc.strands - 1
}

func vrna_idx_row_wise(length uint64) []int {
	var i uint64
	var idx []int = make([]int, length+1)

	for i = 1; i <= length; i++ {
		idx[i] = (((length + 1 - i) * (length - i)) / 2) + length + 1
	}
	return idx
}

func vrna_idx_col_wise(length uint64) []int {
	var i uint64
	var idx []int = make([]int, length+1)

	for i = 1; i <= length; i++ {
		idx[i] = (i * (i - 1)) / 2
	}
	return idx
}

func energy_of_extLoop_pt(fc *vrna_fold_compound_t, i int, pt []int) int {
	var sn []uint64
	var energy, mm5, mm3, bonus, p, q, q_prev, length, dangle_model, n_seq, ss, u,
		start int
	var s, s1 []int
	var S, S5, S3 [][]int
	var a2s [][]uint64
	var P *vrna_param_t
	var md *vrna_md_t
	var sc *vrna_sc_t
	var scs []*vrna_sc_t

	/* initialize vars */
	length = fc.length
	sn = fc.strand_number
	P = fc.params
	md = P.model_details
	dangle_model = md.dangles
	s = fc.sequence_encoding2
	s1 = fc.sequence_encoding
	sc = fc.sc
	S = nil
	S5 = nil
	S3 = nil
	a2s = nil
	n_seq = 1
	scs = nil

	energy = 0
	bonus = 0
	p = 1
	start = 1
	q_prev = -1

	/* seek to opening base of first stem */
	for p <= length && !pt[p] {
		p++
	}

	// vivek: do we ever reach here? sc is not init-ed anywhere. Does C init it automatically?
	/* add soft constraints for first unpaired nucleotides */
	if sc {
		if sc.energy_up {
			bonus += sc.energy_up[start][p-start]
		}
		/* how do we handle generalized soft constraints here ? */
	}

	for p < length {
		var tt int
		/* p must have a pairing partner */
		q = int(pt[p])

		/* get type of base pair (p,q) */
		tt = vrna_get_ptype_md(s[p], s[q], md)
		if (sn[p-1] == sn[p]) && (p > 1) {
			mm5 = s1[p-1]
		} else {
			mm5 = -1
		}

		if (sn[q] == sn[q+1]) && (q < length) {
			mm3 = s1[q+1]
		} else {
			mm3 = -1
		}
		energy += vrna_E_ext_stem(tt, mm5, mm3, P)

		/* seek to the next stem */
		p = q + 1
		q_prev = q
		for p <= length && !pt[p] {
			p++
		}

		/* add soft constraints for unpaired region */
		if sc && (q_prev+1 <= length) {
			if sc.energy_up {
				bonus += sc.energy_up[q_prev+1][p-q_prev-1]
			}
			/* how do we handle generalized soft constraints here ? */
		}

		if p == i {
			break /* cut was in loop */
		}
	}

	return energy + bonus
}

func vrna_E_ext_stem(type_1 uint64, n5d int, n3d int, p *vrna_param_t) int {
	var energy int = 0

	if n5d >= 0 && n3d >= 0 {
		energy += p.mismatchExt[type_1][n5d][n3d]
	} else if n5d >= 0 {
		energy += p.dangle5[type_1][n5d]
	} else if n3d >= 0 {
		energy += p.dangle3[type_1][n3d]
	}

	if type_1 > 2 {
		energy += p.TerminalAU
	}

	return energy
}

func stack_energy(vrna_fold_compound_t *fc, i int, pt []int) int {
	/* recursively calculate energy of substructure enclosed by (i,j) */
	var sn, so, ss []uint64
	var ee, energy, j, p, q int
	var sequence string
	var s []int
	var P *vrna_param_t
	var md *vrna_md_t

	sn = fc.strand_number
	so = fc.strand_order
	ss = fc.strand_start
	s = fc.sequence_encoding2
	P = fc.params
	md = &(P.model_details)
	energy = 0

	j = pt[i]

	sequence = fc.sequence
	if md.pair[s[i]][s[j]] == 0 {
		panic(fmt.Sprintf("bases %v and %v (%v%v) can't pair!",
			i, j,
			sequence[i-1],
			sequence[j-1]))
	}

	p = i
	q = j

	for p < q {
		/* process all stacks and interior loops */
		p++
		q--

		for pt[p] == 0 {
			p++
		}

		for pt[q] == 0 {
			q--
		}

		if (pt[q] != int(p)) || (p > q) {
			break
		}

		ee = 0

		if md.pair[s[q]][s[p]] == 0 {
			// if (verbosity_level > VRNA_VERBOSITY_QUIET) {
			panic(fmt.Sprintf("bases %d and %d (%c%c) can't pair!",
				p, q,
				sequence[p-1],
				sequence[q-1]))
			// }
		}

		ee = vrna_eval_int_loop(fc, i, j, p, q)

		// if (verbosity_level > 0) {
		//   vrna_cstr_print_eval_int_loop(output_stream,
		//                                 i, j,
		//                                 string[i - 1], string[j - 1],
		//                                 p, q,
		//                                 string[p - 1], string[q - 1],
		//                                 (vc.type == VRNA_FC_TYPE_COMPARATIVE) ?
		//                                 (int)ee / (int)vc.n_seq :
		//                                 ee);
		// }

		energy += ee
		i = p
		j = q
	} /* end while */

	/* p,q don't pair must have found hairpin or multiloop */

	if p > q {
		/* hairpin */
		ee = vrna_eval_hp_loop(fc, i, j)
		energy += ee

		// if (verbosity_level > 0) {
		//   vrna_cstr_print_eval_hp_loop(output_stream,
		//                                i, j,
		//                                string[i - 1], string[j - 1],
		//                                (vc.type == VRNA_FC_TYPE_COMPARATIVE) ?
		//                                (int)ee / (int)vc.n_seq :
		//                                ee);
		// }

		return energy
	}

	/* (i,j) is exterior pair of multiloop */
	for p < j {
		/* add up the contributions of the substructures of the ML */
		energy += stack_energy(vc, p, pt)
		p = pt[p]
		/* search for next base pair in multiloop */
		p++
		for pt[p] == 0 {
			p++
		}
	}

	ee = 0

	var ii int = cut_in_loop(i, pt, sn)
	if ii == 0 {
		ee, err = energy_of_ml_pt(vc, i, pt)
		if err != nil {
			panic(err)
		}
	} else {
		ee = energy_of_extLoop_pt(vc, ii, pt)
	}

	energy += ee
	// if (verbosity_level > 0) {
	//   vrna_cstr_print_eval_mb_loop(output_stream,
	//                                i, j,
	//                                string[i - 1], string[j - 1],
	//                                (vc.type == VRNA_FC_TYPE_COMPARATIVE) ?
	//                                (int)ee / (int)vc.n_seq :
	//                                ee);
	// }

	return energy
}

func vrna_eval_int_loop(fc *vrna_fold_compound_t, i, j, k, l int) {
	e := INF

	e = eval_int_loop(fc, i, j, k, l)

	return e
}

func eval_int_loop(fc *vrna_fold_compound_t, i, j, k, l int) int {
	var n_seq, s uint64
	var sn, ss []uint64
	var a2s [][]uint64
	var e, type_1, type2, with_ud int
	var rtype []int
	var S, S2 []int
	var SS, S5, S3 [][]int
	var P *vrna_param_t
	var md *vrna_md_t
	// var domains_up *vrna_ud_t
	// struct sc_int_dat sc_wrapper;

	n_seq = 1
	P = fc.params
	md = &(P.model_details)
	sn = fc.strand_number
	ss = fc.strand_start
	rtype = &(md.rtype[0])
	S = fc.sequence_encoding
	S2 = fc.sequence_encoding2
	SS = nil
	S5 = nil
	S3 = nil
	a2s = nil
	// domains_up  = fc.domains_up;
	// with_ud     = ((domains_up) && (domains_up.energy_cb)) ? 1 : 0;
	e = INF

	// init_sc_int(fc, &sc_wrapper);

	{
		var energy, e5, e3, u1, u2 int

		energy = 0

		type_1 = vrna_get_ptype_md(S2[i], S2[j], md)
		type2 = vrna_get_ptype_md(S2[l], S2[k], md)

		u1 = k - i - 1
		u2 = j - l - 1

		if (sn[i] == sn[k]) && (sn[l] == sn[j]) {
			/* regular interior loop */
			energy = E_IntLoop(u1, u2, type_1, type2, S[i+1], S[j-1], S[k-1], S[l+1], P)
		} else {
			/* interior loop like cofold structure */
			var Si, Sj int
			if sn[i+1] == sn[i] {
				Si = S[i+1]
			} else {
				Si = -1
			}

			if sn[j] == sn[j-1] {
				Sj = S[j-1]
			} else {
				Sj = -1
			}

			energy = E_IntLoop_Co(rtype[type_1], rtype[type2],
				i, j, k, l,
				ss[fc.strand_order[1]], /* serves as cut point substitute */
				Si, Sj,
				S[k-1], S[l+1],
				md.dangles,
				P)
		}

		/* add soft constraints */
		// vivek: do something about sc
		// if (sc_wrapper.pair)
		//   energy += sc_wrapper.pair(i, j, k, l, &sc_wrapper);

		e = energy

		// vivek: do something about sc
		//   if with_ud {
		//     e5, e3 = 0, 0

		//     u1  = k - i - 1
		//     u2  = j - l - 1

		//     if u1 > 0 {
		//       e5 = domains_up.energy_cb(fc,
		//                                  i + 1, k - 1,
		//                                  VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
		//                                  domains_up.data);
		//     }

		//     if (u2 > 0) {
		//       e3 = domains_up.energy_cb(fc,
		//                                  l + 1, j - 1,
		//                                  VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
		//                                  domains_up.data);
		//     }

		//     e = Min(e, energy + e5);
		//     e = Min(e, energy + e3);
		//     e = Min(e, energy + e5 + e3);
		//   }
		// }

		// free_sc_int(&sc_wrapper);
	}

	return e
}

func vrna_get_ptype_md(i, j int, md *vrna_md_t) uint64 {
	var tt uint64 = uint64(md.pair[i][j])

	if tt == 0 {
		return 7
	} else {
		return tt
	}
}

func E_IntLoop(n1, n2, type_1, type_2, si1, sj1, sp1, sq1 int, P *vrna_param_t) int {
	/* compute energy of degree 2 loop (stack bulge or interior) */
	var nl, ns, u, energy int

	if n1 > n2 {
		nl = n1
		ns = n2
	} else {
		nl = n2
		ns = n1
	}

	if nl == 0 {
		return P.stack[type_1][type_2] /* stack */
	}

	if ns == 0 {
		/* bulge */
		if nl <= MAXLOOP {
			energy = P.bulge[nl]
		} else {
			energy = P.bulge[30] + int(P.lxc*Math.Log(float64(nl)/30.0))
		}

		if nl == 1 {
			energy += P.stack[type_1][type_2]
		} else {
			if type_1 > 2 {
				energy += P.TerminalAU
			}

			if type_2 > 2 {
				energy += P.TerminalAU
			}
		}

		return energy
	} else {
		/* interior loop */
		if ns == 1 {
			if nl == 1 {
				/* 1x1 loop */
				return P.int11[type_1][type_2][si1][sj1]
			}

			if nl == 2 {
				/* 2x1 loop */
				if n1 == 1 {
					energy = P.int21[type_1][type_2][si1][sq1][sj1]
				} else {
					energy = P.int21[type_2][type_1][sq1][si1][sp1]
				}
				return energy
			} else {
				/* 1xn loop */
				if nl+1 <= MAXLOOP {
					energy = P.internal_loop[nl+1]
				} else {
					energy = P.internal_loop[30] + int(P.lxc*math.Log((float64(nl)+1.0)/30.0))
				}
				energy += Min(MAX_NINIO, (nl-ns)*P.ninio[2])
				energy += P.mismatch1nI[type_1][si1][sj1] + P.mismatch1nI[type_2][sq1][sp1]
				return energy
			}
		} else if ns == 2 {
			if nl == 2 {
				/* 2x2 loop */
				return P.int22[type_1][type_2][si1][sp1][sq1][sj1]
			} else if nl == 3 {
				/* 2x3 loop */
				energy = P.internal_loop[5] + P.ninio[2]
				energy += P.mismatch23I[type_1][si1][sj1] + P.mismatch23I[type_2][sq1][sp1]
				return energy
			}
		}

		{
			/* generic interior loop (no else here!)*/
			u = nl + ns
			if u <= MAXLOOP {
				energy = P.internal_loop[u]
			} else {
				energy = P.internal_loop[30] + int(P.lxc*math.Log(float64(u)/30.0))
			}

			energy += Min(MAX_NINIO, (nl-ns)*P.ninio[2])

			energy += P.mismatchI[type_1][si1][sj1] + P.mismatchI[type_2][sq1][sp1]
		}
	}

	return energy
}

func ON_SAME_STRAND(I, J, C int) bool {
	if (I >= C) || (J < C) {
		return true
	} else {
		return false
	}
}

func E_IntLoop_Co(type_1, type_2, i, j, p, q, cutpoint, si1, sj1, sp1, sq1, dangles int,
	P *vrna_param_t) int {
	var e, energy, d3, d5, d5_2, d3_2, tmm, tmm_2 int
	var ci, cj, cp, cq bool

	energy = 0
	if type_1 > 2 {
		energy += P.TerminalAU
	}

	if type_2 > 2 {
		energy += P.TerminalAU
	}

	if !dangles {
		return energy
	}

	ci = ON_SAME_STRAND(i, i+1, cutpoint)
	cj = ON_SAME_STRAND(j-1, j, cutpoint)
	cp = ON_SAME_STRAND(p-1, p, cutpoint)
	cq = ON_SAME_STRAND(q, q+1, cutpoint)

	if ci {
		d3 = P.dangle3[type_1][si1]
	}

	if cj {
		d5 = P.dangle5[type_1][sj1]
	}

	if cp {
		d5_2 = P.dangle5[type_2][sp1]
	}

	if cq {
		d3_2 = P.dangle3[type_2][sq1]
	}

	if cj && ci {
		tmm = P.mismatchExt[type_1][sj1][si1]
	} else {
		tmm = d5 + d3
	}

	if cp && cq {
		tmm_2 = P.mismatchExt[type_2][sp1][sq1]
	} else {
		tmm_2 = d5_2 + d3_2
	}

	if dangles == 2 {
		return energy + tmm + tmm_2
	}

	/* now we may have non-double dangles only */
	if p-i > 2 {
		if j-q > 2 {
			/* all degrees of freedom */
			e = Min(tmm, d5)
			e = Min(e, d3)
			energy += e
			e = Min(tmm_2, d5_2)
			e = Min(e, d3_2)
			energy += e
		} else if j-q == 2 {
			/* all degrees of freedom in 5' part between i and p */
			e = Min(tmm+d5_2, d3+d5_2)
			e = Min(e, d5+d5_2)
			e = Min(e, d3+tmm_2)
			e = Min(e, d3+d3_2)
			e = Min(e, tmm_2) /* no dangles on enclosing pair */
			e = Min(e, d5_2)  /* no dangles on enclosing pair */
			e = Min(e, d3_2)  /* no dangles on enclosing pair */
			energy += e
		} else {
			/* no unpaired base between q and j */
			energy += d3 + d5_2
		}
	} else if p-i == 2 {
		if j-q > 2 {
			/* all degrees of freedom in 3' part between q and j */
			e = Min(tmm+d3_2, d5+d3_2)
			e = Min(e, d5+d3_2)
			e = Min(e, d3+d3_2)
			e = Min(e, d5+tmm_2)
			e = Min(e, tmm_2)
			e = Min(e, d5_2)
			e = Min(e, d3_2)
			energy += e
		} else if j-q == 2 {
			/* one possible dangling base between either side */
			e = Min(tmm, tmm_2)
			e = Min(e, d3)
			e = Min(e, d5)
			e = Min(e, d5_2)
			e = Min(e, d3_2)
			e = Min(e, d3+d3_2)
			e = Min(e, d5+d5_2)
			energy += e
		} else {
			/* one unpaired base between i and p */
			energy += Min(d3, d5_2)
		}
	} else {
		/* no unpaired base between i and p */
		if j-q > 2 {
			/* all degrees of freedom in 3' part between q and j */
			energy += d5 + d3_2
		} else if j-q == 2 {
			/* one unpaired base between q and j */
			energy += Min(d5, d3_2)
		}
	}

	return energy
}

/**
 *  @brief Evaluate free energy of a hairpin loop
 *
 *  @ingroup eval
 *
 *  @note This function is polymorphic! The provided #vrna_fold_compound_t may be of type
 *  #VRNA_FC_TYPE_SINGLE or #VRNA_FC_TYPE_COMPARATIVE
 *
 *  @param  fc  The #vrna_fold_compound_t for the particular energy evaluation
 *  @param  i   5'-position of the base pair
 *  @param  j   3'-position of the base pair
 *  @returns    Free energy of the hairpin loop closed by @f$ (i,j) @f$ in deka-kal/mol
 */
func vrna_eval_hp_loop(fc *vrna_fold_compound_t, i, j int) int {
	char * *Ss
	var Ss []string
	var a2s [][]uint64
	var S, S2 []int
	var SS, S5, S3 [][]int
	var sn []uint64
	var u, e, s, type_1, n_seq, en, noGUclosure int
	var P *vrna_param_t
	var md *vrna_md_t
	// //  vrna_ud_t         *domains_up;
	//  struct sc_hp_dat  sc_wrapper;

	P = fc.params
	md = &(P.model_details)
	noGUclosure = md.noGUclosure
	sn = fc.strand_number
	//  domains_up  = fc.domains_up
	e = INF

	if sn[j] != sn[i] {
		return eval_hp_loop_fake(fc, i, j)
	}

	//  vivek: handle sc stuff
	//  init_sc_hp(fc, &sc_wrapper);

	/* regular hairpin loop */

	S = fc.sequence_encoding
	S2 = fc.sequence_encoding2
	u = j - i - 1
	type_1 = vrna_get_ptype_md(S2[i], S2[j], md)

	if !(noGUclosure && ((type_1 == 3) || (type_1 == 4))) {
		e = E_Hairpin(u, type_1, S[i+1], S[j-1], fc.sequence[i-1:], P)
	}

	//  if e != INF {
	// 	 if (sc_wrapper.pair)
	// 		 e += sc_wrapper.pair(i, j, &sc_wrapper);

	// 	 /* consider possible ligand binding */
	// 	 if (domains_up && domains_up->energy_cb) {
	// 		 en = domains_up->energy_cb(fc,
	// 																i + 1, j - 1,
	// 																VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP,
	// 																domains_up->data);
	// 		 if (en != INF)
	// 			 en += e;

	// 		 e = MIN2(e, en);
	// 	 }
	//  }

	//  free_sc_hp(&sc_wrapper);

	return e
}

func eval_hp_loop_fake(fc *vrna_fold_compound_t, i, j int) int {
	var S, S2 []int
	var sn []uint64
	var u, e, ij, type_1, en, noGUclosure int
	var idx []int
	var P *vrna_param_t
	var sc *vrna_sc_t
	var md *vrna_md_t
	// vrna_ud_t     *domains_up;

	idx = fc.jindx
	P = fc.params
	md = &(P.model_details)
	noGUclosure = md.noGUclosure
	sn = fc.strand_number
	// domains_up  = fc.domains_up
	e = INF

	S = fc.sequence_encoding
	S2 = fc.sequence_encoding2
	sc = fc.sc
	u = j - i - 1
	ij = idx[j] + i
	type_1 = vrna_get_ptype_md(S2[j], S2[i], md)

	if noGUclosure && (type_1 == 3 || type_1 == 4) {
		return e
	}

	/* hairpin-like exterior loop (for cofolding) */
	var si, sj int

	if sn[i+1] == sn[i] {
		si = S[i+1]
	} else {
		si = -1
	}

	if sn[j] == sn[j-1] {
		sj = S[j-1]
	} else {
		sj = -1
	}

	switch md.dangles {
	case 0:
		e = vrna_E_ext_stem(type_1, -1, -1, P)

	case 2:
		e = vrna_E_ext_stem(type_1, sj, si, P)

	default:
		e = vrna_E_ext_stem(type_1, -1, -1, P)
		e = Min(e, vrna_E_ext_stem(type_1, sj, -1, P))
		e = Min(e, vrna_E_ext_stem(type_1, -1, si, P))
		e = Min(e, vrna_E_ext_stem(type_1, sj, si, P))
	}

	/* add soft constraints */
	// vivek: handle sc stuff
	// if (sc) {
	// 	if sc->energy_up {
	// 		e += sc->energy_up[i + 1][u];
	// 	}

	// 	if (sc->energy_bp)
	// 		e += sc->energy_bp[ij];

	// 	if (sc->f)
	// 		e += sc->f(i, j, i, j, VRNA_DECOMP_PAIR_HP, sc->data);
	// }

	/* consider possible ligand binding */
	// if (domains_up && domains_up->energy_cb) {
	// 	en = domains_up->energy_cb(fc,
	// 															i + 1, j - 1,
	// 															VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP,
	// 															domains_up->data);
	// 	if (en != INF)
	// 		en += e;

	// 	e = MIN2(e, en);
	// }

	return e
}

/**
 *  @brief Compute the Energy of a hairpin-loop
 *
 *  To evaluate the free energy of a hairpin-loop, several parameters have to be known.
 *  A general hairpin-loop has this structure:<BR>
 *  <PRE>
 *        a3 a4
 *      a2     a5
 *      a1     a6
 *        X - Y
 *        |   |
 *        5'  3'
 *  </PRE>
 *  where X-Y marks the closing pair [e.g. a <B>(G,C)</B> pair]. The length of this loop is 6 as there are
 *  six unpaired nucleotides (a1-a6) enclosed by (X,Y). The 5' mismatching nucleotide is
 *  a1 while the 3' mismatch is a6. The nucleotide sequence of this loop is &quot;a1.a2.a3.a4.a5.a6&quot; <BR>
 *  @note The parameter sequence should contain the sequence of the loop in capital letters of the nucleic acid
 *  alphabet if the loop size is below 7. This is useful for unusually stable tri-, tetra- and hexa-loops
 *  which are treated differently (based on experimental data) if they are tabulated.
 *  @see scale_parameters()
 *  @see vrna_param_t
 *  @warning Not (really) thread safe! A threadsafe implementation will replace this function in a future release!\n
 *  Energy evaluation may change due to updates in global variable "tetra_loop"
 *
 *  @param  size  The size of the loop (number of unpaired nucleotides)
 *  @param  type  The pair type of the base pair closing the hairpin
 *  @param  si1   The 5'-mismatching nucleotide
 *  @param  sj1   The 3'-mismatching nucleotide
 *  @param  string  The sequence of the loop (May be @p NULL, otherwise mst be at least @f$size + 2@f$ long)
 *  @param  P     The datastructure containing scaled energy parameters
 *  @return The Free energy of the Hairpin-loop in dcal/mol
 */
func E_Hairpin(size, type_1, si1, sj1 int, sequence string, P *vrna_param_t) int {
	var energy int

	if size <= 30 {
		energy = P.hairpin[size]
	} else {
		energy = P.hairpin[30] + int(P.lxc*math.Log(float64(size)/30.0))
	}

	if size < 3 {
		return energy /* should only be the case when folding alignments */
	}

	if P.model_details.special_hp == 1 {
		var tl string
		var idx int
		if size == 4 {
			tl = sequence[:6]
			// vivek: this could be a point of failure. Maybe change the above to 7
			// memcpy(tl, sequence, sizeof(char) * 6);
			idx = strings.Index(string(P.Tetraloops), tl)
			if idx != -1 {
				return P.Tetraloop_E[idx/7]
			}
		} else if size == 6 {
			tl = sequence[:8]
			idx = strings.Index(string(P.Hexaloops), tl)
			if idx != -1 {
				return P.Hexaloop_E[idx/9]
			}
		} else if size == 3 {
			tl = sequence[:5]
			idx = strings.Index(string(P.Triloops), tl)
			if idx != -1 {
				return P.Triloop_E[idx/6]
			}

			if type_1 > 2 {
				return energy + P.TerminalAU
			} else {
				return energy
			}
		}
	}

	energy += P.mismatchH[type_1][si1][sj1]

	return energy
}

func cut_in_loop(i int, pt []int, sn []uint64) int {
	/* walk around the loop;  return 5' pos of first pair after cut if
	 * cut_point in loop else 0 */
	var p, j int

	p, j = pt[i]

	for {
		i = pt[p]
		p = i + 1
		for pt[p] == 0 {
			p++
		}

		if !((p != j) && (sn[i] == sn[p])) {
			break
		}
	}
	if sn[i] == sn[p] {
		return 0
	} else {
		return p
	}
}

/**
 *** i is the 5'-base of the closing pair
 ***
 *** since each helix can coaxially stack with at most one of its
 *** neighbors we need an auxiliarry variable  cx_energy
 *** which contains the best energy given that the last two pairs stack.
 *** energy  holds the best energy given the previous two pairs do not
 *** stack (i.e. the two current helices may stack)
 *** We don't allow the last helix to stack with the first, thus we have to
 *** walk around the Loop twice with two starting points and take the minimum
 ***/
func energy_of_ml_pt(vc *vrna_fold_compound_t, i int, pt []int) (int, error) {
	var sn []uint64
	var energy, cx_energy, tmp, tmp2, best_energy, dangle_model, logML, circular, ss, n, n_seq int
	var idx, rtype []int
	best_energy = INF

	var i1, j, p, q, q_prev, q_prev2, u, uu, x, type_1, count, mm5, mm3, tt, ld5, new_cx,
		dang5, dang3, dang int
	var e_stem, e_stem5, e_stem3, e_stem53 int
	var mlintern [NBPAIRS + 1]int
	var s, s1 []int
	var S, S5, S3 [][]int
	var a2s [][]uint64
	var P *vrna_param_t
	var md *vrna_md_t
	var sc *vrna_sc_t
	var scs []*vrna_sc_t

	/* helper variables for dangles == 1|5 case */
	var E_mm5_available int  /* energy of 5' part where 5' mismatch of current stem is available */
	var E_mm5_occupied int   /* energy of 5' part where 5' mismatch of current stem is unavailable */
	var E2_mm5_available int /* energy of 5' part where 5' mismatch of current stem is available with possible 3' dangle for enclosing pair (i,j) */
	var E2_mm5_occupied int  /* energy of 5' part where 5' mismatch of current stem is unavailable with possible 3' dangle for enclosing pair (i,j) */

	n = vc.length
	sn = vc.strand_number
	P = vc.params
	md = &(P.model_details)
	idx = vc.jindx

	circular = md.circ
	dangle_model = md.dangles
	logML = md.logML
	rtype = &(md.rtype[0])
	s = vc.sequence_encoding2
	sc = vc.sc
	s1 = vc.sequence_encoding
	S = nil
	S5 = nil
	S3 = nil
	a2s = nil
	n_seq = 1
	scs = nil

	bonus = 0

	if i >= pt[i] {
		return INF, errors.New("energy_of_ml_pt: i is not 5' base of a closing pair!")
	}

	if i == 0 {
		j = n + 1
	} else {
		j = int(pt[i])
	}

	if i != 0 {
		/* (i,j) is closing pair of multibranch loop, add soft constraints */
		if sc {
			if sc.energy_bp {
				bonus += sc.energy_bp[idx[j]+i]
			}
		}
	}

	/* init the variables */
	energy = 0
	u = 0 /* the total number of unpaired nucleotides */
	p = i + 1
	q_prev = i - 1
	q_prev2 = i

	for x = 0; x <= NBPAIRS; x++ {
		mlintern[x] = P.MLintern[x]
	}

	/* seek to opening base of first stem */
	for p <= j && !pt[p] {
		p++
	}

	/* add bonus energies for first stretch of unpaired nucleotides */

	u += p - i - 1
	if sc {
		if sc.energy_up {
			bonus += sc.energy_up[i+1][u]
		}
	}

	switch dangle_model {
	case 0:
		for p < j {
			/* p must have a pairing partner */
			q = int(pt[p])
			/* get type of base pair (p,q) */
			tt = vrna_get_ptype_md(s[p], s[q], md)

			energy += E_MLstem(tt, -1, -1, P)

			/* seek to the next stem */
			p = q + 1
			q_prev = q
			q_prev2 = q
			for p < j && !pt[p] {
				p++
			}
			u += p - q - 1 /* add unpaired nucleotides */

			if sc {
				if sc.energy_up {
					bonus += sc.energy_up[q+1][p-q-1]
				}
			}
		}

		/* now lets get the energy of the enclosing stem */
		if i > 0 {
			/* actual closing pair */
			tt = vrna_get_ptype_md(s[j], s[i], md)

			energy += E_MLstem(tt, -1, -1, P)
		} else {
			/* virtual closing pair */
			energy += E_MLstem(0, -1, -1, P)
		}

	case 2:
		for p < j {
			/* p must have a pairing partner */
			q = int(pt[p])
			/* get type of base pair (p,q) */
			tt = vrna_get_ptype_md(s[p], s[q], md)

			if sn[p-1] == sn[p] {
				mm5 = s1[p-1]
			} else {
				mm5 = -1
			}

			if sn[q] == sn[q+1] {
				mm3 = s1[q+1]
			} else {
				mm3 = -1
			}

			energy += E_MLstem(tt, mm5, mm3, P)

			/* seek to the next stem */
			p = q + 1
			q_prev = q
			q_prev2 = q
			for p < j && !pt[p] {
				p++
			}
			u += p - q - 1 /* add unpaired nucleotides */

			if sc {
				if sc.energy_up {
					bonus += sc.energy_up[q+1][p-q-1]
				}
			}
		}
		if i > 0 {
			/* actual closing pair */
			tt = vrna_get_ptype_md(s[j], s[i], md)

			if sn[j-1] == sn[j] {
				mm5 = s1[j-1]
			} else {
				mm5 = -1
			}

			if sn[i] == s1[i+1] {
				mm3 = s1[i+1]
			} else {
				mm3 = -1
			}

			energy += E_MLstem(tt, mm5, mm3, P)
		} else {
			/* virtual closing pair */
			energy += E_MLstem(0, -1, -1, P)
		}

	case 3: /* we treat helix stacking different */
		for count = 0; count < 2; count++ {
			/* do it twice */
			ld5 = 0 /* 5' dangle energy on prev pair (type) */
			if i == 0 {
				j = uint64(pt[0] + 1)
				type_1 = 0 /* no pair */
			} else {
				j = uint64(pt[i])
				type_1 = vrna_get_ptype_md(s[j], s[i], md)

				/* prime the ld5 variable */
				if sn[j-1] == sn[j] {
					ld5 = P.dangle5[type_1][s1[j-1]]
					p = uint64(pt[j-2])

					if p == 1 && sn[j-2] == sn[j-1] {
						if P.dangle3[md.pair[s[p]][s[j-2]]][s1[j-1]] < ld5 {
							ld5 = 0
						}
					}
				}
			}

			i1 = i
			p = i + 1
			u = 0
			energy = 0
			cx_energy = INF

			for {
				/* walk around the multi-loop */
				new_cx = INF

				/* hop over unpaired positions */
				for p <= uint64(pt[0]) && pt[p] == 0 {
					p++
				}

				/* memorize number of unpaired positions */
				u += p - i1 - 1

				if sc {
					if sc.energy_up {
						bonus += sc.energy_up[i1+1][p-i1-1]
					}
				}

				/* get position of pairing partner */
				if p == uint64(pt[0])+1 {
					q = 0
					tt = 0 /* virtual root pair */
				} else {
					q = uint64(pt[p])
					/* get type of base pair P->q */
					tt = vrna_get_ptype_md(s[p], s[q], md)
				}

				energy += mlintern[tt]
				cx_energy += mlintern[tt]

				dang5 = 0
				dang3 = 0
				if (sn[p-1] == sn[p]) && p > 1 {
					dang5 = P.dangle5[tt][s1[p-1]] /* 5'dangle of pq pair */
				}

				if (sn[i1] == sn[i1+1]) && i1 < uint64(s[0]) {
					dang3 = P.dangle3[type_1][s1[i1+1]] /* 3'dangle of previous pair */
				}

				switch p - i1 - 1 {
				case 0:
					/* adjacent helices */
					if i1 != 0 {
						if sn[i1] == sn[p] {
							new_cx = energy + P.stack[rtype[type_1]][rtype[tt]]
							/* subtract 5'dangle and TerminalAU penalty */
							new_cx += -ld5 - mlintern[tt] - mlintern[type_1] + 2*mlintern[1]
						}

						ld5 = 0
						energy = Min(energy, cx_energy)
					}

				case 1: /* 1 unpaired base between helices */
					dang = Min(dang3, dang5)
					energy = energy + dang
					ld5 = dang - dang3
					/* may be problem here: Suppose
					* cx_energy>energy, cx_energy+dang5<energy
					* and the following helices are also stacked (i.e.
					* we'll subtract the dang5 again */
					if cx_energy+dang5 < energy {
						energy = cx_energy + dang5
						ld5 = dang5
					}

					new_cx = INF /* no coax stacking with mismatch for now */
				default: /* many unpaired base between helices */
					energy += dang5 + dang3
					energy = Min(energy, cx_energy+dang5)
					new_cx = INF /* no coax stacking possible */
					ld5 = dang5
				}
				type_1 = tt
				cx_energy = new_cx
				i1 = q
				p = q + 1

				if q == i {
					break
				}
			}
			best_energy = Min(energy, best_energy) /* don't use cx_energy here */
			/* skip a helix and start again */
			for pt[p] == 0 {
				p++
			}
			if i == uint64(pt[p]) {
				break
			}

			i = uint64(pt[p])
		} /* end doing it twice */

		energy = best_energy

	default:
		E_mm5_available = INF
		E2_mm5_available = INF
		E_mm5_occupied = 0
		E2_mm5_occupied = 0
		for p < j {
			/* p must have a pairing partner */
			q = int(pt[p])
			/* get type of base pair (p,q) */
			tt = vrna_get_ptype_md(s[p], s[q], md)

			if q_prev+2 < p {
				E_mm5_available = Min(E_mm5_available, E_mm5_occupied)
				E_mm5_occupied = E_mm5_available
			}

			if q_prev2+2 < p {
				E2_mm5_available = Min(E2_mm5_available, E2_mm5_occupied)
				E2_mm5_occupied = E2_mm5_available
			}

			if (sn[p-1] == sn[p]) && !pt[p-1] {
				mm5 = s1[p-1]
			} else {
				mm5 = -1
			}

			if (sn[q] == sn[q+1]) && !pt[q+1] {
				mm3 = s1[q+1]
			} else {
				mm3 = -1
			}

			e_stem = E_MLstem(tt, -1, -1, P)
			e_stem5 = E_MLstem(tt, mm5, -1, P)
			e_stem3 = E_MLstem(tt, -1, mm3, P)
			e_stem53 = E_MLstem(tt, mm5, mm3, P)

			tmp = E_mm5_occupied + e_stem3
			tmp = Min(tmp, E_mm5_available+e_stem53)
			tmp = Min(tmp, E_mm5_available+e_stem3)
			tmp2 = E_mm5_occupied + e_stem
			tmp2 = Min(tmp2, E_mm5_available+e_stem5)
			tmp2 = Min(tmp2, E_mm5_available+e_stem)

			E_mm5_occupied = tmp
			E_mm5_available = tmp2

			tmp = E2_mm5_occupied + e_stem3
			tmp = Min(tmp, E2_mm5_available+e_stem53)
			tmp = Min(tmp, E2_mm5_available+e_stem3)
			tmp2 = E2_mm5_occupied + e_stem
			tmp2 = Min(tmp2, E2_mm5_available+e_stem5)
			tmp2 = Min(tmp2, E2_mm5_available+e_stem)

			E2_mm5_occupied = tmp
			E2_mm5_available = tmp2

			/* seek to the next stem */
			p = q + 1
			q_prev = q
			q_prev2 = q
			for p < j && !pt[p] {
				p++
			}
			u += p - q - 1 /* add unpaired nucleotides */

			if sc {
				if sc.energy_up {
					bonus += sc.energy_up[q+1][p-q-1]
				}
			}
		}

		if i > 0 {
			/* actual closing pair */
			type_1 = vrna_get_ptype_md(s[j], s[i], md)

			if (sn[j-1] == sn[j]) && !pt[j-1] {
				mm5 = s1[j-1]
			} else {
				mm5 = -1
			}

			if (sn[i] == sn[i+1]) && !pt[i+1] {
				mm3 = s1[i+1]
			} else {
				mm3 = -1
			}

			if q_prev+2 < p {
				E_mm5_available = Min(E_mm5_available, E_mm5_occupied)
				E_mm5_occupied = E_mm5_available
			}

			if q_prev2+2 < p {
				E2_mm5_available = Min(E2_mm5_available, E2_mm5_occupied)
				E2_mm5_occupied = E2_mm5_available
			}

			e_stem = E_MLstem(type_1, -1, -1, P)
			e_stem5 = E_MLstem(type_1, mm5, -1, P)
			e_stem3 = E_MLstem(type_1, -1, mm3, P)
			e_stem53 = E_MLstem(type_1, mm5, mm3, P)
		} else {
			/* virtual closing pair */
			e_stem = E_MLstem(0, -1, -1, P)
			e_stem5 = e_stem
			e_stem3 = e_stem
			e_stem53 = e_stem
		}

		/* now lets see how we get the minimum including the enclosing stem */
		energy = E_mm5_occupied + e_stem
		energy = Min(energy, E_mm5_available+e_stem5)
		energy = Min(energy, E_mm5_available+e_stem)
		energy = Min(energy, E2_mm5_occupied+e_stem3)
		energy = Min(energy, E2_mm5_occupied+e_stem)
		energy = Min(energy, E2_mm5_available+e_stem53)
		energy = Min(energy, E2_mm5_available+e_stem3)
		energy = Min(energy, E2_mm5_available+e_stem5)
		energy = Min(energy, E2_mm5_available+e_stem)
	} /* end switch dangle_model */

	energy += P.MLclosing

	/*
	* logarithmic ML loop energy if logML
	* does this work for comparative predictions as well?
	 */
	if logML && (u > 6) {
		energy += 6*P.MLbase + int(P.lxc*math.Log(float64(u)/6.0))
	} else {
		energy += u * P.MLbase
	}

	return energy + bonus
}

/**
 *  @def E_MLstem(A,B,C,D)
 *  <H2>Compute the Energy contribution of a Multiloop stem</H2>
 *  This definition is a wrapper for the E_Stem() funtion.
 *  It is substituted by an E_Stem() funtion call with argument
 *  extLoop=0, so the energy contribution returned reflects a
 *  stem introduced in a multiloop.<BR>
 *  As for the parameters B (si1) and C (sj1) of the substituted
 *  E_Stem() function, you can inhibit to take 5'-, 3'-dangles
 *  or mismatch contributions to be taken into account by passing
 *  -1 to these parameters.
 *
 *  @see    E_Stem()
 *  @param  A The pair type of the stem-closing pair
 *  @param  B The 5'-mismatching nucleotide
 *  @param  C The 3'-mismatching nucleotide
 *  @param  D The datastructure containing scaled energy parameters
 *  @return   The energy contribution of the introduced multiloop stem
 */
func E_MLstem(type_1, si1, sj1, P *vrna_param_t) int {
	var energy int = 0

	if si1 >= 0 && sj1 >= 0 {
		energy += P.mismatchM[type_1][si1][sj1]
	} else if si1 >= 0 {
		energy += P.dangle5[type_1][si1]
	} else if sj1 >= 0 {
		energy += P.dangle3[ttype_1ype][sj1]
	}

	if type_1 > 2 {
		energy += P.TerminalAU
	}

	energy += P.MLintern[type_1]

	return energy
}
