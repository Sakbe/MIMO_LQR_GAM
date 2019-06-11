
//******************************************************************************
//
//		LQR.cpp - LQR CONTROL AND KALMAN FILTER
//		CORONA - CIANCIULLI 05/2019
//
//******************************************************************************


#include "LQR.h"

//OBJECTLOADREGISTER(IPID,"$Id: IPID.cpp,v 1.0 29/4/2011 14:22:36 ivoc Exp $")
// if cycle time is supplied
//CONSTRUCTOR AND DESTRUCTOR
LQR::LQR(){
	//this-> N_state = 10;
	//this-> N_input = 2;
	//this-> N_output = 2;
	
	this-> X_LQR = (float[N_state]){0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	
	this-> x_dot_pos = (float[N_state]){0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	
	this-> x_pos = (float[N_state]){0.0983539734804260, 0.480564466137287, -0.204246862173158, 0.0605647026605041, -0.369796209977501, 0.111380460346364, 0.0194246756583234, 0.144179620264490, -0.0534417070457820, 0.215640033094750}; 
	
	//this-> K_LQR_pos =(float[N_input * N_state]) {-39.8874724584811, -70.6769666120137, 354.869808121735, 370.159896447712, -106.569303855515, -10.6727031034049, -29.8250480200422, 75.3583282160950, -135.655340474040, -27.8586258177635,
	//						4.85034537803885, 16.3699444462248, -103.189045436383, -88.9493312103841, 36.1506690550597, 19.0337453517795, 10.4943236824918, -28.3848265454111, 45.6326833691854, 8.77044127460252};
	
	this-> A_est_pos= (float[N_state * N_state]) {1.13163023314199, -0.343597594340211, -0.213879462160118, 0.0179916185487389, -0.150369187380823, 0.0269299858719678, -0.0433137874319706, 0.179179577153557, -0.207578023558474, -0.0352702122199913,
																0.287600842583737, 0.540556228880367, -0.204154977999538, 0.0102852640208614, -0.209928800555951, 0.0360517621845844, 0.0666388339639277, 0.178351267657655, -0.171995165509242, -0.127345902924508,
																0.0754707536107985, -0.238668989358296, 0.818854849839767, -0.164733743815283, -0.0217690782606367, 0.0146591987555307, -0.0324172794802642, 0.0606121920942799, -0.179686741383566, -0.0508405095031112,
																0.131925929817064, -0.372248680592072, 0.0276052171048853, 1.02179796923870, -0.257351929706168, 0.0179101866506223, -0.0840360783373023, 0.171227557561418, -0.287524674661170, 0.000517459687887203,
																0.177657933921153, -0.444331155707191, -0.340898331624331, 0.221398005868031, 0.680207481216171, 0.414718779863749, 0.0495581812996120, -0.0146824110541471, -0.256160870498049, -0.0879751273508828,
																0.00428166489028195, 0.00468605824926672, -0.118841222744660, -0.0775930411311325, -0.470877574197134, 0.736438110295846, 0.292243145826319, -0.106667041723583, 0.184354045805486, -0.0352856736421628,
																0.0211465169994229, 0.0920804800367413, 0.0775101736851033, -0.0871076469823266, -0.00274664766965684, 0.301191024811585, -0.853677736106622, -0.338194597621530, 0.0903540072889142, 0.126857877037311,
																-0.101039074797006, 0.172115106521145, 0.126989777403869, -0.113977526679146, 0.0303493831666241, 0.225517804060444, 0.120367523058471, -0.182536009695298, -0.773100535530777, 0.197640938910675,
																-0.0125352928445362, -0.0356027448693943, -0.00986410114649150, 0.0104414601859412, -0.00612203568636584, -0.193714169340406, 0.0971839573892438, -0.508057156805493, 0.207904559974678, 0.494301295754186,
																-0.0681141671749747, 0.0631731918311591, 0.0391916236690946, -0.0709854208561612, -0.00992959690135519, -0.0263008458827116, 0.103215896765230, -0.0365368432532828, -0.136301181420319, 0.0809720948251243};
	//B_est_pos =[vertical; horizonal; rc; zc]
	
	this-> B_est_pos=(float[N_state * ( N_input +  N_output)])  {-1.16125013575191*1e-05, 0.000163758480113681, 0.000257613772728220, 0.000278483718207995, -0.000369222690232994, 0.000391024178414270, 2.19403868131743*1e-05, -0.000254502475896278, -0.000152028009309686, -0.00118431170706074,
																					2.59526548636673*1e-05, 0.000134280648824691, -0.000550079593989017, 7.06765132488954*1e-05, 0.00132298996363479, 0.00167244079903906, -0.000296832394696530, -0.000922111433587482, -0.00107176586288692, -0.00101009010026278,
																					0.876479961841367, 1.01242478492612, 0.574210915670357, 0.964761414967758, 1.23892804780780, -0.0130977878038005, -0.797987484859047, -0.416006708588100, 0.196839888372422, -0.161860817345070, 
																					-2.30413733637204, -2.86691484603969, -1.51143528695349, -2.25683449771732, -3.09391877161166, -0.208771202482571, 0.765181314330032, 1.36998563604846, -0.376017766040073, 0.592294754381924};
	
	
	this-> C_est_pos=(float[N_state * (N_state +  N_output)])  {0.172597562096771, -0.0406841025944168, -0.0424424164905047, 0.000820643421739336, -0.0326833771976484, -0.00365572309023441, -0.000568188978630238, 0.0249183302168008, -0.0404898910313346, -0.00958025363872561,
																					0.0903364751248640, -0.0663772722506334, -0.0471458612175972, 0.00301323332985550, -0.0332791814731680, -0.000992351893563895, -0.00420353785437508, 0.0333815426885655, -0.0451543116864300, -0.0101756896082312,
																					1.18884773817693, -0.358654960655223, -0.220162298985963, 0.0189011564774812, -0.148452934033422, 0.00241780432962949, -0.0278396433717354, 0.168989127864764, -0.211269724974846, -0.0464470550704881,
																					0.229509553160214, 0.604796934780455, -0.245072346293721, 0.0206397526839626, -0.165825390855866, 0.00210741024042335, -0.0303096856161090, 0.187024152831786, -0.235139879194527, -0.0517909904013656,
																					0.167787245300399, -0.317882162249978, 0.804819138644354, 0.0167488520649557, -0.131619210631909, 0.00213233078661748, -0.0246677698542292, 0.149793541644763, -0.187296669169996, -0.0411784786310142,
																					0.167164609215762, -0.351089470722029, -0.213473006599433, 1.01865739890114, -0.143466499241793, 0.00282692400818382, -0.0275556055487798, 0.164751369201408, -0.204878512085112, -0.0449625071255927,
																					0.200175764880664, -0.396987340766533, -0.242669570947558, 0.0209988010405207, 0.836608787354835, 0.00290642839112818, -0.0309668006232276, 0.186713604156288, -0.232881861304842, -0.0511586671634703,
																					0.178051672848331, -0.265978025575213, -0.167663221297088, 0.0136843925931207, -0.114075245525449, 1.00080512191246, -0.0199948167120765, 0.126767288948896, -0.160831250018297, -0.0355289006437349,
																					0.00383589029017545, -0.150145747826621, -0.0834758634882206, 0.00857120657092464, -0.0542641781559624, 0.00296767662899450, 0.987056687456786, 0.0678842564372292, -0.0802227318163792, -0.0172990494480662,
																					0.108985599301217, -0.174203909575089, -0.108930345887065, 0.00902947312509098, -0.0739145819293058, 0.000725681099243316, -0.0132264431869919, 1.08273682441609, -0.104503342109657, -0.0230522482456343,
																					0.165699441469261, -0.268112493366233, -0.167416316986419, 0.0139148512567483, -0.113546403204613, 0.00116983872218807, -0.0203913602166371, 0.127260516498298, 0.839384444871921, -0.0354209830451735,
																					0.134952441212368, -0.216695226800387, -0.135429073714409, 0.0112373141887950, -0.0918789918794712, 0.000918716714456840, -0.0164631624775202, 0.102894303993053, -0.129926097424931, 0.971342499342926};
	this-> D_est_pos=(float[(N_state + N_output) * (N_input + N_output)]) { 0, 							0, 					0, 					0, 				0,					 0, 			0, 					0, 					0, 					0,			 0, 					0,
																			0, 							0, 					0, 					0, 				0, 					 0, 			0, 					0, 					0, 					0,			 0, 					0,
																			0.0837226767933073, -0.201930464264813, 0.932065758525171, 0.882254761320045, 0.823346631120921, 1.03206071769100, 1.09155104677677, 0.434307167651101, 0.898729664509743, 0.336036710307567, 0.530958079293443, 0.422170361983903,
																			-0.201930464264813, 0.581221747348804, -2.40253906173113, -2.63736312680753, -2.12922174119470, -2.36012581017436, -2.66345378002374, -1.76397174672629, -1.04091890922158, -1.15888747773842, -1.78456363151675, -1.44184795296202};


//	this-> N_BAR_pos = (float[N_input +  N_output]){-16738.8719868375, 3479.82442345434,
//							488.651395096507, -28.3958381254954};


	
	//more slow
	this-> N_BAR_pos = (float[N_input +  N_output]){-10292.2693526553, 2133.74129809333,
							-944.140351815484, 269.810568562941};
	this-> K_LQR_pos =(float[N_input * N_state]) {-0.677239930551175, -5.56445051715661, 29.1530400398665, 10.5689835482429, -6.43042836144154, 0.784006911127944, -1.67848331436731, 6.89611601694075, -11.6198347798690, -4.41024916164279,
							-0.253637909222695, 0.998896201615173, -6.91291373508615, -4.34730845911924, 2.02654430640833, 0.350198601401546, 0.540286005706202, -1.95918796444674, 3.18376611395890, 0.992300688479676};
	


	this-> x_neg=(float[N_state])  {-0.0494766508450061, -0.111034664506845, 0.0242441680811973, 0.00586343293525677, -0.00854660070320271, 0.000970048477408260, -0.000888992834197872, 0.00166157894650906, 0.000997825822154208, 0.00231148829292225};
	this-> x_dot_neg=(float[N_state])  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
//	this-> K_LQR_neg=(float[N_input * N_state])  {-1.10243700688381, -26.6888168111469, 248.208386299799, 2272.69314089294, 341.435985384170, 116.974185358079, -1571.00753489878, -331.351679250269, -1271.06969797964, -4441.09340489016,
//							-0.233414809301534, -21.1374545809135, 237.120488501101, 3170.68925311282, -143.026532612836, 258.466914847382, -1451.21045848359, -298.057584698670, -925.765520002125, -4247.87411829111};
	this-> A_est_neg=(float[N_state * N_state])  {0.503943697958563, -0.703733733757264, -1.20816317477897, 0.428531919950512, 1.23129699755514, -0.561852362663766, 0.338388473423280, -0.438150205967471, 0.289597040968814, -0.0480026201041203,
												0.0165579491300655, 0.757869954005823, -0.552904529698174, -0.281944731967268, 0.356868065244389, -0.0494763712769074, -0.255634206854070, -0.191293066044900, 0.131948522189369, -0.227979368614917,
												-0.00855669962600962, -0.278452952597036, 0.587563230437082, -0.155539594933950, 0.905217586514629, -0.335361092999183, 0.355018445302572, -0.198224607540228, 0.221100552953902, -0.553699095984438,
												-0.120098896877763, -0.159716003758573, -0.234950952213230, 1.08154444189073, 0.268636173784699, -0.122520411580655, -0.137925333022668, -0.220787505809466, 0.00256955633241329, -0.629883316799911,
												-0.0838005900304720, -0.160344184303423, -0.318499811799575, 0.186491951640922, 0.884720365627487, -0.507855653942687, 0.484217572469569, 0.00989922031363996, -0.550380188727969, 0.0440915577571055,
												0.0915067248521410, 0.0625476554578888, 0.108898125714359, -0.0363552767766716, 0.0566053869752655, -0.224517765317217, -0.843288628891406, -0.0911441079948653, 0.817987619265961, 0.293664777556748,
												0.127851733466733, -0.346629185846461, -0.369646152710029, -0.258214338111236, 0.564907037971381, 0.317764110980288, 0.159404262220738, 0.211403951847015, 0.382461602028459, -0.409769213831037,
												-0.0491967501924090, 0.0649212718847929, 0.0570751713799720, 0.0412564432537669, 0.000896657382601702, 0.361210554206455, -0.118385677655563, -0.812871753541494, -0.263437437325031, 0.542236840982197,
												-0.0182123844229220, -0.0135169925960958, 0.0180740969222822, 0.0652218701092334, 0.00398035478490197, -0.0141984591587385, -0.354436816769253, 0.304360358719350, -0.347047883770401, 0.375325178023404,
												0.0517900539802355, -0.0802994317903856, -0.0895635516356528, -0.0174992027900167, -0.0162819302261394, -0.372208196109113, 0.229358029899087, -0.226568642205041, 0.228543644273500, 0.445259535912489};
	this-> B_est_neg=(float[N_state * ( N_input +  N_output)]) {-4.66761235181270*1e-06, 4.44690659267741*1e-06, 2.29092031626160*1e-05, 2.82352037279260*1e-06, -8.49578049481197*1e-06, 4.65936717166927*1e-06, -6.20676767469788*1e-06, 1.96481773957388*1e-06, -2.29259339772744*1e-06, 1.49424416715292*1e-06,
																-2.38310588446242*1e-05, -1.27920234547367*1e-05, -6.97286578383591*1e-07, 6.67115753894253*1e-06, -1.98963150726010*1e-05, -1.03834671682820*1e-05, 7.07537647138095*1e-07, -4.38712299400625*1e-06, 5.28217038799443*1e-06, -3.29207149631720*1e-06,
																-2.10897005759671, 		-0.568945831911984,		 -0.735563919109126,		 -0.473256656052601, 	-0.472693429841058, 		0.213606910229934, 	-0.884477437212619, 		0.155047714906280, 		-0.0390239568417088, 	-0.192211895530816,
																1.27810843768912, 		-0.137655568424764,		 -0.0808648703360026,		 0.315333616407226, 		0.199917506965288, 		-0.261884796632875, -0.533478687078687, 		0.178569850729397, 		0.0525630937085863, 		-0.193056813631261};
	
	this-> C_est_neg=(float[N_state * (N_state +  N_output)]) {-0.0491389528881855, -0.0690239621676708, -0.110783527259196, 0.0151963866422061, 0.118048216983325, -0.0577296202729347, 0.0329931294883742, -0.0423543191424276, 0.0278840340937851, 0.00239413345051835,
																0.115523694907464, 0.0292595696392667, 0.0918664165271154, -0.105518601963745, -0.107725393285383, 0.0297392946178057, -0.0836961390943443, 0.0190337432883174, -0.0607817916127189, -0.123356138821500,
																0.519043294802380, -0.628054978186631, -1.02407483557455, 0.173676675377220, 1.09474339463005, -0.527169023212924, 0.325117110993426, -0.385771360086981, 0.271214851684188, 0.0655010427214920,
																-0.00721117338615529, 0.677708233745577, -0.411893504329264, -0.161562968445698, 0.415822804477084, -0.257193446785139, -0.00954593358591677, -0.195230341962687, 0.0152922515186002, -0.275938577849891,
																-0.0175447015150904, -0.162995876131127, 0.785098463993159, -0.0671678342802382, 0.218764320197990, -0.130846123182222, 0.00540328174411952, -0.0988942577413839, 0.0149197756688195, -0.121597750135144,
																-0.0489030153656803, -0.308708196828379, -0.414448404905380, 0.889186128617969, 0.423879674246704, -0.248689633149408, 0.0217722944178242, -0.187480780186624, 0.0363618064163303, -0.210051314152982,
																-0.00165609996536642, -0.203419826251490, -0.258600447568362, -0.105002540378134, 1.26068974604691, -0.162170877894514, -0.00815653801951600, -0.123189879519053, 0.00815490659440868, -0.177903950435155,
																0.0551419433619698, -0.202504643494264, -0.230506228528832, -0.163949536777454, 0.224921902121195, 0.841717777100605, -0.0499271985243076, -0.121988207145135, -0.0212455766994582, -0.250475384017951,
																0.0326201171457825, -0.207061541588684, -0.246961779345526, -0.142776211247045, 0.244458702691437, -0.163165827156559, 0.966442596005071, -0.125004176326625, -0.00943708026354638, -0.225411065389052,
																0.0638199900517775, -0.153185621787812, -0.163883868530924, -0.147151372843481, 0.156676026902523, -0.118503592825000, -0.0540425554449654, 0.907973466197798, -0.0275021389127691, -0.218035707002714,
																0.135027025265622, -0.0879217581154608, -0.0474009561400222, -0.187410303169312, 0.0299846952219131, -0.0625423393179604, -0.103454782030276, -0.0516973718678192, 0.933338350755660, -0.252269173551680,
																0.0704663385888843, -0.149979440921235, -0.156668563901329, -0.152423080034649, 0.148534350781777, -0.115579293248765, -0.0587876032940097, -0.0900094147420284, -0.0310536953649612, 0.776215139449839};
	this-> D_est_neg=(float[(N_state + N_output) * (N_input + N_output)]){0,					 		0, 					0, 				0, 						0, 				0, 					0, 					0, 						0, 				0, 						0, 					0,
																		  0,						 	0, 					0, 				0, 						0, 				0, 					0, 					0, 						0, 				0, 						0,					0,
																		  0.792748001070510, 0.127200492028566, -1.89986323015272, -0.875378185180234, -0.448488414111140, -0.855933161747806, -0.551306520752944, -0.525229532468674, -0.546922069627300, -0.388126656171121, -0.181883164704487, -0.376686561697634,
																		  0.127200492028566, 0.649838802214520, 1.26326480026836, -0.101299827329797, -0.00775356577001699, 0.0343497717108373, -0.0729949361324825, -0.250330990815615, -0.181624658268750, -0.258525806346280, -0.456209422157503, -0.278086171487326};
//	this-> N_BAR_neg=(float[N_input +  N_output]) {6692.27395164254, -960.635617013133,
	//						-2256.41280370701, -489.068332481912};


	//more slow
	this-> N_BAR_neg = (float[N_input +  N_output]){6577.20274262693, -953.521247320254,
							-2537.02933847542, -432.837113235285};
	this-> K_LQR_neg =(float[N_input * N_state]) {-0.0143959540956107, -0.325150394017629, 3.13416013810849, 26.5396916000373, 4.19324806144036, 1.55806909562497, -18.7468162690424, -3.85545880019961, -15.3160343140834, -53.4859040919299,
							-0.00458150757613837, -0.243881430718078, 2.88347865234915, 35.3882054725928, -2.01002787695297, 2.81135503841810, -16.8859632308927, -3.47894536276546, -10.3473795912101, -48.0323686610169};
	
	
}

LQR::~LQR(){
}

//////////////////////////// Controllers /////////////////////////////////////////////////////////


LQRouputs  LQR::MIMO_CONTROL_POSITIVE(float R_ref, float Z_ref, float R_real, float Z_real, float I_vertical, float I_horizontal){
	float u_Nbar_0 = 0;
	float u_Nbar_1 = 0;
	float u_vertical=0;
	float u_horizontal=0;
	int i = 0;
	int j = 0;
	float temp = 0.0;
	float *outputs;
	outputs =(float[2]) {0, 0};
	LQRouputs Outputs={0,0};
	
	
	
	//Nbar(2x2)
	u_Nbar_0 = (R_ref* this->N_BAR_pos[0]) + (Z_ref* this->N_BAR_pos[1]);
	u_Nbar_1 = (R_ref * this->N_BAR_pos[2]) + (Z_ref * this->N_BAR_pos[3]);
	
	//u=u-k_lqr*x
	for (i = 0; i < N_state; i++) {
		temp += this->X_LQR[i] * this->K_LQR_pos[i];
		}
	u_vertical = u_Nbar_0 - temp;
	temp = 0;		
	
	for (i = 0; i < N_state; i++) {
		temp += this->X_LQR[i] * this->K_LQR_pos[i+N_state];
		
	}
	u_horizontal = u_Nbar_1 - temp;
	temp = 0;
	
	Outputs.Ivert = u_vertical;
	Outputs.Ihor = u_horizontal;
	
	return Outputs;
}


LQRouputs  LQR::MIMO_CONTROL_NEGATIVE(float R_ref, float Z_ref, float R_real, float Z_real, float I_vertical, float I_horizontal){
	float u_Nbar_0 = 0;
	float u_Nbar_1 = 0;
	float u_vertical=0;
	float u_horizontal=0;
	int i = 0;
	int j = 0;
	float temp = 0.0;
	float *outputs;
	outputs =(float[2]) {0, 0};
	LQRouputs Outputs={0,0};
	
	
	
	//Nbar(2x2)
	u_Nbar_0 = (R_ref* this->N_BAR_neg[0]) + (Z_ref* this->N_BAR_neg[1]);
	u_Nbar_1 = (R_ref * this->N_BAR_neg[2]) + (Z_ref * this->N_BAR_neg[3]);
	
	//u=u-k_lqr*x
	for (i = 0; i < N_state; i++) {
		temp += this->X_LQR[i] * this->K_LQR_neg[i];	
	}
	u_vertical = u_Nbar_0 - temp;
	temp = 0.0;	
	
	for (i = 0; i < N_state; i++) {
		temp += this->X_LQR[i] * this->K_LQR_neg[i+N_state];
		
	}
	u_horizontal = u_Nbar_1 - temp;
		temp = 0.0;
		
	Outputs.Ivert = u_vertical;
	Outputs.Ihor = u_horizontal;
	
	return Outputs;
}




///////////////////////////////////////////// KALMAN FILTERS////////////////////////////////////////////////////////





Kalman LQR:: KALMAN_FILTER_POS(float R_real, float Z_real, float I_vertical, float I_horizontal, int sign){
	//the same as before
	//if we use this in mimo class we can comment from 102 to 115
	//i putted because we can copy this in IPID class
	float temp=0.0;
	int j=0;
	int i=0;
	int m=0;
	int n=0;
	// buffers time
	
	float buff=0.0;
	float buffer=0.0;
	
	
	float* y_est;
	y_est =(float[2]){0, 0};
	float* X_est;
	X_est =(float[10]){0, 0,0,0,0,0,0,0,0,0};
	
	
	Kalman Outputs={0.085,0.085,X_est};
	
	for(i=0; i< N_state; i++){
		buff = this->x_pos[i];
		this->X_LQR[i] = buff;
	}
	buff=0.0;	
			
	if(sign == 1){	
		
		for(j = 0; j < N_state; j++){
			buff = 0.0;
			for (i = 0; i < N_state; i++) {
				m=i + j*(N_state);
				buff = buff + this->A_est_pos[m] * this->x_pos[i];
			}
			this->x_dot_pos[j] = buff;
		}
		
		for(i = 0; i < N_state; i++){
			j= i + (N_state);
			m=i + 2*(N_state);
			n=i + 3*(N_state);
			temp =  this->B_est_pos[i] * I_vertical + this->B_est_pos[j] * I_horizontal + this->B_est_pos[m] * R_real + this->B_est_pos[n] * Z_real;
			this->x_dot_pos[i] += temp;
			
		}
		temp = 0.0;
		
		//x^(k+1) = A_est * x(k)+ B_est * u(k)
				
		//[y_est x_est] = C_est x_est + D_est [u_real y_real]
		for(j = 0; j < N_output; j++){
			buff = 0.0;
			for (i = 0; i < N_state; i++) {
				m=i + j*(N_state);
				 buff += this->C_est_pos[m] * this->x_pos[i];
			}
			y_est[j] = buff;
		}
		
	

			buff =  this->D_est_pos[24] * R_real + this->D_est_pos[36] * Z_real;
			y_est[0] += buff;
			buff = 0.0;
			buff =  this->D_est_pos[25] * R_real + this->D_est_pos[37] * Z_real;
			y_est[1] += buff;
			buff = 0.0;
			
			

		

		
			//x(k)=x(k+1) for the next step

		
		Outputs.Kalman_R= y_est[0];
		Outputs.Kalman_Z= y_est[1];
		Outputs.X_est= this-> x_pos;
		
		buffer = 0.0;
		for (i = 0; i < N_state; i++) {
			buffer =this->x_dot_pos[i];
			this->x_pos[i]=buffer;
		}
		
		}else{
			Outputs.Kalman_R= 0.085;
			Outputs.Kalman_Z= 0.085;
				for (i = 0; i < N_state; i++) {
						X_est[i]=0.0;}
										
			Outputs.X_est=X_est;
			}
	

	return Outputs;
	
}

Kalman LQR:: KALMAN_FILTER_NEG(float R_real, float Z_real, float I_vertical, float I_horizontal, int sign){
	//the same as before
	//if we use this in mimo class we can comment from 102 to 115
	//i putted because we can copy this in IPID class
	float temp=0.0;
	int j=0;
	int i=0;
	int m=0;
	int n=0;
	// buffers time
	
	float buff=0.0;
	float buffer=0.0;
	

	float* y_est;
	y_est =(float[2]){0, 0};
	float* X_est;
	X_est =(float[10]){0, 0,0,0,0,0,0,0,0,0};
	Kalman Outputs={0.085,0.085,X_est};
	
	for(i=0; i< N_state; i++){
		buff = this->x_neg[i];
		this->X_LQR[i] = buff;;
	}
	buff=0.0;	
	if(sign == 1){	
		
		for(j = 0; j < N_state; j++){
			buff = 0.0;
			for (i = 0; i < N_state; i++) {
				m=i + j*(N_state);
				buff = buff + this->A_est_neg[m] * this->x_neg[i];
			}
			this->x_dot_neg[j] = buff;
		}
		
		for(i = 0; i < N_state; i++){
			j= i + (N_state);
			m=i + 2*(N_state);
			n=i + 3*(N_state);
			temp =  this->B_est_neg[i] * I_vertical + this->B_est_neg[j] * I_horizontal + this->B_est_neg[m] * R_real + this->B_est_neg[n] * Z_real;
			this->x_dot_neg[i] += temp;
			
		}
		temp = 0.0;

		
		
		//[y_est x_est] = C_est x_est + D_est [u_real y_real]
		for(j = 0; j < N_output; j++){
			buff = 0.0;
			for (i = 0; i < N_state; i++) {
				m=i + j*(N_state);
				 buff += this->C_est_neg[m] * this->x_neg[i];
			}
			y_est[j] = buff;
		}
		


			buff =  this->D_est_neg[24] * R_real + this->D_est_neg[36] * Z_real;
			y_est[0] += buff;
			buff = 0.0;
			buff =  this->D_est_neg[25] * R_real + this->D_est_neg[37] * Z_real;
			y_est[1] += buff;
			buff = 0.0;
		
			

		
			//x(k)=x(k+1) for the next step

		
		Outputs.Kalman_R= y_est[0];
		Outputs.Kalman_Z= y_est[1];
		Outputs.X_est=this->x_neg;
		
		buffer = 0.0;
		for (i = 0; i < N_state; i++) {
			buffer =this->x_dot_neg[i];
			this->x_neg[i]=buffer;
		}
		
		}else{
			Outputs.Kalman_R= 0.085;
			Outputs.Kalman_Z= 0.085;
				for (i = 0; i < N_state; i++) {
						X_est[i]=0.0;}
										
			Outputs.X_est=X_est;
			}
	

	return Outputs;
	
}
int LQR:: erase() {
	int i = 0;
	float* state0_pos;
	float* state0_neg;
	state0_pos = (float[N_state]){0.0983539734804260, 0.480564466137287, -0.204246862173158, 0.0605647026605041, -0.369796209977501, 0.111380460346364, 0.0194246756583234, 0.144179620264490, -0.0534417070457820, 0.215640033094750}; 
	state0_neg=(float[N_state])  {-0.0494766508450061, -0.111034664506845, 0.0242441680811973, 0.00586343293525677, -0.00854660070320271, 0.000970048477408260, -0.000888992834197872, 0.00166157894650906, 0.000997825822154208, 0.00231148829292225};
	
	for (i = 0; i < N_state; i++) {
			this->x_dot_pos[i]=0.0;
		}
	for (i = 0; i < N_state; i++) {
			this->x_dot_neg[i]=0.0;
		}
	for (i = 0; i < N_state; i++) {
			this->x_pos[i]=state0_pos[i];
		}
	for (i = 0; i < N_state; i++) {
			this->x_neg[i]=state0_neg[i];
		}

	
	return 1;
	}


	/* POSITIVE
	/////normal
	this-> N_BAR_pos = (float[N_input +  N_output]){-16738.8719868375, 3479.82442345434,
							488.651395096507, -28.3958381254954};	
	this-> K_LQR_pos =(float[N_input * N_state]) {-39.8874724584811, -70.6769666120137, 354.869808121735, 370.159896447712, -106.569303855515, -10.6727031034049, -29.8250480200422, 75.3583282160950, -135.655340474040, -27.8586258177635,
							4.85034537803885, 16.3699444462248, -103.189045436383, -88.9493312103841, 36.1506690550597, 19.0337453517795, 10.4943236824918, -28.3848265454111, 45.6326833691854, 8.77044127460252};
	/////slow
	this-> N_BAR_pos = (float[N_input +  N_output]){-11086.5128004769, 2304.91025662465,
							-624.913612903589, 203.085650822326};
	this-> K_LQR_pos =(float[N_input * N_state]) {-8.01681293888861, -26.8325522087463, 134.585605910785, 71.3907687002046, -31.1613639245902, 1.06653823136754, -8.67641023817368, 30.6278817417476, -52.6867168563030, -17.6236782619040,
							-0.359505258997209, 5.27976739349132, -34.7284975804598, -24.9860321020505, 10.5336478520904, 3.20301260401835, 2.98116013872604, -9.62426558998255, 15.6891220336232, 4.25997674316981};
	
	
	//more slow
	this-> N_BAR_pos = (float[N_input +  N_output]){-10292.2693526553, 2133.74129809333,
							-944.140351815484, 269.810568562941};
	this-> K_LQR_pos =(float[N_input * N_state]) {-0.677239930551175, -5.56445051715661, 29.1530400398665, 10.5689835482429, -6.43042836144154, 0.784006911127944, -1.67848331436731, 6.89611601694075, -11.6198347798690, -4.41024916164279,
							-0.253637909222695, 0.998896201615173, -6.91291373508615, -4.34730845911924, 2.02654430640833, 0.350198601401546, 0.540286005706202, -1.95918796444674, 3.18376611395890, 0.992300688479676};
	
	*/
	/* NEGATIVE
	//////normal
	this-> K_LQR_neg=(float[N_input * N_state])  {-1.10243700688381, -26.6888168111469, 248.208386299799, 2272.69314089294, 341.435985384170, 116.974185358079, -1571.00753489878, -331.351679250269, -1271.06969797964, -4441.09340489016,
							-0.233414809301534, -21.1374545809135, 237.120488501101, 3170.68925311282, -143.026532612836, 258.466914847382, -1451.21045848359, -298.057584698670, -925.765520002125, -4247.87411829111};
	this-> N_BAR_neg=(float[N_input +  N_output]) {6692.27395164254, -960.635617013133,
							-2256.41280370701, -489.068332481912};
		
	/////slow
	this-> N_BAR_pos = (float[N_input +  N_output]){6589.72680834879, -954.006515986846,
							-2507.09039980905, -438.508741661756};
	this-> K_LQR_pos =(float[N_input * N_state]) {-0.140009763339089, -3.18557215603840, 30.5815322306256, 261.182404835037, 41.0701964065595, 15.1457297490992, -184.071391002032, -37.9502637307121, -150.255931786410, -524.709987410739,
							-0.0430382829459303, -2.40364558484570, 28.2516724393487, 349.930303853398, -19.3580296640594, 27.8929619552656, -166.244389576828, -34.2340029306290, -102.365591295273, -474.430220412449};
	
	
	//more slow
	this-> N_BAR_pos = (float[N_input +  N_output]){6577.20274262693, -953.521247320254,
							-2537.02933847542, -432.837113235285};
	this-> K_LQR_pos =(float[N_input * N_state]) {-0.0143959540956107, -0.325150394017629, 3.13416013810849, 26.5396916000373, 4.19324806144036, 1.55806909562497, -18.7468162690424, -3.85545880019961, -15.3160343140834, -53.4859040919299,
							-0.00458150757613837, -0.243881430718078, 2.88347865234915, 35.3882054725928, -2.01002787695297, 2.81135503841810, -16.8859632308927, -3.47894536276546, -10.3473795912101, -48.0323686610169};
	
	*/
