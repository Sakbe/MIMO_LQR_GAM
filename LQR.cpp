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
	
	this-> x_dot_pos = (float[N_state]){0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	
	this-> x_pos = (float[N_state]){0.0983539734804260, 0.480564466137287, -0.204246862173158, 0.0605647026605041, -0.369796209977501, 0.111380460346364, 0.0194246756583234, 0.144179620264490, -0.0534417070457820, 0.215640033094750}; 
	this-> K_LQR_pos =(float[N_input * N_state]) {-88.5447855580499, 72.5449707404182, 299.919285932049, 97.3986210388653, -126.922485598113, -19.2782249920040, -72.3603047217459, -35.5046768020497, -23.4791895209002, 34.4552016998678,
											    47.4794910727297, -62.4889455506822, -139.015357704676, -5.27359167725381, 55.1119649348921, 73.2619061481024, 64.3020177469677, 40.0122247290657, 12.5025905939338, -48.7322901014367};
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
	this-> N_BAR_pos = (float[4]){-9541.53743579577, 1313.89696343387
								  -2013.80067654175, 940.107633664446};
	this-> error_pos =(float[N_output]) {0, 0};
								  
	this-> x_neg=(float[N_state])  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	this-> x_dot_neg=(float[N_state])  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	this-> K_LQR_neg=(float[N_state])  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	this-> A_est_neg=(float[N_state])  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	this-> B_est_neg=(float[N_state]) {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	this-> N_BAR_neg=(float[N_state]) {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	this-> error_neg=(float[N_output]) {0, 0};
	

	
}

LQR::~LQR(){
}

float*  LQR::MIMO_CONTROL_POSITIVE(float R_ref, float Z_ref, float R_real, float Z_real, float I_vertical, float I_horizontal){
	float u_vertical = 0;
	float u_horizontal = 0;
	int i = 0;
	int j = 0;
	float temp = 0;
	float *outputs;
	outputs =(float[2]) {0, 0};
	this->error_pos[0] = R_ref - R_real;
	this->error_pos[1] = Z_ref - Z_real;
	//Nbar(4x4)*error(2x1)
	u_vertical = (this->error_pos[0] * this->N_BAR_pos[0]) + (this->error_pos[0] * this->N_BAR_pos[2]);
	u_horizontal = (this->error_pos[1] * this->N_BAR_pos[1]) + (this->error_pos[1] * this->N_BAR_pos[3]);
	
	//u=u-k_lqr*x
	for (i = 0; i < N_state; i++) {
		temp = this->x_pos[i] * this->K_LQR_pos[i];
		u_vertical = u_vertical - temp;
		temp = 0;		
	}
	for (i = 0; i < N_state; i++) {
		temp = this->x_pos[i] * this->K_LQR_pos[i+N_state];
		u_horizontal = u_horizontal - temp;
		temp = 0;		
	}
	
	//x(k+1)=Ax(k)
	for(j = 0; j < N_state; j++){
		for (i = 0; i < N_state; i++) {
			this->x_dot_pos[j] = this->A_est_pos[i + j*(N_state-1)] * this->x_pos[i];
		}
	}
	//B_est_pos =[vertical; horizonal; rc; zc]
	for(i = 0; i < N_state; i++){
		temp =  this->B_est_pos[i] * I_vertical + this->B_est_pos[ i + (N_state)] * I_horizontal + this->B_est_pos[ i + 2*(N_state)] * R_real + this->B_est_pos[ i + 3*(N_state)] * Z_real;
		this->x_dot_pos[i] = this->x_dot_pos[i] + temp;
		temp = 0;
	}
	//x(k)=x(k+1) for the next step
	for (i = 0; i < N_state; i++) {
		this->x_pos[i] = this->x_dot_pos[i];
	}
	outputs[0] = u_vertical;
	outputs[1] = u_horizontal;
	return outputs;
}

Kalman LQR:: KALMAN_FILTER(float R_real, float Z_real, float I_vertical, float I_horizontal, int sign){
	//the same as before
	//if we use this in mimo class we can comment from 102 to 115
	//i putted because we can copy this in IPID class
	float temp;
	int j=0;
	int i=0;
	Kalman Outputs={0.085,0.085};
	if(sign == 1){
	
	for(j = 0; j < N_state; j++){
		for (i = 0; i < N_state; i++) {
			this->x_dot_pos[j] = this->x_dot_pos[j]+this->A_est_pos[i + j*(N_state-1)] * x_pos[i];
		}
	}
	for(i = 0; i < N_state; i++){
		temp =  this->B_est_pos[i] * I_vertical + this->B_est_pos[ i + (N_state)] * I_horizontal + this->B_est_pos[ i + 2*(N_state)] * R_real + this->B_est_pos[ i + 3*(N_state)] * Z_real;
		x_dot_pos[i] = x_dot_pos[i] + temp;
		temp = 0;
	}

	//x^(k+1) = A_est * x(k)+ B_est * u(k)
	//i need also to copy her C_est and D_est to have y^(k) that we can compare with the real y
	float* y_est;
	//float* x_est;
	y_est =(float[2]) {0, 0};
	
	//[y_est x_est] = C_est x_est + D_est [u_real y_real]
	for(j = 0; j < N_output; j++){
		for (i = 0; i < N_state; i++) {
			y_est[j] = y_est[j] + this->C_est_pos[i + j*(N_output-1)] * this->x_pos[i];
		}}
	for (i = 0; i < N_output; i++) {
		y_est[i] = y_est[i] + this->D_est_pos[24 + i] * I_vertical + this->D_est_pos[36 + i] * I_horizontal;
	}
	/*//x_est
	for(j = this->N_output; j < this->N_output + this->N_state; j++){
		for (i = 0; i < this->N_state; i++) {
			this->x_est[j] = this->x_est[j] + this->C_est_pos[i + j*(this->N_output-1)] * x_pos[i];
		}}
	//x_est = D*u
	//float I_vertical, float I_horizontal
	for (i = 0; i < this->N_state; i++) {
		this->x_pos[i] = this->x_pos[i] + this->D_est_pos[26 + i] * I_vertical + this->D_est_pos[38 + i] * I_horizontal;
	}*/
	
		//x(k)=x(k+1) for the next step
	for (i = 0; i < N_state; i++) {
		this->x_pos[i] = this->x_dot_pos[i];
	}
				Outputs.Kalman_R=y_est[0];
				Outputs.Kalman_R=y_est[1];
		}
		
		else{
	
				Outputs.Kalman_R=0.085;
				Outputs.Kalman_Z=0.085;}
	

	return Outputs;
	
}

int LQR:: nothing() {
	return 1;
	}
