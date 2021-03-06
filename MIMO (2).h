//******************************************************************************
//
//		MIMO.h - LQR CONTROL AND KALMAN FILTER
//		CORONA CIANCIULLI 05/2019
//
//******************************************************************************

#if !defined (MIMO_H)
#define MIMO_H


#include "Level0.h"
#include "Level1.h"
#include "Level2.h"
//OBJECT_DLL(IPID)

/** MIMO CONTROLLER */
class MIMO {

//OBJECT_DLL_STUFF(IPID)
private:
	float N_state;
	float N_input;
	float N_output;
	
	float *x_dot_pos;
	float *x_pos;
	float *K_LQR_pos;
	float *A_est_pos;
	float *B_est_pos;
	float *C_est_pos;
	float *D_est_pos;
	float *N_BAR_pos;
	float *error_pos;
	
	float *x_dot_neg;
	float *x_neg;
	float *K_LQR_neg;
	float *A_est_neg;
	float *B_est_neg;
	float *C_est_neg;
	float *D_est_neg;
	float *N_BAR_neg;
	float *error_neg;	

public:	
	MIMO ();
	float* MIMO_CONTROL_POSITIVE(float R_ref, float Z_ref, float R_real, float Z_real, float I_vertical, float I_horizontal);
	bool MIMO_CONTROL_NEGATIVE(float R_ref, float Z_ref, float R_real, float Z_real, float I_vertical, float I_horizontal);
	float KALMAN_FILTER(float R_real, float Z_real, float I_vertical, float I_horizontal);
	
	~MIMO();

private:

	bool SortWaveform();
	bool RemoveRepeatedValues();

};


#endif











