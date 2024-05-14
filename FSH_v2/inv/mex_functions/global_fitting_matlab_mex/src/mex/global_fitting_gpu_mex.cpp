// global_fitting_gpu_mex.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>

// local_fitting_cuda_mex.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include "mex.h"
#include <iostream>
#include <complex>
using namespace std;


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

//require the range compressed data arranged in (Ntrack,Nazimuth,Nrange)

#define My_PI  3.14159265358979323846
#define max(a,b) (a)>(b)?(a):(b)
#define min(a,b) (a)<(b)?(a):(b)

inline double rad2deg(double radian)
{
	return radian * 180 / My_PI;
}

inline double deg2rad(double deg)
{
	return deg * My_PI / 180;
}
inline double sqr(const double &x)
{
	return x * x;
}


extern "C" int SearchKb(const float * coh, const float * chm, int nObs, float S_min, float S_max, int nS, float C_min, float C_max, int nC,
	const float * sinc_LUT_xx, const float * sinc_LUT_yy, int nLut, double *return_val, double height_upper);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	if (nrhs != 11)
	{
		mexErrMsgTxt("11 Inputs Required!");
		return;
	}

	//reading input parameters

    //coherence input 
	if (!mxIsSingle(prhs[0]))
	{
		mexErrMsgTxt("The 1th input parameters must be a single array!");
		return;
	}


	// for matlab version >= 2017
	float * coh_vec = (float *)mxGetPr(prhs[0]);

	int nobs = max(mxGetM(prhs[0]), mxGetN(prhs[0]));


	// change


	if (!mxIsSingle(prhs[1]))
	{
		mexErrMsgTxt("The second input parameters must be a float array!");
		return;
	}

	float * chm_vec = (float *)mxGetPr(prhs[1]);

	if (!mxIsScalar(prhs[2]))
	{
		mexErrMsgTxt("The third input parameters must be a scalar!");
		return;
	}
	float Smin = mxGetScalar(prhs[2]);


	if (!mxIsScalar(prhs[3]))
	{
		mexErrMsgTxt("The 4th input parameters must be a scalar!");
		return;
	}

	float Smax = mxGetScalar(prhs[3]);


	if (!mxIsScalar(prhs[4]))
	{
		mexErrMsgTxt("The 5th input parameters must be a scalar!");
		return;
	}

	int nS = mxGetScalar(prhs[4]);

	if (!mxIsScalar(prhs[5]))
	{
		mexErrMsgTxt("The 6th input parameters must be a scalar!");
		return;
	}

	float Cmin = mxGetScalar(prhs[5]);


	if (!mxIsScalar(prhs[6]))
	{
		mexErrMsgTxt("The 7th input parameters must be a scalar!");
		return;
	}

	float Cmax = mxGetScalar(prhs[6]);

	if (!mxIsScalar(prhs[7]))
	{
		mexErrMsgTxt("The 8th input parameters must be a scalar!");
		return;
	}

	int nC = mxGetScalar(prhs[7]);


	if (!mxIsSingle(prhs[8]))
	{
		mexErrMsgTxt("The 9th input parameters must be single vector!");
		return;
	}

	float * lut_xx = (float *)mxGetPr(prhs[8]);
	int nlut = max(mxGetM(prhs[8]), mxGetN(prhs[8]));

	if (!mxIsSingle(prhs[9]))
	{
		mexErrMsgTxt("The 10th input parameters must be single vector!");
		return;
	}

	float * lut_yy = (float *)mxGetPr(prhs[9]);
	if (!mxIsScalar(prhs[10]))
	{
		mexErrMsgTxt("The 11th input parameters must be a scalar!");
		return;
	}
	double height_upper = mxGetScalar(prhs[10]);



	//double * h_resi = new double[nvalid_gedi];
	//double * S_est = new double[nvalid_gedi];
	//double * C_est = new double[nvalid_gedi];

	//memset(h_resi, 0, nvalid_gedi * sizeof(double));
	//memset(S_est, 0, nvalid_gedi * sizeof(double));
	//memset(C_est, 0, nvalid_gedi * sizeof(double));

	//int * Lines_central = new int[nvalid_gedi];
	//int * Pixel_central = new int[nvalid_gedi];

	//memset(Lines_central, 0, nvalid_gedi * sizeof(int));
	//memset(Pixel_central, 0, nvalid_gedi * sizeof(int));







	/**********************************************
	Define the output array
	**********************************************/
	if (nlhs != 2)
	{
		mexErrMsgTxt("Only six output parameter!");
	}



	double S_global = 0.0f;
	double C_global = 0.0f;


	






	//plhs[0] = mxCreateDoubleMatrix(nFocRow, nFocCol, mxCOMPLEX);





	//test reading array 
 //   mexPrintf("coh_input(500,500): %lf\n", coh_input[499 * col_ + 499]);
	//mexPrintf("chm_input(500,500): %lf\n", chm_input[499 * col_ + 499]);
	//mexPrintf("row_:%d\n", row_);
	//mexPrintf("col_:%d\n", col_);
	//mexPrintf("nvalid_gedi:%d\n", nvalid_gedi);
	//mexPrintf("S_min:%f\n", S_min);
	//mexPrintf("S_max:%f\n", S_max);
	//mexPrintf("nS:%d\n", nS);
	//mexPrintf("C_min:%f\n", C_min);
	//mexPrintf("C_max:%f\n", C_max);
	//mexPrintf("nC:%d\n", nC);
	//mexPrintf("sinc_LUT_xx[717]:%f\n", sinc_LUT_xx[716]);
	//mexPrintf("sinc_LUT_yy[717]:%f\n", sinc_LUT_yy[716]);
	//mexPrintf("nLut:%d\n", nLut);
	//mexPrintf("nwinx:%d\n", nwinx);
	//mexPrintf("nwiny:%d\n", nwiny);
	//mexPrintf("weight_input(central, central):%f\n", weight_input[16 * (nwinx+1)+ 16]);
	//mexPrintf("chm_valid_max:%f\n", chm_valid_max);
	//mexPrintf("chm_valid_min:%f\n", chm_valid_min);
	//mexPrintf("coh_valid:%f\n", coh_valid);

	//return;
	 
	//mexPrintf("OvsFac:%lf\n", OvsFactor);
	//

	//mexPrintf("coh_vec[500]:%lf\n", coh_vec[499]);
	//mexPrintf("chm_vec[500]:%lf\n", chm_vec[499]);
	//mexPrintf("nobs:%d\n", nobs);
	//mexPrintf("Smin:%lf\n", Smin);
	//mexPrintf("Smax:%lf\n", Smax);
	//mexPrintf("nS:%d\n", nS);
	//mexPrintf("Cmin:%lf\n", Cmin);
	//mexPrintf("Cmax:%lf\n", Cmax);
	//mexPrintf("nC:%d\n", nC);
	//mexPrintf("lut_xx[500]:%lf\n", lut_xx[499]);
	//mexPrintf("lut_yy[500]:%lf\n", lut_yy[499]);
	//mexPrintf("nlut:%d\n", nlut);
	//mexPrintf("h eight_upper:%f\n", height_upper);
	    
	int valid_num = 0.0f;


	double return_val[2] = { 0.0 };
	float duration = 0.f;
	duration = SearchKb(coh_vec, chm_vec, nobs, Smin, Smax, nS, Cmin, Cmax, nC,
		lut_xx, lut_yy, nlut, return_val, height_upper);



	if (duration < 0)
	{
		mexPrintf("error code:%d\n",-duration);
		mexErrMsgTxt("A non-zero code is returned during the GPU computing. Please debug it!");

	}
	else
	{
		mexPrintf("GPU computation time: %f ms.\n",
			duration);
	}

	plhs[0] = mxCreateDoubleScalar(return_val[0]);
	plhs[1] = mxCreateDoubleScalar(return_val[1]);

	return;

}











