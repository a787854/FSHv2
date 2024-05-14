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


float Local_search(const float * coh, const float * chm, int Lines, int Pixels, int nvalid_gedi, float S_min, float S_max, int nS, float C_min, float C_max, int nC,
	const float * sinc_LUT_xx, const float * sinc_LUT_yy, int nLut, int nwinx, int nwiny, float * nweight, double * h_resi, double * S_est, double * C_est, int * Lines_central,
	int *Pixel_central, double chm_valid_max, double chm_valid_min, double coh_valid, int n_least_valid, int &n_valid);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	if (nrhs != 18)
	{
		mexErrMsgTxt("18 Inputs Required!");
		return;
	}

	//reading input parameters

	//reading Scene in Complex number 
	if (!mxIsSingle(prhs[0]))
	{
		mexErrMsgTxt("The 1th input parameters must be a single array!");
		return;
	}



	// for matlab version >= 2017
	float * coh_mat = (float *)mxGetPr(prhs[0]);

	int row_ = mxGetM(prhs[0]);
	int col_  = mxGetN(prhs[0]);


	// change 

	float * coh_input = new float[row_ * col_];
	for (int i = 0; i < row_; i++)
	{
		for (int j = 0; j < col_; j++)
		{
			
			coh_input[i * col_ + j] = coh_mat[j* row_ + i];
		}
	}
	


	if (!mxIsSingle(prhs[1]))
	{
		mexErrMsgTxt("The second input parameters must be a float array!");
		return;
	}

	float * chm_mat = (float *)mxGetPr(prhs[1]);


	if (mxGetM(prhs[1]) != row_ || mxGetN(prhs[1]) != col_)
	{
		mexErrMsgTxt("The second array should be in the same shape with the first array!");
		return;
	}


	float * chm_input = new float[row_* col_];
	for (int i = 0; i < row_; i++)
	{
		for (int j = 0; j < col_; j++)
		{

			chm_input[i * col_ + j] = chm_mat[j* row_ + i];
		}
	}

	if (!mxIsScalar(prhs[2]))
	{
		mexErrMsgTxt("The third input parameters must be a scalar!");
		return;
	}
	int nvalid_gedi = mxGetScalar(prhs[2]);


	if (!mxIsScalar(prhs[3]))
	{
		mexErrMsgTxt("The 4th input parameters must be a scalar!");
		return;
	}

	float S_min = mxGetScalar(prhs[3]);

	
	if (!mxIsScalar(prhs[4]))
	{
		mexErrMsgTxt("The 5th input parameters must be a scalar!");
		return;
	}

	float S_max = mxGetScalar(prhs[4]);

	if (!mxIsScalar(prhs[5]))
	{
		mexErrMsgTxt("The 6th input parameters must be a scalar!");
		return;
	}

	int nS = mxGetScalar(prhs[5]);


	if (!mxIsScalar(prhs[6]))
	{
		mexErrMsgTxt("The 7th input parameters must be a scalar!");
		return;
	}

	float C_min = mxGetScalar(prhs[6]);

	if (!mxIsScalar(prhs[7]))
	{
		mexErrMsgTxt("The 8th input parameters must be a scalar!");
		return;
	}

	float C_max = mxGetScalar(prhs[7]);

	
	if (!mxIsScalar(prhs[8]))
	{
		mexErrMsgTxt("The 9th input parameters must be a scalar!");
		return;
	}

	int nC = mxGetScalar(prhs[8]);

	if (!mxIsSingle(prhs[9]))
	{
		mexErrMsgTxt("The 10th input parameters must be a float array!");
		return;
	}

	if (mxGetM(prhs[9]) != 1 && mxGetN(prhs[9]) != 1)
	{
		mexErrMsgTxt("The 10th input parameters must be a float vector!");
		return;
	}

	int nLut = max(mxGetM(prhs[9]), mxGetN(prhs[9]));


	float *sinc_LUT_xx = (float*)mxGetPr(prhs[9]);

	if (!mxIsSingle(prhs[10]))
	{
		mexErrMsgTxt("The 11th input parameters must be a float array!");
		return;
	}

	float *sinc_LUT_yy = (float*)mxGetPr(prhs[10]);


	if (!mxIsScalar(prhs[11]))
	{
		mexErrMsgTxt("The 12th input parameters must be a scalar!");
		return;
	}

	int nwinx = mxGetScalar(prhs[11]);

	if (!mxIsScalar(prhs[12]))
	{
		mexErrMsgTxt("The 13th input parameters must be a scalar!");
		return;
	}

	int nwiny = mxGetScalar(prhs[12]);

	
	if (!mxIsSingle(prhs[13]))
	{
		mexErrMsgTxt("The 14th input parameters must be a single matrix!");
		return;
	}

	float * weight_local = (float*) mxGetPr(prhs[13]);


	float * weight_input = new float[(nwinx+1) * (nwiny+1)];

	for (int i = 0; i < nwiny + 1; i++)
	{
		for (int j = 0; j < nwinx + 1; j++)
		{
			weight_input[i* (nwinx + 1) + j] = weight_local[j* (nwiny + 1) + i];
			//mexPrintf("weight_local(%d, %d):%f\n", i,j,weight_input[i* (nwinx + 1) + j]);
		}

	}

	if (!mxIsScalar(prhs[14]))
	{
		mexErrMsgTxt("The 15th input parameters must be a scalar!");
		return;
	}

	float chm_valid_max = mxGetScalar(prhs[14]);

	if (!mxIsScalar(prhs[15]))
	{
		mexErrMsgTxt("The 16th input parameters must be a scalar!");
		return;
	}

	float chm_valid_min = mxGetScalar(prhs[15]);


	if (!mxIsScalar(prhs[16]))
	{
		mexErrMsgTxt("The 17th input parameters must be a scalar!");
		return;
	}

	float coh_valid = mxGetScalar(prhs[16]);


	if (!mxIsScalar(prhs[17]))
	{
		mexErrMsgTxt("The 18th input parameters must be a scalar!");
		return;
	}
	int n_least_valid = mxGetScalar(prhs[17]);


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


	int ntotal_num = 0;
	
	




	/**********************************************
	Define the output array
	**********************************************/
	if (nlhs != 6)
	{
		mexErrMsgTxt("Only six output parameter!");
	}


	const size_t nD = 1;
	size_t Dimension[nD] = { nvalid_gedi };


	plhs[0] = mxCreateNumericArray(nD, Dimension, mxDOUBLE_CLASS, mxREAL);


	//double * output_r = mxGetPr(plhs[0]);
	//double * output_i = mxGetPi(plhs[0]);
	double * h_resi = mxGetPr(plhs[0]);

	//plhs[0] = mxCreateDoubleMatrix(nFocRow, nFocCol, mxCOMPLEX);

	//complex<float>* output = new complex<float>[nfocx*nfocy*nfocz];

	plhs[1] = mxCreateNumericArray(nD, Dimension, mxDOUBLE_CLASS, mxREAL);
	double * S_par = mxGetPr(plhs[1]);

	plhs[2] = mxCreateNumericArray(nD, Dimension, mxDOUBLE_CLASS, mxREAL);
	double * C_par = mxGetPr(plhs[2]);


	plhs[3] = mxCreateNumericArray(nD, Dimension, mxINT32_CLASS, mxREAL);
	int * Lines_central = (int*)mxGetPr(plhs[3]);


	plhs[4] = mxCreateNumericArray(nD, Dimension, mxINT32_CLASS, mxREAL);
	int * Pixels_central = (int*)mxGetPr(plhs[4]);



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

	int valid_num = 0;

	float duration = Local_search(coh_input, chm_input, row_, col_, nvalid_gedi, S_min, S_max, nS, C_min, C_max, nC,
		sinc_LUT_xx, sinc_LUT_yy, nLut, nwinx, nwiny, weight_input, h_resi, S_par, C_par, Lines_central,
		Pixels_central, chm_valid_max, chm_valid_min, coh_valid, n_least_valid, valid_num);

	plhs[5] = mxCreateDoubleScalar(valid_num);

	if (duration < 0)
	{
		mexErrMsgTxt("A non-zero code is returned during the GPU computing. Please debug it!");

	}
	else
	{
		mexPrintf("GPU computation time: %f s.\n",
			duration );
	}


	delete[] weight_input;
	delete[] coh_input;
	delete[] chm_input;

	return;

}










