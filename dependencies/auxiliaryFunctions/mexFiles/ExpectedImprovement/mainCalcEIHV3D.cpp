#include "mex.h"   
//#include "ehvi_calculations.h"
//#include "ehvi_sliceupdate.h"
//#include "ehvi_montecarlo.h"
#include "ehvi_multi.h"
#include <vector>
#include <cstring>
using namespace std;

//Checks if p dominates P. Removes points dominated by p from P and return the number of points removed.
int checkdominance(vector<individual> & P, individual p){
	int nr = 0;
	for (int i = P.size() - 1; i >= 0; i--){
		if (p.f[0] >= P[i].f[0] && p.f[1] >= P[i].f[1] && p.f[2] >= P[i].f[2]){
			printf("Individual %i is dominated or the same as another point; removing.\n", (i + 1));
			P.erase(P.begin() + i);
			nr++;
		}
	}
	return nr;
}


vector<double>  loadtestcase_interface(int d_r_exi, const double* data0, int d_r_eval, double* data2, double* data1, vector<individual> & testcase,												 vector<double> r){
	int inds = 0;
	for (int iRow = 0; iRow<d_r_exi; iRow++){

		testcase[iRow].f[0] = *(data0 + iRow);
		testcase[iRow].f[1] = *(data0 + iRow + d_r_exi);
		testcase[iRow].f[2] = *(data0 + iRow + 2 * d_r_exi);

		//checkdominance(testcase, testcase[iRow]);
	}

	//r[0] = *(data1);
	//r[1] = *(data1 + 1);
	//r[2] = *(data1 + 2);
	/*printf("Ref: ");*/

	//vector<mus> pdf(1, new mus);
	mus tempMus;
	vector<mus> pdf(d_r_eval, tempMus);
	for (int iReference = 0; iReference<d_r_eval; iReference++){
		
		pdf[iReference].mu[0] = *(data2 + iReference);
		pdf[iReference].mu[1] = *(data2 + iReference + 1 * d_r_eval);
		pdf[iReference].mu[2] = *(data2 + iReference + 2 * d_r_eval);
		//printf("pdf: %i\n", iReference);

		pdf[iReference].s[0] = *(data2 + iReference + 3 * d_r_eval);
		pdf[iReference].s[1] = *(data2 + iReference + 4 * d_r_eval);
		pdf[iReference].s[2] = *(data2 + iReference + 5 * d_r_eval);
		//printf("pdf: %i\n", iReference);
	}

	
	vector<double> ehvi_2d = ehvi3d_sliceupdate(testcase, r, pdf);

	return ehvi_2d;
}
  
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{  

	int d_r_exi, d_c_exi, d_r_eval, d_c_eval, d_r_r, d_c_r;
	
	vector<double> r(DIMENSIONS,0);
	double mu[DIMENSIONS];
	double s[DIMENSIONS];

	const  double* data0 = (const double*)mxGetData(prhs[0]);
	const  double* dataN = (const double*)mxGetData(prhs[0]);
	d_r_exi = mxGetM(prhs[0]);
	d_c_exi = mxGetN(prhs[0]);

	d_r_r = mxGetM(prhs[1]);
	d_c_r = mxGetN(prhs[1]);
	double * data1 = (double*)mxGetData(prhs[1]);
	r[0] = *(data1);
	r[1] = *(data1 + 1);
	r[2] = *(data1 + 2);

	vector<individual> testcase(d_r_exi, individual());

	// Evaluation points
	d_r_eval = mxGetM(prhs[2]);
	d_c_eval = mxGetN(prhs[2]);
	double * data2 = (double*)mxGetData(prhs[2]);

	vector<double> answer = loadtestcase_interface(d_r_exi, dataN, d_r_eval, data2, data1, testcase, r);


	// Define output
	plhs[0] = mxCreateNumericMatrix(d_r_eval, 1, mxDOUBLE_CLASS, mxREAL);
	// Save Output values
	double* outputBuffer1 = (double*)mxGetData(plhs[0]);
	for (int i = 0; i<d_r_eval; i++){
		outputBuffer1[i] = answer[i];
	}
 
}  
