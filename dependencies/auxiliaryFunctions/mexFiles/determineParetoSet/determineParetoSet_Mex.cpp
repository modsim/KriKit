#include <vector>
#include "determineParetoPoints.h"
#include "mex.h"   
using namespace std;


 
  
void mexFunction(int nlhs, mxArray *plhs[],int nrhs,const mxArray *prhs[])  
{  
	
	if (nrhs==0){
		printf("[paretoSet]=determineParetoSet(objectiveMatrix)\n");
		printf("objectiveMatrix ... Matrix of size nSamples X nObjectives\n");
		return;
	}
	// Input declaration
		// Number of Points in the grid for testing
	double * objectiveMatrix = (double*)mxGetData(prhs[0]);

    // Inititalization of characterstic numbers
	size_t nSamples = (int)mxGetM(prhs[0]);
	size_t nObjectives = (int)mxGetN(prhs[0]);
	vector<bool> paretoOptimalBool(nSamples, false);

	// If Input
	if (nSamples == 0){ printf("Matrix is empty\n"); return; }
	if (nObjectives == 1){ printf("Only One Objective ... Use rather min or max\n"); return; }
	

	//printf("nSamples: %i - nObjectives:%i \n", nSamples, nObjectives);
	class determineParetoPoints determinationObj(nSamples, nObjectives, objectiveMatrix);

	determinationObj.doCalculation(paretoOptimalBool);
	size_t nNonDominatedPoints = accumulate(paretoOptimalBool.begin(), paretoOptimalBool.end(), 0);
	
   // Define output
	plhs[0] = mxCreateNumericMatrix(nNonDominatedPoints, nObjectives, mxDOUBLE_CLASS, mxREAL);
		// Save Output values
   double* output = (double*)mxGetData(plhs[0]);
   size_t iParetoPoint = 0;
   for (size_t iSample = 0; iSample < nSamples; iSample++){
	   // Only save non-dominated samples
	   if (paretoOptimalBool[iSample]){
		   for (size_t iObective = 0; iObective < nObjectives; iObective++){
			   output[iParetoPoint + iObective*nNonDominatedPoints] = *(objectiveMatrix + iSample + iObective*nSamples);
		   }
		   // Updtae Running Index
		   iParetoPoint++;
	   }
	   
   }

   //return;
}  