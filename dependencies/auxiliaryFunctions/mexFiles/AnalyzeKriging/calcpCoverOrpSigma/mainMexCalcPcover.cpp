#include <vector>
#include "calculatePcoverOrPsigma.h"
#include "mex.h"   
using namespace std;


 
  
void mexFunction(int nlhs, mxArray *plhs[],int nrhs,const mxArray *prhs[])  
{  

	// Input Definition
	#define nMemberParetoCurveMexArray prhs[2];
	#define totalParetoSetMexArray prhs[4];
	#define gridMexArray prhs[5];

	// Input declaration
		// Number of Points in the grid for testing
	double nGridPoints = *((double*)mxGetData(prhs[0]));
		// Number of realization done by conditional simulation
	double nRealizations = *((double *)mxGetData(prhs[1]));
		// Number of pareto optimal points found in each realization
	double * nMemberParetoCurve = (double*)mxGetData(prhs[2]);
		// Number of objects of investigation (Pareto optimality)
	double nObjects = *((double*)mxGetData(prhs[3]));
		// Array containing the pareto optimal points of all realizations
	double* totalParetoSetArray = (double*)mxGetData(prhs[4]);
		// Decide if pCover or pSigma should be calculated
	bool* calcPCoverAndNotPSigma = (bool*)mxGetData(prhs[6]);
		// In case of pSigma: Save the points of expected pareto curve
	double * expectedParetoPoints = NULL;
		// In case of pSigma: Save the points of expected pareto curve
	double * pCoverArrayAPriori = NULL;
		// Number of points in the expected Pareto curve
    size_t nExpectedParetoPoints = 0;
		// COrrect initial values if needed
	if (!(*calcPCoverAndNotPSigma)){
		expectedParetoPoints = (double*)mxGetData(prhs[7]);
		pCoverArrayAPriori = (double*)mxGetData(prhs[8]);
		nExpectedParetoPoints = (int)mxGetM(prhs[7]);
	}
	
    // Inititalization of characterstic numbers
	size_t nRowsTotalParetoSet;
	size_t nColumnsTotalParetoSet;
	size_t nElementsMemberParetoCurve;

	// Declaration
	nRowsTotalParetoSet = (int)mxGetM(prhs[4]);
	nColumnsTotalParetoSet = (int)mxGetN(prhs[4]);

	// Count total number of pareto points over all relizations
	nElementsMemberParetoCurve = (int)mxGetM(prhs[2]);
	int iSum = 0;
	for (size_t iElement = 0; iElement < nElementsMemberParetoCurve; iElement++){
		iSum += *(nMemberParetoCurve + iElement);
	}

	// Declaration of calculation object
    class calculatePcoverOrPsigma calculationObj((int)nRealizations, (int)nGridPoints, (int)nObjects,
            *calcPCoverAndNotPSigma, iSum);
	
	
	// Correct missig values
		// Import for pSigma calculation
    calculationObj.nExpectedParetoPoints = nExpectedParetoPoints;
	calculationObj.ExpectedParetoPoints = expectedParetoPoints;
	calculationObj.pCoverArrayAPriori = pCoverArrayAPriori;
		// Save nMemberParetoCurve array
   for (size_t iElement = 0; iElement < nElementsMemberParetoCurve; iElement++){
	   calculationObj.nMemberParetoCurve[iElement] = *(nMemberParetoCurve + iElement);
   }
		// Save totalParetoSetArray
   for (size_t iRow = 0; iRow < nRowsTotalParetoSet; iRow++){
	   for (size_t iCol = 0; iCol < nColumnsTotalParetoSet; iCol++){
		   calculationObj.totalParetoSet[iRow][iCol] = *(totalParetoSetArray + (iCol*nRowsTotalParetoSet + iRow));
	   }
   }

   // Do Actual Calculation
   double * pCoverArray = (double *)mxCalloc(nGridPoints, sizeof(double));
   calculationObj.doCalculation((double*)mxGetData(prhs[5]), pCoverArray);


   // Define output
   plhs[0] = mxCreateNumericMatrix(nGridPoints, 1, mxDOUBLE_CLASS, mxREAL);
		// Save Output values
   double* output = (double*)mxGetData(plhs[0]);
   for (size_t iGridPoint = 0; iGridPoint < nGridPoints; iGridPoint++){
	   output[iGridPoint] = *(pCoverArray + iGridPoint);
   }

   return;
}  