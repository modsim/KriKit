#include "calculatePcoverOrPsigma.h"


calculatePcoverOrPsigma::calculatePcoverOrPsigma(size_t nRealizationsInput)
	:nMemberParetoCurve(nRealizationsInput, 1){
}

calculatePcoverOrPsigma::calculatePcoverOrPsigma(int nRealizationsInput, int nGridPointsInput, int nObjectsInput, bool calcPCoverAndNotPSigmaInput, int iSumOfMembers)
{
	// 	mexPrintf("calculatePcoverOrPsigma\n");

	calcPCoverAndNotPSigma = calcPCoverAndNotPSigmaInput;
	nGridPoints = nGridPointsInput;
	nRealizations = nRealizationsInput;
	nObjects = nObjectsInput;
	nMemberParetoCurve = vector<int>(nRealizationsInput, 0);
	sum_nMemberPareto = iSumOfMembers;
	nExpectedParetoPoints = 0;
	ExpectedParetoPoints = NULL;
	totalParetoSet = vector<vector<double>>(sum_nMemberPareto, vector<double>(nObjects, 0.0));
}

calculatePcoverOrPsigma::~calculatePcoverOrPsigma()
{
}

void calculatePcoverOrPsigma::doCalculation(double * gridX, double* outputVec){
	// Declaration
	vector<bool> dominatedBool(nRealizations, false);
	int max_nMemberPareto = *(max_element(nMemberParetoCurve.begin(), nMemberParetoCurve.end()));
	vector<vector<bool>> smallerThanEqual(max_nMemberPareto, vector<bool>(nObjects, false));
	vector<vector<bool>> smallerThanEqualTotal(sum_nMemberPareto, vector<bool>(nObjects, false));
	vector<vector<bool>> smallerThanEqualExpected(nExpectedParetoPoints, vector<bool>(nObjects, false));
	bool dominatedByExpected_Bool;
	bool dominatedByRealization_Bool;
	

	// Actual Calculation (determine value for each point in the grid)
	for (size_t iGridPoint = 0; iGridPoint < nGridPoints; iGridPoint++){
		// Initialization for each round
		// Matrix for saving if grid point is dominated or not
		for (size_t iRealization = 0; iRealization < nRealizations; iRealization++){
			dominatedBool[iRealization] = false;
		}
		// 
		for (size_t iRow = 0; iRow < smallerThanEqual.size(); iRow++){
			for (size_t iCol = 0; iCol < smallerThanEqual[0].size(); iCol++){
				smallerThanEqual[iRow][iCol] = false;
			}
		}

		// Check in case of 
		if ((!calcPCoverAndNotPSigma) && (*(pCoverArrayAPriori + iGridPoint) == 0 || *(pCoverArrayAPriori + iGridPoint) == 1)){
			*(outputVec + iGridPoint) = 0;
		}
		else{

			// Check if point is dominated by points of RNP(random non - dominated
			// points) associated with realizations (Do it once for all)
			checkIfDominatedByRNP(smallerThanEqualTotal, gridX, iGridPoint);

			// Check if point is dominated by points of expected Pareto front
			if (!calcPCoverAndNotPSigma){
				checkIfDominatedByExpectedCurve(smallerThanEqualExpected, gridX, iGridPoint);
			}

			// Check for each realization
			for (size_t iRealization = 0; iRealization < nRealizations; iRealization++){

				if (calcPCoverAndNotPSigma){
					// Find points which are dominated by realization
					dominatedBool[iRealization] = checkIfRealizationDominatesGridPoint(smallerThanEqualTotal, iRealization);
				}
				else{
					// Find points which are dominated by expected pareto curve
					dominatedByExpected_Bool = checkIfDominatedByAtLeatOneExpectedPoint(smallerThanEqualExpected);
					// Find points which are dominated by realization
					dominatedByRealization_Bool = checkIfRealizationDominatesGridPoint(smallerThanEqualTotal, iRealization);


					if ((dominatedByRealization_Bool&&!dominatedByExpected_Bool) || (!dominatedByRealization_Bool&&dominatedByExpected_Bool)){
						dominatedBool[iRealization] = true;
					}
					else{
						dominatedBool[iRealization] = false;
					}
				}

			}

			// Final Value = Ratio of realizations
			*(outputVec + iGridPoint) = accumulate(dominatedBool.begin(), dominatedBool.end(), 0) / ((double)nRealizations);
		}
	}
}

void calculatePcoverOrPsigma::checkIfDominatedByRNP(vector<vector<bool>> &smallerThanEqualTotal, double* gridPoints, size_t iGridPoint){

	for (size_t iRow = 0; iRow< smallerThanEqualTotal.size(); iRow++){
		for (size_t iCol = 0; iCol < smallerThanEqualTotal[0].size(); iCol++){
			smallerThanEqualTotal[iRow][iCol] = totalParetoSet[iRow][iCol] <= *(gridPoints + iGridPoint + iCol*nGridPoints);
		}
	}
}

void calculatePcoverOrPsigma::checkIfDominatedByExpectedCurve(vector<vector<bool>> &smallerThanEqualExpected, double* gridPoints, size_t iGridPoint){
    
	size_t nRows = smallerThanEqualExpected.size();
	for (size_t iRow = 0; iRow< nRows; iRow++){
		for (size_t iCol = 0; iCol < smallerThanEqualExpected[0].size(); iCol++){
			smallerThanEqualExpected[iRow][iCol] = *(ExpectedParetoPoints + iRow + iCol*nRows) <= *(gridPoints + iGridPoint + iCol*nGridPoints);
		}
	}
}

bool calculatePcoverOrPsigma::checkIfRealizationDominatesGridPoint(vector<vector<bool>> const &smallerThanEqualTotal, size_t iRealization){
	int testValueObjectDominated = 0;
	size_t nConsideredMembersBeginning = 0;
	size_t nConsideredMembersEnd = sum_nMemberPareto;

	// Determin which entries are associated which the current realization
	if (iRealization > 0){
		nConsideredMembersBeginning = accumulate(nMemberParetoCurve.begin(), nMemberParetoCurve.begin() + iRealization, 0);
	}
	
	if (iRealization < nRealizations - 1){
		nConsideredMembersEnd = accumulate(nMemberParetoCurve.begin(), nMemberParetoCurve.begin() + iRealization + 1, 0);
	}

	bool checkResult = true;
	for (size_t iRow = nConsideredMembersBeginning; iRow < nConsideredMembersEnd; iRow++)
	{
		checkResult = true;
		//printf("iRealization: %i - iRow %i:  ", iRealization, iRow);
		for (size_t iObj = 0; iObj < nObjects; iObj++){
			checkResult = checkResult&smallerThanEqualTotal[iRow][iObj];
		}
		if (checkResult){ return true; }
	}
		
	return false;
}


bool calculatePcoverOrPsigma::checkIfDominatedByAtLeatOneExpectedPoint(vector<vector<bool>> const &smallerThanEqualExpected){

	bool checkResult = true;
	for (size_t iRow = 0; iRow < smallerThanEqualExpected.size(); iRow++)
	{
		checkResult = true;

		for (size_t iObj = 0; iObj < nObjects; iObj++){
			checkResult = checkResult&smallerThanEqualExpected[iRow][iObj];
		}
		if (checkResult){ return true; }
	}

	return false;
}
