#pragma once
#include "mex.h"   
#include <cstdio>
#include <vector>
#include <algorithm>
#include <numeric>

using namespace std;

class calculatePcoverOrPsigma
{
public:
	calculatePcoverOrPsigma();
	calculatePcoverOrPsigma(size_t nRealizationsInput);
	calculatePcoverOrPsigma(int nRealizationsInput, int nGridPointsInput, int nObjectsInput, bool calcPCoverAndNotPSigma, int iSumOfMembers);
	~calculatePcoverOrPsigma();


	void doCalculation(double * gridX, double* outputVec);

	void checkIfDominatedByRNP(vector<vector<bool>> &smallerThanEqualTotal, double* gridPoints, size_t iGridPoint);
	void checkIfDominatedByExpectedCurve(vector<vector<bool>> &smallerThanEqualExpected, double* gridPoints, size_t iGridPoint);
	bool checkIfRealizationDominatesGridPoint(vector<vector<bool>> const &smallerThanEqualTotal, size_t iRealization);
	bool checkIfDominatedByAtLeatOneExpectedPoint(vector<vector<bool>> const &smallerThanEqualExpected);
	// Variable 
	int nGridPoints;
	int nRealizations;
	vector<int> nMemberParetoCurve;
	int nObjects;
	int sum_nMemberPareto;
	bool calcPCoverAndNotPSigma;
	double * ExpectedParetoPoints;
	double * pCoverArrayAPriori;
	size_t nExpectedParetoPoints;
	vector<vector<double>> totalParetoSet;
};

