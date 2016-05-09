#pragma once
#include <iostream>
#include <vector>
#include <numeric>
#include "mex.h"
using namespace std;

class determineParetoPoints
{
public:
	determineParetoPoints(size_t nSamples, size_t nObjectives, double *objectiveMatrix);
	~determineParetoPoints();

	void doCalculation(vector<bool> &paretoOptimalBool);
	bool checkIfPointIsParetooptimal(size_t currentIndex);
	const size_t nSamples;
	const size_t nObjectives;
	double *ObjectiveMatrix;
};

