#include "determineParetoPoints.h"


determineParetoPoints::determineParetoPoints(size_t nSamplesInput, size_t nObjectivesInput, double * objectiveMatrix) :nSamples(nSamplesInput), nObjectives(nObjectivesInput)
{
	ObjectiveMatrix = objectiveMatrix;
}


determineParetoPoints::~determineParetoPoints()
{
}


void determineParetoPoints::doCalculation(vector<bool> &paretoOptimalBool){
	bool isNotDominated = false;
	for (size_t iSample = 0; iSample < nSamples; iSample++){
		paretoOptimalBool[iSample] = checkIfPointIsParetooptimal(iSample);
	}
}

bool determineParetoPoints::checkIfPointIsParetooptimal(size_t currentIndex){
	// Initialization
	bool testIsWorseInAllObjectives = false;

	for (size_t iSample = 0; iSample < nSamples; iSample++){
		if (iSample != currentIndex){
			
			// Assume yes in the beginning
			testIsWorseInAllObjectives = true;
			for (size_t iObjective = 0; iObjective < nObjectives;iObjective++){
				// Is current Point worse in the considered objective than all other point?
				testIsWorseInAllObjectives = testIsWorseInAllObjectives&(*(ObjectiveMatrix + currentIndex + iObjective*nSamples))>(*(ObjectiveMatrix + iSample + iObjective*nSamples));
			}
			// Return false if point is dominated by at least one other point in the data set
			if (testIsWorseInAllObjectives){
				return false;
			}
		}
	}

	// No point was found in the data set which dominates the current considered points
	return true;
}