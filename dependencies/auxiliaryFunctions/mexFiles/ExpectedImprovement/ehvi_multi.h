//Include this if you want to calculate the EHVI of multiple individuals at the same
//time. This is more efficient than repeatedly calling an EHVI function on the same
//population.
#include "helper.h"
#include "mex.h"
#include "ehvi_consts.h"
#include <vector>
#include <math.h>
//#include "ehvi_hvol.h"
//#include <deque>
#include <algorithm>
//#include <iostream>
//#include <array>
using namespace std;

struct mus{ //Holds mean & variance for a Gaussian distribution
  //double mu[DIMENSIONS];
  // //double s[DIMENSIONS];
  //std::array<double, DIMENSIONS> mu;
  //std::array<double, DIMENSIONS> s;
	std::vector<double> mu;
	std::vector<double> s;
 
	mus()
		:mu(DIMENSIONS, 0), s(DIMENSIONS, 0){};
};

//vector<double> ehvi3d_5term(deque<individual*> P, double r[], vector<mus *> & pdf);
vector<double> ehvi3d_sliceupdate(vector<individual> P, vector<double> r, vector<mus > & pdf);
