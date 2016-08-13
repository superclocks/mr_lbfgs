#ifndef _ONEDIMOPT_
#define _ONEDIMOPT_
#include <algorithm>  
#include <vector>
#include <math.h>
#include "mr_lbfgs.h"
using namespace std;  
class OneDimOpt  
{  
	public:  
		//参考《Practical Optimization》中4.4节
    	vector<double> GoldenSectionSearch(double (*p)(double x, vector<double>& xk1, vector<double>& pk1),
    			double& l, double& u, double& x, vector<double>& xk, vector<double>& pk, double tol = 1e-6);
    	//参考《Practical Optimization》中4.5节
    	vector<double> QuadraticInterpolationSearch(double (*p)(double x, vector<double>& xk1, vector<double>& pk1),
    			double& l, double& u, double& x, vector<double>& xk, vector<double>& pk, double tol = 1e-6);
};  
#endif  
