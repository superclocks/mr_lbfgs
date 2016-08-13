#include "one_dim_opt.h" 


vector<double> OneDimOpt::QuadraticInterpolationSearch(double (*p)(double x, vector<double>& xk, vector<double>& pk), \
		double& l, double& u, double& x, vector<double>& xk, vector<double>& pk,double tol)
{
	int mm = 0;
	//1)
	double x1 = l;
	double x3 = u;
	double x0_aver = 1e99;

	//2)
	double x2 = 0.5 * (x1 + x3);
	double f1 = p(x1, xk, pk);
	double f2 = p(x2, xk, pk);
	double f3 = p(x3, xk, pk);

	//3)
	vector<double> r;
	while(true)
	{
		mm++;
		double x_aver = ((x2*x2 - x3*x3)*f1 + (x3*x3 - x1*x1)*f2 + (x1*x1 - x2*x2)*f3 )/(2*((x2 - x3)*f1 + (x3 - x1)*f2 + (x1 - x2)*f3));
		double f_aver = p(x_aver, xk, pk);
		if(fabs(x_aver - x0_aver) < tol)
		{
			x = x_aver;
			r.push_back(f_aver);
			r.push_back(x);
			return r;
		}
		//4)
		if(x1 < x_aver && x_aver < x2)
		{
			if(f_aver <= f2)
			{
				x3 = x2;
				f3 = f2;
				x2 = x_aver;
				f2 = f_aver;
			}
			else
			{
				x1 = x_aver;
				f1 = f_aver;
			}
		}
		else if(x2 < x_aver && x_aver < x3)
		{
			if(f_aver <= f2)
			{
				x1 = x2;
				f1 = f2;
				x2 = x_aver;
				f2 = f_aver;
			}
			else
			{
				x3 = x_aver;
				f3 = f_aver;
			}
		}
		x0_aver = x_aver;
	}

}


vector<double> OneDimOpt::GoldenSectionSearch(double (*p)(double x, vector<double>& xk, vector<double>& pk), \
		double& l, double& u,double& x, vector<double>& xk, vector<double>& pk,double tol)
{
	//1)
	double xlk = l;
	double xuk = u;
	//2)
	double I1 = xuk - xlk;
	double K = 1.618034;
	double Ikp1 = I1 / K;
	double xak = xuk - Ikp1;
	double xbk = xlk + Ikp1;
	double fak = p(xak, xk, pk);
	double fbk = p(xbk, xk, pk);
	//3)
	vector<double> r;
	while(true)
	{
		double Ikp2 = Ikp1 / K;
		if(fak >= fbk)
		{
			xlk = xak;
			//xuk = xuk;

			xak = xbk;
			xbk = xlk + Ikp2;

			fak = fbk;
			fbk = p(xbk, xk, pk);
		}
		else
		{
			//xlk = xlk;
			xuk = xbk;

			xak = xuk - Ikp2;
			xbk = xak;

			fbk = fak;
			fak = p(xak, xk, pk);
		}
		Ikp1 = Ikp2;
		//4)
		if(Ikp2 < tol || xak > xbk)
		{
			if(fak > fbk)
				x = 0.5 * (xbk + xuk);
			if(fak == fbk)
				x = 0.5 * (xak + xbk);
			if(fak < fbk)
				x = 0.5 * (xlk + xak);

			l = xlk;
			u = xuk;
			r.push_back(p(x, xk, pk));
			r.push_back(x);
			return r;
		} 
	}
}

