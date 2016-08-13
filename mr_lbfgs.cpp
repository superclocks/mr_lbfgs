/*************************************************************************
	> File Name: mr_lbfgs.h
	> Author: 
	> Mail: 
	> Created Time: 2015年09月28日 星期一 22时08分48秒
 ************************************************************************/

#include "mr_lbfgs.h"

MrLbfgs::MrLbfgs(double (*obj)(vector<double>& x), enum LossType loss_type, int m, int n,string file_path, double lamda, int itera, double tol) :
MathUnit(n+1 ), _obj(obj),_loss_type(loss_type),_file_path(file_path), _lamda(lamda),_m(1),_vectors_num(m), _N(n + 1), _itera(itera), _tol(tol)
{
	_compute_type = string("batch");
	vector<double> t(_N, 0);
	_w = t;
}
MrLbfgs::~MrLbfgs()
{}

void MrLbfgs::VectorDot()
{
	_dot_matrix.clear();
    vector<vector<double> >::iterator it1;
    vector<vector<double> >::iterator it2;
    for(it1 = _b.begin(); it1 != _b.end(); it1++)
    {
        vector<double> tmp;
        for(it2 = it1; it2 != _b.end(); it2++)
        {
            double dot = Dot(*it1, *it2);
            tmp.push_back(dot);
        }
        _dot_matrix.push_back(tmp);
    }
}

vector<double> MrLbfgs::Grad(vector<double>& x, double delta)
{
    vector<double> g;
    for(int i = 0; i < _N; i++)
    {
        double f0 = _obj(x);
        x[i] += delta; 
        double f1 = _obj(x);
        x[i] -= delta;
        g.push_back((f1 - f0) / delta);
    }
    return g;
}
double MrLbfgs::Obj(vector<double>& w)
{
	ifstream reader(_file_path.c_str(), ios::in);
	string line;
	vector<int> _id;
	vector<double> _xi;
	double _yi;
	vector<double> aver_g(_N ,0);
	double loss_val = 0.0;
	int N = 0;
	while(getline(reader, line))
	{
				if(line.compare("") == 0)
					break;
				Feature::Default(line, _id, _xi, _yi);
	            double wTx = Dot(_id, _xi, w);

				if(_loss_type == LogLoss_)
					loss_val += LogLoss::loss(wTx, _yi);
				else if(_loss_type == LogLoss01_)
					loss_val += LogLoss01::loss(wTx, _yi);
				_xi.clear();
				_id.clear();
				N++;
	}
	reader.close();
	return loss_val / N;
}
void MrLbfgs::L1(double lamda, vector<double>& w, vector<double>& g)
{
	vector<double>::iterator g_it = g.begin();
	for(vector<double>::iterator it = w.begin(); it != w.end(); it++)
	{
		double t = (*it) > 0 ? 1.0 : 0.0;
		*g_it += lamda * t;
		g_it++;
	}
}
 vector<double> MrLbfgs::BatchGrad(vector<double>& w) //批量计算目标函数
{
	   ifstream reader(_file_path.c_str(), ios::in);
		string line;
		vector<int> _id;
		vector<double> _xi;
		double _yi;
		vector<double> aver_g(_N ,0);
		int N = 1;
		while(getline(reader, line))
		{
					if(line.compare("") == 0)
						break;
					Feature::Default(line, _id, _xi, _yi);
		            double wTx = Dot(_id, _xi, w);
					double loss_val = 0.0;
					if(_loss_type == LogLoss_)
						loss_val = LogLoss::dloss(wTx, _yi);
					else if(_loss_type == LogLoss01_)
						loss_val = LogLoss01::dloss(wTx, _yi);
					vector<double> g = Multiply(loss_val, _id,  _xi);

					vector<double> t1 = Subtract(g, aver_g);
					vector<double> t2 = Multiply(1.0 / N, t1);
					aver_g = Add(aver_g, t2);
					N++;
					_xi.clear();
					_id.clear();
		}
		reader.close();
		L1(_lamda ,w, aver_g);
		return aver_g;
}

vector<double> MrLbfgs::InitGrad(vector<double>& x0)
{
    vector<double> sk(_N, 1e-7);
    vector<double> xk = x0;
    vector<double> xk_1;
    for(int i = 0; i < _m; i++)
    {
    	vector<double> t1 =  Multiply(i + 1, sk); //临时变量
        xk_1 = Add(xk, t1);
        vector<double> g2;
        vector<double> g1;
        if(_compute_type.compare("batch") == 0)
        {
        	 g2 = BatchGrad(xk_1);
        	 g1 = BatchGrad(xk);
        }
        else
        {
        	g2 = Grad(xk_1);
        	g1 = Grad(xk);
        }
        vector<double> g21 = Subtract(g2, g1);
        if(Sum(g21) == 0)
        {
            for(int j = 0; j < 10; j++)
            {
            	vector<double> t2 = Multiply(j + 1, sk); //临时变量
                vector<double> xk_1 = Add(xk, t2);
                vector<double> gt2;
                vector<double> gt1;
                if(_compute_type.compare("batch") == 0)
                {
                	gt1 = BatchGrad(xk_1);//临时变量
                	gt2 = BatchGrad(xk);//临时变量
                }
                else
                {
                	gt1 = Grad(xk_1);//临时变量
                	gt2 = Grad(xk);//临时变量
                }
                vector<double> t = Subtract(gt1, gt2);
                if(Sum(t) != 0)
                {
                    g21 = t;
                    break;
                }
            }
        }
        vector<double> xk_del = Subtract(xk_1, xk);
        _s_array.push_back(xk_del);
        _y_array.push_back(g21);
        xk = xk_1;
    }
    if(_compute_type.compare("batch") == 0)
    	_g = BatchGrad(xk_1);
    else
    	_g = Grad(xk_1);
    _b.assign(_s_array.begin(), _s_array.end());
    _b.insert(_b.end(), _y_array.begin(), _y_array.end());
    _b.push_back(_g);

    return xk_1;
}

vector<double> MrLbfgs::CallPK()
{
    vector<double> pk(_N, 0.0);
    for(int i = 0; i < _m * 2 + 1; i++)
    {
    	vector<double> t = Multiply(_delta[i], _b[i]); //临时变量
        pk = Add(pk, t);
    }
    return pk;
}

void MrLbfgs::VfLbfgsMR()
{
    //line 1 - 3
	_delta.clear();
    _delta.assign(_m * 2, 0);
    _delta.push_back(-1);
    //line 4 - 8
    vector<double> alfa;
    for(int i = _m - 1; i >= 0; i--)
    {
        int j = i;
        double s1 = 0.0;
        for(int k = 0; k < _m * 2 + 1; k++)
        {
            if(k > j)
                s1 += _delta[k] * _dot_matrix[j][k - j];
            else
                s1 += _delta[k] * _dot_matrix[k][j - k];
        }
        double s2 = 0.0;
        if(j > _m - 1 + j)
            s2 = _dot_matrix[_m - 1 + j][j - (_m - 1 + j)];
        else
            s2 = _dot_matrix[j][_m - 1 + j - j];
        double alfa_i = s1 / s2;
        alfa.push_back(alfa_i);
        _delta[_m - 1 + j] = _delta[_m - 1 + j] - alfa_i;
    }
    reverse(alfa.begin(), alfa.end());
    //line 9 - 11
    double numerator = 0.0;
    if ((_m - 1) > (_m * 2 - 1))
    	numerator = _dot_matrix[_m * 2 - 1][_m - 1 - (_m * 2 - 1)];
    else
    	numerator = _dot_matrix[_m - 1][_m * 2 - 1 - (_m - 1)];
    double denominator = _dot_matrix[_m * 2 - 1][0];
    double rat = numerator / denominator;
    for(int i = 0; i < _m * 2 + 1; i++)
        _delta[i] = rat * _delta[i];
    //line 12 - 16
    for(int i = 0; i < _m; i++)
    {
        int j = i;
        double s1 = 0.0;
        for(int k = 0; k < _m * 2 + 1; k++)
        {
            s1 += _m - 1 + j > k ? _delta[k] * _dot_matrix[k][_m - 1 + j - k] : \
        _delta[k] * _dot_matrix[_m - 1 + j][k - (_m - 1 + j)];
        }
        double s2 = j > _m - 1 + j ? _dot_matrix[_m - 1 + j][j - (_m - 1 + j)] : \
        _dot_matrix[j][_m - 1 + j - j];
        double beta = s1 / s2;
        _delta[j] = _delta[j] + alfa[i] - beta;
    }

}
double MrLbfgs::SearchObj(double x, vector<double>& xk, vector<double>& pk)
{
	vector<double> t = Multiply(x, pk);
	vector<double> t1 = Add(xk, t);
	double f;
	if(_compute_type.compare("batch") == 0)
		f = Obj(t1);
	else
		f = _obj(t1);
  // cout << f <<"," << x<< endl;
    return f;
}

void MrLbfgs::Optimization(vector<double>& x0)
{
	vector<double> xk;
	if(x0.size() == 0)
		xk = InitGrad(_w);
	else
		xk = InitGrad(x0);

    for(int i = 0; i < _itera; i++)
    {
    	 if(i + 1 >= _vectors_num)
    	        	_m = _vectors_num;
    	 else
    	        	_m = i + 1;
    	cout << "itera = " << i << " : "<< Obj(xk) <<endl;
        _x_list.push_back(xk);
        if(_compute_type.compare("batch") == 0)
        	_y_list.push_back(Obj(xk));
        else
        	_y_list.push_back(_obj(xk));
        VectorDot();
        vector<double> gk;
        if(_compute_type.compare("batch") == 0)
        	gk = BatchGrad(xk);
        else
        	gk = Grad(xk);
        VfLbfgsMR();
        vector<double> pk = CallPK();
        double x_opt;
        double l = -10.0;
        double u = 10.0;
        //vector<double> f = GoldenSectionSearch(&MrLbfgs::SearchObj, l, u, x_opt, xk, pk);
        vector<double> f1 = QuadraticInterpolationSearch(&MrLbfgs::SearchObj, l, u, x_opt, xk, pk);
        if(fabs(x_opt) < 1e-10)
        	break;
        vector<double> t2 = Multiply(x_opt, pk); //临时变量
        xk = Add(xk, t2);

     	vector<double> gk_1;
        if(_compute_type.compare("batch") == 0)
        	gk_1 = BatchGrad(xk);
        else
        	gk_1 = Grad(xk);
        vector<double> sk = Multiply(x_opt, pk);
        vector<double> yk = Subtract(gk_1, gk);
        if(i + 1 >= _vectors_num)
        {
        	_b.erase(_b.begin());
        	_b.insert(_b.begin() + _m - 1, sk);
        	_b.erase(_b.begin() + _m);
        	_b[_m * 2 - 1] = yk;
        	_b.push_back(gk_1);
        }
        else
        {
        	_b.insert(_b.begin() + _m, sk);
        	_b.insert(_b.begin() + _m * 2 + 1, yk);
        	_b[_m * 2 + 2] = gk_1;
        }


        if(Norm(gk) < _tol)
            break;
    }
    	_x_list.push_back(xk);
    	  if(_compute_type.compare("batch") == 0)
    		  _y_list.push_back(Obj(xk));
    	  else
    	      _y_list.push_back(_obj(xk));
        _w = xk;
}

void MrLbfgs::Print()
{
	cout << "Optimization infos:\n";
	int i = 0;
	for(vector<vector<double> >::iterator it1 = _x_list.begin(); it1 != _x_list.end(); it1++)
	{
		cout << "iter " << i ;
		cout << " y = " << _y_list[i];
		cout << " x = ";
		for(vector<double>::iterator it2 = (*it1).begin(); it2 != (*it1).end(); it2++)
		{
			cout << "," << *it2;
		}
		cout << endl;
		i++;
	}
}
vector<double> MrLbfgs::GoldenSectionSearch(double (MrLbfgs::*p)(double x, vector<double>& xk, vector<double>& pk), \
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
	double fak = (this->*p)(xak, xk, pk);
	double fbk =  (this->*p)(xbk, xk, pk);
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
			fbk = (this->*p)(xbk, xk, pk);
		}
		else
		{
			//xlk = xlk;
			xuk = xbk;
			xak = xuk - Ikp2;
			xbk = xak;
			fbk = fak;
			fak =  (this->*p)(xak, xk, pk);
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
			r.push_back( (this->*p)(x, xk, pk));
			r.push_back(x);
			return r;
		}
	}
}

vector<double> MrLbfgs::QuadraticInterpolationSearch(double (MrLbfgs::*p)(double x, vector<double>& xk, vector<double>& pk), \
		double& l, double& u, double& x, vector<double>& xk, vector<double>& pk,double tol)
{
	int mm = 0;
	//1)
	double x1 = l;
	double x3 = u;
	double x0_aver = 1e99;

	//2)
	double x2 = 0.5 * (x1 + x3);
	double f1 =  (this->*p)(x1, xk, pk);
	double f2 =  (this->*p)(x2, xk, pk);
	double f3 =  (this->*p)(x3, xk, pk);

	//3)
	vector<double> r;
	while(true)
	{
		mm++;
		double x_aver = ((x2*x2 - x3*x3)*f1 + (x3*x3 - x1*x1)*f2 + (x1*x1 - x2*x2)*f3 )/(2*((x2 - x3)*f1 + (x3 - x1)*f2 + (x1 - x2)*f3));
		double f_aver =  (this->*p)(x_aver, xk, pk);
		if(fabs(x_aver - x0_aver) < tol || mm > 50)
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
double MrLbfgs::Predict(vector<int>& id, vector<double>& x)
{
    if(_loss_type == LogLoss_)
    {
        double z = Dot(id, x, _w);
        double pi = LogLoss::deci(z);
        return pi;
        /*if(pi > 0.5)
            return 1;
        else
            return 0;*/
    }
    else if(_loss_type == LogLoss01_)
    {
        double z = Dot(id, x, _w);
        double pi = LogLoss01::deci(z);
        return pi;
        /*if(pi > 0.5)
            return 1;
        else
            return -1;*/
    }
}
void MrLbfgs::SaveModel(const char* path)
{
	ofstream writer(path, ios::out);
	int id = 0;
	for(vector<double>::iterator it = _w.begin(); it != _w.end(); it++)
	{
		writer << id++ << "\t" << *it << endl;
	}
	writer.close();
}
void MrLbfgs::Evaluate(const char* f)
{
    ifstream reader(f, ios::in);
    if(!reader.is_open())
    {
        cout<< "Open file failed." << endl;
        assert(-1);
    }
    string s;
    vector<int> ids;
    vector<double> x;
    int k = 0;
    int count = 0;
    while(getline(reader, s))
    {
    	if(s.compare("") == 0)
    		break;
        /*vector<string> ele;
        split(ele, s, is_any_of(" "));
        vector<string>::iterator it = ele.begin();
        int y = lexical_cast<int>(ele[0]);
        it++;
    	x.clear();
    	ids.clear();
        for(; it != ele.end(); it++)
        {
            vector<string> tmp;
            split(tmp, *it, is_any_of(":"));
            int id = lexical_cast<int>(tmp[0]);
            double val =  lexical_cast<double>(tmp[1]);
            ids.push_back(id);
            x.push_back(val);
        }*/
        double y;
        ids.clear();
        x.clear();
        Feature::Default(s, ids, x, y);
        int y_pred;
        double pi = Predict(ids, x);
        if(_loss_type == LogLoss_)
            y_pred = pi > 0.5 ? 1 : -1;
        else if(_loss_type == LogLoss01_)
            y_pred = pi > 0.5 ? 1 : 0;
//cout<<  "pi = " << pi <<  "y_pred = " << y_pred << " y = " << y << endl;
        if(y_pred != y)
            k++;
        count++;
    }
    cout << "error rating = " << (double)k / count << endl;
}

