/*************************************************************************
	> File Name: mr_lbfgs.h
	> Author: 
	> Mail: 
	> Created Time: 2015年09月28日 星期一 22时08分48秒
 ************************************************************************/

#ifndef _MR_LBFGS_
#define _MR_LBFGS_
//

#include "math_unit.h"
#include "one_dim_opt.h"
#include "feature.h"
#include "loss_func.h"
#include <stdio.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <iostream>
#include <assert.h>
#include <fstream>
#include <boost/lexical_cast.hpp>
using namespace std;
using namespace boost;
//#define
enum LossType
{
    LogLoss_,
    LogLoss01_,
    HingeLoss_,
    SquaredHingeLoss_,
    SmoothHingeLoss_
};

enum BasisFunc
{
    Default_,
    Gaussian_

};
class MrLbfgs : public MathUnit
{
    private:
        int _itera;
        double _tol;
        vector<vector<double> > _s_array;
        vector<vector<double> >_y_array;
        vector<vector<double > > _dot_matrix;
        vector<vector<double> > _b;
        vector<double> _g;
        vector<double> _delta;
        vector<double> _w;
        double _lamda;
        vector<vector<double> > _x_list;
        vector<double> _y_list;
        int _m; //计算近似矩阵所需要的向量个数
        int _vectors_num;
        int _N; //特征个数
        string _compute_type;
        string _file_path; //训练样本的数据文件
        enum LossType _loss_type; //损失函数类型
       double (*_obj)(vector<double>& x); //损失函数指针
    public:
        MrLbfgs( double (*obj)(vector<double>& x), enum LossType loss_type ,int m, int n, string file_path = "" , double lamda = 0.1, int itera = 50, double tol = 1e-6);
        ~MrLbfgs();

        void VectorDot();
        vector<double> InitGrad(vector<double>& x0);
        //void VfLbfgs();
        void VfLbfgsMR();
        vector<double> CallPK();
        vector<double> Grad(vector<double>& x, double delta = 1e-9); //数值梯度，计算不可解析的目标函数
        vector<double> BatchGrad(vector<double>& x); //批量计算解析梯度，计算逻辑回归，SVM的绞链函数，线性回归等函数
        double Obj(vector<double>& x); //批量计算目标函数
        void Optimization(vector<double>& x0);
        double SearchObj(double x, vector<double>& xk, vector<double>& pk);
        void Print();
        vector<double> GoldenSectionSearch(double (MrLbfgs::*p)(double x, vector<double>& xk1, vector<double>& pk1),
            			double& l, double& u, double& x, vector<double>& xk, vector<double>& pk, double tol = 1e-6);
        vector<double> QuadraticInterpolationSearch(double (MrLbfgs::*p)(double x, vector<double>& xk1, vector<double>& pk1),
            			double& l, double& u, double& x, vector<double>& xk, vector<double>& pk, double tol = 1e-6);
        void L1(double lamda, vector<double>& w, vector<double>& g);
        void SaveModel(const char* path);
        double Predict(vector<int>& id, vector<double>& x);
        void Evaluate(const char* f);
};

#endif
