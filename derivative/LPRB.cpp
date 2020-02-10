#include "LPRB.h"
#include <C:\Users\andrew-pc\programmation\C\eigen-3.3.7\eigen-3.3.7\Eigen\Dense>
#include <math.h>
#include <iostream>

using namespace Eigen;

LPRB::LPRB(double Yn[], MatrixXd Xn, int nn, float hn)
:Xp(Xn),h(hn),n(nn)
{
    Y = Map<VectorXd>(Yn,nn);
}

double LPRB::value(double x)
{
    double *weights = new double[n];
    double* dif_squared = new double;
    for(int i = 0 ; i < n ; i++){
        *dif_squared = pow((x - Xp(i,1)),2);
        *(weights+i) = 2.0/sqrt(3.14159)*(*dif_squared)/(pow(h,2))*exp(-1*(*dif_squared)/pow(h,2)); //Use of binomial kernel
    }
    
    VectorXd V = Map<VectorXd>(weights,n);
    MatrixXd W = V.asDiagonal();
    B = ((Xp.transpose()*W*Xp).inverse())*(Xp.transpose()*W*Y); // formula for weighted least squares resolution


    delete dif_squared;
    delete [] weights;

    return B(0)+B(1)*x+B(2)*x*x+B(3)*x*x*x;
}

double LPRB::RSS(float hn)
{
    h = hn;
    double sum(0);
    for(int i=0 ; i<n ; i++){
        sum = sum + pow(LPRB::value(Xp(i,1))-Y[i],2);
    }
    return sum;
}

void LPRB::setH(float hn){
    h = hn;
}

LPRB::~LPRB()
{
    //dtor
}
