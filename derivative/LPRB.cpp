#include "LPRB.h"
#include <C:\Users\andrew-pc\programmation\C\eigen-3.3.7\eigen-3.3.7\Eigen\Dense>
#include <math.h>
#include <iostream>
#include <ctime>
#include "annex.h"

using namespace Eigen;

//Constructor
LPRB::LPRB(double Yn[], MatrixXd Xn, int nn, float hn)
:Xp(Xn),h(hn),n(nn)
{
    Y = Map<VectorXd>(Yn,nn);
}

double LPRB::value(double x) //Method to give prediction at data point x
{

    double *weights = new double[n];
    double* dif_squared = new double;
    for(int i = 0 ; i < n ; i++){ //weights array is filled with the kernel values for each data point.
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

double LPRB::RSS(float hn) //Method to compute the residual sum of squares of the regression with bandwidth hn
{
        
    h = hn;
    double sum(0);
    for(int i=0 ; i<n ; i++){
        if (i % int((double)n / (double)RSS_PRECISION) == 0) {
            sum += pow(this->value(Xp(i, 1)) - Y[i], 2);
        }        
    }
    sum += pow((this->value(0) - Y[0]),2);  // These two lines are used to prevent extreme misfitting on the sides (that can happen when the bandwidth is small)
    sum += pow((this->value(1) - Y[n-1]),2);
    return sum;
}

//Destructor
LPRB::~LPRB()
{
    //dtor
}
