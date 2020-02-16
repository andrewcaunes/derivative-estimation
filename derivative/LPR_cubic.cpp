/*
The LPR_cubic class is used to smooth some discrete data by locally fitting a cubic polynomial to it (Local Polynomial Regression method).
It functions by first trying multiple bandwidths in some range with a binomial kernel. It then adapts this bandwidth to a Gaussian kernel.
The fitting is done in real time for any value when the value(double x) method.

The constructor LPR_cubic(double Yn[], MatrixXd Xn, int nn, float hn) executes this method.
One can call the smoothed function through the value(double x) method, after having created an object of the class.
*/





#include "LPR_cubic.h"
#include <eigen-3.3.7\Eigen\Dense>
#include <math.h>
#include <iostream>
#include "LPRB.h"
#include "annex.h"
#include <ctime>

using namespace Eigen;

//Constructor
LPR_cubic::LPR_cubic(double Yn[], MatrixXd Xn, int nn, float hn)
:Xp(Xn),h(hn),n(nn)
{
    Y = Map<VectorXd>(Yn, nn); // Initialization

    // We will now search for the best bandwidth possible with a binomial kernel

    LPRB *lprb = new LPRB(Yn, Xn, nn, hn);      // LPRB object is a local cubic regressor using a certain bandwidth
    double* min = new double;
    *min = 100000;
    double* argmin = new double;
    *argmin = hn;
    const float adaptation(1.01431);
    double* tmp = new double;
    *tmp = 0;

    for (int i = 0; i < SIZE; i++) {
        checkBandwidthRange(min, argmin, tmp, BANDWIDTH_ORDERS[i], lprb); //This function is defined in annex.cpp. Feel free to modify the constant BANDWIDTH_ORDERS in annex.h.
    }
    
    h = *argmin;  // The bandwidth of the LPR_cubic object is set to the one minimizing the RSS
    h = h / adaptation; // The bandwidth is then adapted to a Gaussian Kernel
    
    delete argmin;
    delete min;
    delete tmp;
    delete lprb;

}

double LPR_cubic::value(double x)
{
    
    double* weights = new double[n];
    double* dif_squared = new double;

    for(int i = 0 ; i < n ; i++){ //weights array is filled with the kernel values for each data point.
        *dif_squared = pow((x - Xp(i,1)),2);
        *(weights+i) = 1.0/sqrt(2*3.14159)*exp(-1*(*dif_squared)/2/pow(h,2)); //Use of Gaussian kernel
    }

    VectorXd V = Map<VectorXd>(weights,n);
    MatrixXd W = V.asDiagonal();

    B = ((Xp.transpose()*W*Xp).inverse())*(Xp.transpose()*W*Y); // weighted least squares resolution formula
    

    delete dif_squared; //dynamic memory is released
    delete [] weights;
    

    return B(0)+B(1)*x+B(2)*x*x+B(3)*x*x*x; //The 
}

//Destructor
LPR_cubic::~LPR_cubic()
{
    //dtor
} 
