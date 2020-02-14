/*
The LPR_cubic class is used to smooth some discrete data by fitting a cubic polynomial to it (Local Polynomial Regression method).
It functions by first trying multiple bandwidths in some range with a binomial kernel. It then adapts this bandwidth to a Gaussian kernel.

The constructor LPR_cubic(double Yn[], MatrixXd Xn, int nn, float hn) executes this method.
One can call the smoothed function through the value(double x) non static member function, after having created an object of the class.
*/





#include "LPR_cubic.h"
#include <C:\Users\andrew-pc\programmation\C\eigen-3.3.7\eigen-3.3.7\Eigen\Dense>
#include <math.h>
#include <iostream>
#include "LPRB.h"
#include "annex.h"

using namespace Eigen;

LPR_cubic::LPR_cubic(double Yn[], MatrixXd Xn, int nn, float hn)
:Xp(Xn),h(hn),n(nn)
{
    Y = Map<VectorXd>(Yn, nn); //simple initialization of Y

    // We will now search for the best bandwidth possible with a binomial kernel

    LPRB obj = LPRB(Yn, Xn, nn, hn);      // LPRB object is a local cubic regressor using a certain bandwidth
    double *d = new double;
    double* min = new double;
    *min = 100000;
    float *argmin = new float;
    *argmin = hn;
    const float adaptation(1.01431);

    for (int i = 1; i < 10; i++) {      // Here we try a range of bandwiths and compute the residual sum of square of the regression for each
        *d = obj.RSS(0.01 * i);
        if (*d < *min && *d>10) {
            *min = *d;
            *argmin = 0.01 * i;
        }
    }
    for (int i = 1; i < 10; i++) {      // Here we try a range of bandwiths and compute the residual sum of square of the regression for each
        *d = obj.RSS(0.1 * i);
        if (*d < *min && *d>10) {
            *min = *d;
            *argmin = 0.1 * i;
        }
    }
    for (int i = 1; i < 10; i++) {      // Here we try a range of bandwiths and compute the residual sum of square of the regression for each
        *d = obj.RSS(0.5 * i);
        if (*d < *min && * d>10) {
            *min = *d;
            *argmin = 0.5 * i;
        }
    }
    for (int i = 1; i < 10; i++) {      // Here we try a range of bandwiths and compute the residual sum of square of the regression for each
        *d = obj.RSS(0.001 * i);
        if (*d < *min && *d>10) {
            *min = *d;
            *argmin = 0.001 * i;
        }
    }
    h = *argmin;  // The bandwidth of the LPR_cubic object is set to the one minimizing the RSS
    h = h / adaptation; // The bandwidth is then adapted to a Gaussian Kernel

    std::cout << "\nh is " << h << " and the min is " << *min << std::endl;
    
    delete argmin;
    delete min;
    delete d;
}

double LPR_cubic::value(double x)
{
    waits = vector<double>();
    double* weights = new double[n];
    double* dif_squared = new double;

    for(int i = 0 ; i < n ; i++){
        *dif_squared = pow((x - Xp(i,1)),2);
        *(weights+i) = 1.0/sqrt(2*3.14159)*exp(-1*(*dif_squared)/2/pow(h,2)); //Use of Gaussian kernel
        waits.push_back(*(weights+i));
    }

    VectorXd V = Map<VectorXd>(weights,n);
    MatrixXd W = V.asDiagonal();

    B = ((Xp.transpose()*W*Xp).inverse())*(Xp.transpose()*W*Y); // weighted least squares resolution

    delete dif_squared; //dynamic memory is released
    delete [] weights;

    return B(0)+B(1)*x+B(2)*x*x+B(3)*x*x*x;
}


LPR_cubic::~LPR_cubic()
{
    //dtor
} 
