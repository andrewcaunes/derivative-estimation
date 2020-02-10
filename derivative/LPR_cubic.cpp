#include "LPR_cubic.h"
#include <C:\Users\andrew-pc\programmation\C\eigen-3.3.7\eigen-3.3.7\Eigen\Dense>
#include <math.h>
#include <iostream>
#include <dlib/optimization.h>
#include "LPRB.h"
#include "praxis.hpp"
#include "annex.h"

using namespace Eigen;

LPR_cubic::LPR_cubic(double Yn[], MatrixXd Xn, int nn, float hn)
:Xp(Xn),h(hn),n(nn)
{
    Y = Map<VectorXd>(Yn,nn);
    /*double starting_point[1];
    starting_point[0] = 1;
    typedef double (*funct)(double *x, int n);

        
    praxis(0.01, 1, 1, 1, starting_point, funct);*/

    LPRB obj = LPRB(Yn, Xn, nn, hn);

    double *d = new double;
    double* min = new double;
    *min = 10;
    float *argmin = new float;
    *argmin = 1;
    for (int i = 1; i < 5; i++) {
        *d = obj.RSS(0.2 * i);
        //cout << "\n" << *d << std::endl;
        cout << " for h = " << 0.2 * i << " , rss is : " << *d << endl;
        if (*d < *min) {
            *min = *d;
            *argmin = 0.2 * i;
            
        }
    }
    h = *argmin;
    cout << "\nh is " << h << " and the min is " << *min << endl;
    LPRB o = LPRB(Yn, Xn, nn, h);
    cout << " 10 : " << o.value(10) <<" and h : " << h;
    delete argmin;
    delete min;
    delete d;
    
    /*dlib::find_min(dlib::bfgs_search_strategy(), dlib::objective_delta_stop_strategy(0.01), obj.RSS, starting_point, -1);*/

}

double LPR_cubic::value(double x)
{
    double* weights = new double[n];
    double* dif_squared = new double;
    for(int i = 0 ; i < n ; i++){
        *dif_squared = pow((x - Xp(i,1)),2);
//        std::cout <<"\n i is " << i <<" and diff is " << *dif_squared <<" . " << Xp(i,1) <<"\n"; TEST
        *(weights+i) = 2.0/sqrt(3.14159)*(*dif_squared)/(pow(h,2))*exp(-1*(*dif_squared)/pow(h,2)); //Use of binomial kernel
    }
    VectorXd V = Map<VectorXd>(weights,n);
    MatrixXd W = V.asDiagonal();
    B = ((Xp.transpose()*W*Xp).inverse())*(Xp.transpose()*W*Y); // formula for weighted least squares resolution
//    B = (Xp.transpose()*W*Y); TEST

//    B = Xp.transpose()*Y; TEST
//    std::cout << Xp <<std::endl; TEST
//    std::cout << "\n voila W: \n" << W << "\n . \n";
//    std::cout << "\n voila B : \n" << B << "\n . \n";
    delete dif_squared;
    delete [] weights;
    return B(0)+B(1)*x+B(2)*x*x+B(3)*x*x*x;
}

//double LPR_cubic::RSS()
//{
//
//}

LPR_cubic::~LPR_cubic()
{
    //dtor
} 
