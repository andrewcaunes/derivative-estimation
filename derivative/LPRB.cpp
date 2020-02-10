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
//        std::cout <<"\n i is " << i <<" and diff is " << *dif_squared <<" . " << Xp(i,1) <<"\n"; TEST
        *(weights+i) = 2.0/sqrt(3.14159)*(*dif_squared)/(pow(h,2))*exp(-1*(*dif_squared)/pow(h,2)); //Use of binomial kernel
    }
    int i(27), j(28), k(29);
    //std::cout << "\nweights : " << "w"<<x<<","<<i<<" : "<< *(weights+i) << "  w" << x << "," << j << " : " << *(weights + j) << "  w" << x << "," << k << " : " << *(weights + k) << std::endl;
    i = 0;
    j = 1;
    k = 2;
    /*if ((int)x % 30 == 0) {
        std::cout << "\n" << x << "\n";
    }*/
    //std::cout << "\nweights : " << "w" << x << ","<<i<<" : "<< *(weights+i) << "  " << x << "," << j << " : " << *(weights + j) << "  w" << x << "," << k << " : " << *(weights + k) << std::endl;
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

double LPRB::RSS(float hn)
{
    h = hn;
    double sum(0);
    for(int i=0 ; i<n ; i++){
        //std::cout << "\nRSS: i = " << i << ", value(Xp(i,1))= " << value(Xp(i, 1)) << ", Y[i]=" << Y[i] << std::endl;
        sum = sum + pow(LPRB::value(Xp(i,1))-Y[i],2);
    }
    //std::cout << "\n RSS is " << sum << std::endl;
    return sum;
}

void LPRB::setH(float hn){
    h = hn;
}

LPRB::~LPRB()
{
    //dtor
}
