// derivative.cpp : Defines the entry point for the application.

/* This program is a C++ implementation of a method for derivative estimation in random design.
The methode is thoroughly described in "Derivative Estimation in Random Design", from Yu Liu, Kris De Brabanter.
Details about the method are given in annex.cpp, along with an explanation of each function.
*/

/*
This file contains the main function of the program.
It uses multiple functions from annex.cpp.
It can be used to (optionally) generate random data and then compute the derivative of the transformation and then (optinally) store the data in a csv file.
*/

#include "derivative.h"


using Eigen::MatrixXd;
using namespace std;


int main()
{
    // INITIALIZATION

    
    double *X = new double[n];
    double *Y = new double[n];
    double *U = new double[n];
    double B(10);
    double sigma2(0);
    int k(10);
    MatrixXd* Up = new MatrixXd;
    double* Yd = new double[n];
    double* Yn = new double[n - 2];

    LPR_cubic* final;
    dimensions* weights = new dimensions[n];
    KDE* kde = new KDE();
    

    // START 
    
    // Random data is generated with the createData function. The type of distribution and the transformation between X and Y 
    // can be chosen in annex.cpp by changing these manually (the library used is <random>).
    cout << "\n Generating data..\n.\n." << endl;
    createData(X, Y, n);
    cout << "\n Data generated successfully." << endl;

    // Data is transformed so that it is initially uniformly distributed. It is also ordered.
    cout << "\n Probability integral transform of data..\n  .\n  ." << endl;
    makeKernelDensityEstimator(kde, X, n);
    makeU(U, X, Y, kde, n);
    makeXp(*Up, U, n);
    cout << "\n Probability integral transform successful." << endl;
    
    // The optimal number of quotients to use in the computation of derivatives is computed here.
    // The computation of B can take some time since the data is smoothed first, as indicated in the paper.
    cout << "\n Computing k (nbr of quotients to use for derivative computation)..\n  .\n  ." << endl;
    computeB(Up, Y, n, &B);
    computeSigma2(&sigma2, Y, n);    
    computeK(&k, B, sigma2, n);
    cout << "\n Optimal k computed successfully: k = " << k << endl;


    // The weights for each quotient is computed.
    cout << "\n Computing weights for quotients..\n  .\n  ." << endl;
    computeWeights(weights, n, k, U);
    cout << "\n Weights computed successfully." << endl;

    // The derivative data is computed.
    cout << "\n Computing derivative estimates..\n  .\n  ." << endl;
    computeDerivative(Yd, Y, U, weights, n, k);
    cout << "\n Derivative estimates computed successfully." << endl;

    
    delete[] Y;
    delete[] weights;
    // The following lines remove boundary data which is often very imprecise
    removeRow(*Up, n-1);
    removeRow(*Up, 0);
    for (int i = 0; i < n - 2; i++) {
        Yn[i] = Yd[i + 1];
        U[i] = U[i + 1];
        X[i] = X[i + 1];
    }
    delete[] Yd;

    // The derivative data is smoothed with a cubic polynomial local regression
    cout << "\n Smoothing data..\n  .\n  ." << endl;
    final = new LPR_cubic(Yn, *Up, n-2, 1);
    cout << "\n Data smoothed successfully" << endl;

    // The derivative data is written, before and after smoothing, in a csv file at the location specified in annex.h (by constants at the top of the file).
    cout << "\n Writing down data in csv file..\n  .\n  ." << endl;
    writeData(Yn, U, X, kde, final);
    cout << "\n Data written down in csv file successfully." << endl;

    // END
    
    delete final;
    delete kde;
    delete[] X;

    return 0;
}