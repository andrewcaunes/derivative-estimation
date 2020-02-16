/*
This file defines the necessary global constants.
Their respective use are documented below.
Feel free to modify them as you need.
*/

#ifndef ANNEX_H_INCLUDED
#define ANNEX_H_INCLUDED
#include <vector>
#include "kde.h"
#include <eigen-3.3.7\Eigen\Dense>
#include "LPRB.h"
#include "LPR_cubic.h"

using Eigen::MatrixXd;
const int n(300);  //Number of data point generated.
const int RSS_PRECISION(30); //Number of points evaluated for RSS computation.
const int SIZE(3); // Used for next line
const double BANDWIDTH_ORDERS[SIZE] = {0.01,0.1,1}; // This array stores the orders of bandwidth to be tested. Example below :
// For BANDWIDTH_ORDER = {0.01,0.1,1}, the RSS will be compared for the following bandwidths : 0.01, 0.02, ... , 0.09, 0.1, 0.2 , ... , 0.9, 1, 2, ..., 9
// The possible orders are 0.001, 0.01, 0.1, 1.
// In theory, the more data one have, the smaller could the bandwidth be, however, I have found that the best results are usually obtained for bandwidth around 0.1

const int NB_POINTS(300); // Specifies the number of data points in the written csv file for smoothed data. See writeData(..) function in annex.cpp.
const double RANGE(20); // Specifies the range for data points in the written csv file for smoothed data. See writeData(..) function in annex.cpp.
const double START(0); // Specifies the lowest data point in the written csv file for smoothed data. See writeData(..) function in annex.cpp.

const string file_name("C:/Users/andrew-pc/derivatives.csv"); //Specifies the name and location of the csv file for raw derivatives data.
const string file_name2("C:/Users/andrew-pc/smoothed_derivatives.csv"); //Specifies the name and location of the csv file for smoothed derivatives data.




typedef double dimensions[n];
void makeXp(MatrixXd& Xp, double X[], int n);
void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove);
bool sort_map(const std::vector<double>& u, const std::vector<double>& y);
void createData(double X[], double Y[], const int n);
void makeKernelDensityEstimator(KDE* kde, double X[], const int n);
void makeU(double U[], double X[], double Y[], KDE* kde, const int n);
void computeB(MatrixXd *Up, double Y[], const int n, double* B);
void computeSigma2(double* sigma2, double Y[], const int n);
void computeK(int *k, double B, double sigma2, const int n);
void computeWeights(dimensions* weights, const int n, double k, double U[]);
void computeDerivative(double Yd[], double Y[], double U[], dimensions* weights, const int n, double k);
void checkBandwidthRange(double* min, double* argmin, double* tmp, double order, LPRB* lprb);
void writeData(double Yn[], double U[], double X[], KDE* kde, LPR_cubic* final);

#endif // ANNEX_H_INCLUDED
