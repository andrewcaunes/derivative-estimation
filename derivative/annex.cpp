/* 
This file contains most of the functions used in the main function.
annex.h defines the necessary constants.
*/

#include <iostream>
#include <eigen-3.3.7\Eigen\Dense>
#include <vector>
#include <random>
#include "kde.h"
#include <cmath>
#include "LPR_cubic.h"
#include "annex.h"
#include "LPRB.h"
#include <fstream>

 using namespace Eigen;
 using namespace std;

 /**
  * \brief function to make the matrix of the powers of the independent variables
  * \param Xp  Reference to matrix
  * \param X  Independent variable data
  * \param n  Number of data entries
  * \return None
  */
 void makeXp(MatrixXd& Xp, double X[], int n) {
     Xp.resize(n, 4);
     for (int i = 0; i < n; i++) {
         for (int j = 0; j < 4; j++) {
             Xp(i, j) = pow(X[i], j);
         }
     }
 }

 /**
  * \brief defines how to compare two vector<double> object
  * \param u  Reference to first object
  * \param u  Reference to second object
  * \return boolean from comparison of first element of both objects
  */
bool sort_map(const std::vector<double> &u,const std::vector<double> &y){
    return u[0] < y[0];
}

/**
  * \brief Removes a row from a MatrixXd object
  * \param matrix  Matrix
  * \param rowToRemove  number of row to remove
  * \return None
  */
void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove)
{
    unsigned int numRows = matrix.rows() - 1;
    unsigned int numCols = matrix.cols();

    if (rowToRemove < numRows)
        matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) = matrix.block(rowToRemove + 1, 0, numRows - rowToRemove, numCols);

    matrix.conservativeResize(numRows, numCols);
}

/**
  * \brief Generates random data
  * \param X  Array of doubles receiving the data generated from a specific random distribution
  * \param Y  Array of doubles receiving the transformed data after a specific transformation and optionally adding some random errors
  * \param n  Number of data entries
  * \return None
  */
void createData(double X[], double Y[], const int n) {
    // This function generates some random data by using the <random> library.
    // You can easily modify the type of generated data by changing the class and parameters of "distribution1" (independent variable X data distribution)
    // , the class of "distribution2" (errors), and the transformation between Y and X.
    // Initially, the independent variable data is uniformly distributed in [0,20] and the transformation is merely a cosine.

    default_random_engine generator;
    uniform_real_distribution<double> distribution1(0,20);
    normal_distribution<double> distribution2(0, 0.2 * 0.2);
    for (int i = 0; i < n; i++) {
        *(X + i) = distribution1(generator);
        *(Y + i) = cos(*(X + i)) + distribution2(generator);
    }
}

/**
  * \brief Estimates the density of a distribution
  * \param kde  KDE object from the kde library (kernel density estimation library) by Tim Nugent (c) 2014
  * \param X  Array of doubles containing the data generated from a specific random distribution
  * \param n  Number of data entries
  * \return None
  */
void makeKernelDensityEstimator(KDE* kde, double X[], const int n) {
    kde->set_kernel_type(1);
    kde->set_bandwidth_opt_type(1);

    for (int i = 0; i < n; i++) {
        kde->add_data(*(X + i));
    }
}

/**
  * \brief Uses probability integral transform to get equivalent data with uniform distribution in [0,1] and orders U,X and Y arrays
  * \param U  Array of doubles receiving the uniformly distributed data
  * \param X  Array of doubles containing the data generated from a specific random distribution
  * \param Y  Array of doubles containing the transformed data
  * \param kde KDE object from the kde library (kernel density estimation library) by Tim Nugent (c) 2014
  * \param n  Number of data entries
  * \return None
  */
void makeU(double U[], double X[], double Y[], KDE* kde, const int n) {

    for (int i = 0; i < n; i++) {
        U[i] = kde->cdf(*(X + i));
    }

    vector<vector<double>> map;
    for (int i = 0; i < n; i++) {
        map.push_back(vector<double>(0));
        map[i].push_back(U[i]);
        map[i].push_back(*(Y + i));
        map[i].push_back(X[i]);
    };

    sort(map.begin(), map.end(), sort_map);

    for (int i = 0; i < n; i++) {
        U[i] = map[i][0]/map[n-1][0];
        *(Y + i) = map[i][1];
        X[i] = map[i][2];
    }

}

/**
  * \brief Computes the sup|r2(u)| for r2 the second derivative of a transformation and u between 0 and 1.
  * \param Up  Pointer to matrix of powers of independent variable uniformly distributed in [0,1]
  * \param Y  Array of doubles containing the transformed data
  * \param n  Number of data entries
  * \param B  Pointer to double receiving the value of the sup
  * \return None
  */
void computeB(MatrixXd* Up, double Y[], const int n, double* B) {
    
    LPR_cubic* approx = new LPR_cubic(Y, *Up, n, 1);
    double* max = new double;
    *max = 0;
    double* tmp = new double;
    *tmp = 0;
    for (int i = 0; i < n - 1; i++) {
        if (i % 5 == 0) {
            approx->value((double(i) / double(n)));
            *tmp = abs((approx->value((double(i) / double(n))) - approx->value((double(i + 1) / double(n)))));
            *tmp /= (1 / double(n));
            if (*tmp > * max) {
                *max = *tmp;
            }
        }
    }

    *B = *max;

    delete approx;
    delete max;
    delete tmp;
}

/**
  * \brief computes the standard deviation of errors using Hall’s sqrt(n)-consistent estimator
  * \param sigma2  Pointer to double receiving the value
  * \param Y  Array of doubles containing the transformed data
  * \param n  Number of data entries
  * \return None
  */
void computeSigma2(double* sigma2, double Y[], const int n) {
    for (int i = 0; i < n - 2; i++) {
        *sigma2 += pow((0.809 * *(Y + i) - 0.5 * *(Y + i + 1) - 0.309 * *(Y + i + 2)), 2);
    }
    *sigma2 *= (1 / ((double)n - 2));
}

/**
  * \brief Computes the optimal number of quotients to be using in the derivative estimation
  * \param k  Pointer to int receiving the value
  * \param B  Parameter of transform (supposed to be computed with computeB(..))
  * \param sigma2  Parameter of transform (supposed to be computed with computeSigma2(..))
  * \param n  Number of data entries
  * \return None
  */
void computeK(int *k, double B, double sigma2, const int n) {
    double* tmp = new double;
    double* min = new double;
    *min = 10000;
    int* argmin = new int;
    *argmin = 10;
    tmp = new double;

    for (int i = 1; i < (int)(((double)n - 1) / 2); i++) {
        *tmp = pow(B * 3 * i * ((double)i + 1) / 4 / ((double)n + 1) / (2 * (double)i + 1), 2);
        *tmp += 3 * sigma2 * ((double)n + 1) * ((double)n + 1) / i / ((double)i + 1) / (2 * (double)i + 1);
        if (*tmp < *min) {
            *min = *tmp;
            *argmin = i;
        }
    }
    *k = *argmin;

    delete argmin;
    delete min;
    delete tmp;
}

/**
  * \brief  Computes the weights to be applied to the quotients in the derivative estimation.
  * \param  weights Bidimensional array of doubles receiving the weights
  * \param  n  Number of data entries
  * \param  k  optimal number of quotients to be using in the derivative estimation (should be found using computeK(..))
  * \param U  Array of doubles containing the uniformly distributed data in [0,1] (should be found using makeU(..))
  * \return  None
  */
void computeWeights(dimensions *weights, const int n, double k, double U[]) {
    double* sum = new double;
    int* ki = new int;
    *ki = 0;
    *sum = 0;

    for (int i = 1; i < n; i++) {

        for (int j = 0; j < k; j++) {
            if (i <= n - k && i >= k + 1) {

                *sum = 0;
                for (int l = 1; l <= k; l++) {
                    *sum += pow(U[i + l - 1] - U[i - l - 1], 2);
                }
                weights[i][j] = pow((U[i + j - 1] - U[i - j - 1]), 2) / *sum;
            }
            else if (i < k + 1) {
                *ki = i - 1;
                *sum = 0;
                if (j <= *ki) {
                    for (int l = 1; l <= *ki; l++) {
                        *sum += pow(U[i + l - 1] - U[i - l - 1], 2);
                    }
                    for (int l = *ki + 1; l <= k; l++) {
                        *sum += pow(U[i + l - 1] - U[i - 1], 2);
                    }
                    weights[i][j] = pow(U[i + j - 1] - U[i - j - 1], 2) / *sum;
                }
                else if (j > * ki) {
                    for (int l = 1; l <= *ki; l++) {
                        *sum += pow(U[i + l - 1] - U[i - l - 1], 2);
                    }
                    for (int l = *ki + 1; l <= k; l++) {
                        *sum += pow(U[i + l - 1] - U[i - 1], 2);
                    }
                    weights[i][j] = pow(U[i + j - 1] - U[i - 1], 2) / *sum;
                }
            }
            else if (i > n - k) {
                *ki = n - i;
                *sum = 0;
                if (j <= *ki) {
                    for (int l = 1; l <= *ki; l++) {
                        *sum += pow(U[i + l - 1] - U[i - l - 1], 2);
                    }
                    for (int l = *ki + 1; l <= k; l++) {
                        *sum += pow(U[i - 1] - U[i - l - 1], 2);
                    }
                    weights[i][j] = pow(U[i + j - 1] - U[i - j - 1], 2) / *sum;
                }
                else if (j > * ki) {
                    for (int l = 1; l <= *ki; l++) {
                        *sum += pow(U[i + l - 1] - U[i - l - 1], 2);
                    }
                    for (int l = *ki + 1; l <= k; l++) {
                        *sum += pow(U[i - 1] - U[i - l - 1], 2);
                    }
                    weights[i][j] = pow(U[i - 1] - U[i - j - 1], 2) / *sum;
                }
            }
        }
    }

    delete sum;
    delete ki;
}

/**
  * \brief  Computes the derivative estimates of transformed data.
  * \param Yd  Array of doubles receiving the derivative estimates
  * \param Y  Array of doubles containing the transformed data
  * \param U  Array of doubles containing the uniformly distributed data in [0,1] (should be found using makeU(..))
  * \param weights Bidimensional array of doubles containing the weights (should be found using computeWeights(..))
  * \param n  Number of data entries
  * \param k  optimal number of quotients to be using in the derivative estimation (should be found using computeK(..))
  * \return None
  */
void computeDerivative(double Yd[], double Y[], double U[], dimensions* weights, const int n, double k) {
    
    int* ki = new int;
    *ki = 0;

    for (int i = 1; i < n; i++) {
        Yd[i - 1] = 0;
        if (i <= n - k && i >= k + 1) {
            for (int j = 1; j <= k; j++) {
                Yd[i - 1] += weights[i - 1][j - 1] * (*(Y + i + j - 1) - *(Y + i - j - 1)) / (U[i + j - 1] - U[i - j - 1]);

            }
        }
        else if (i < k + 1) {
            *ki = i - 1;
            if (i > 1) {
                for (int j = 1; j <= *ki; j++) {
                    if (i == 1) {
                    }
                    Yd[i - 1] += weights[i - 1][j - 1] * (*(Y + i + j - 1) - *(Y + i - j - 1)) / (U[i + j - 1] - U[i - j - 1]);

                }
            }
            for (int j = *ki + 1; j <= k; j++) {
                Yd[i - 1] += weights[i - 1][j - 1] * (*(Y + i + j - 1) - *(Y + i - 1)) / (U[i + j - 1] - U[i - 1]);

            }

        }
        else if (i > n - k) {
            *ki = n - i;
            for (int j = 1; j <= *ki; j++) {
                Yd[i - 1] += weights[i - 1][j - 1] * (*(Y + i + j - 1) - *(Y + i - j - 1)) / (U[i + j - 1] - U[i - j - 1]);
            }
            for (int j = *ki + 1; j <= k; j++) {
                Yd[i - 1] += weights[i - 1][j - 1] * (*(Y + i - 1) - *(Y + i - j - 1)) / (U[i - 1] - U[i - j - 1]);
            }
        }
    }
    delete ki;

}

/**
  * \brief Tries multiple bandwidth ranges for Local Cubic Polynomial Regression, depending on constant array BANDWIDTH_ORDERS
  * \param min  pointer to min of RSS
  * \param argmin  pointer to argmin of RSS
  * \param tmp  temporary pointer
  * \param order  order of range
  * \param lprb  Local Polynomial Regression object (defined in LPRB.h and LPRB.cpp)
  * \return None
  */
void checkBandwidthRange(double* min, double* argmin, double *tmp, double order, LPRB* lprb) {
    auto result = std::find(BANDWIDTH_ORDERS, BANDWIDTH_ORDERS + SIZE, order);

    if (result != BANDWIDTH_ORDERS + SIZE) {
        for (int i = 1; i < 10; i++) {      // Here we try a range of bandwiths and compute the residual sum of square of the regression for each
            *tmp = lprb->RSS(order * i);
            if (*tmp < *min) {
                *min = *tmp;
                *argmin = order * i;
            }
        }
    }
}

/**
  * \brief Writes csv files of derivative estimates. See annex.h for more information.
  * \param Yn  Array of doubles containing derivative estimates
  * \param U  Array of doubles containing Uniformly distributed data in [0,1]
  * \param X  Array of doubles containing original data
  * \param kde  KDE object from the kde library (kernel density estimation library) by Tim Nugent (c) 2014
  * \param final LPR_cubic object for containing the Local Polynomial Regressor object from derivative estimates
  * \return None
  */
void writeData(double Yn[], double U[], double X[], KDE* kde, LPR_cubic* final) {
    ofstream file(file_name.c_str());
    ofstream file2(file_name2.c_str());
    if (file) {
        for (int i = 0; i < n - 2; i++) {
           file << X[i] << "," << kde->pdf(X[i]) * Yn[i] << endl;
        }
    }
    else {
        cout << "\nerror opening file for raw data" << endl;
    }

    delete[] Yn;

    if (file2) {
        for (int i = 0; i < NB_POINTS; i++) {
            file2 << ((double)i) / ((double)NB_POINTS) * (RANGE)+(double)START << "," << kde->pdf(((double)i) / ((double)NB_POINTS)) * final->value(((double)i) / ((double)NB_POINTS)) << endl;
        }
    }
    else {
        cout << "\nerror opening file for smoothed data" << endl;
    }
}


