// derivative.cpp : Defines the entry point for the application.
//

#include "derivative.h"


using Eigen::MatrixXd;
using namespace std;

int main()
{
    // SMOOTHING |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    //cout << "\n Making data : __________________________" << endl;
    //default_random_engine generator;
    //gamma_distribution<double> distribution1(2.0,2.0);
    //normal_distribution<double> distribution2(0, 0.1);

    ////
    ////    for(int i = 0; i<10 ; i++){
    ////        cout << distribution(generator) << endl;
    ////    }


    //    //TEST LPRB
    ////const string file_name("C:/Users/andrew-pc/Desktop/data.txt");
    ////const string file_name2("C:/Users/andrew-pc/Desktop/data2.txt");
    ////ofstream file(file_name.c_str());
    ////ofstream file2(file_name2.c_str());
    //const int n(200);
    //double X[n];
    //double Y[n];//TEST
    //for (int i = 0; i < 200; i++) {
    //    X[i] = distribution1(generator);
    //    Y[i] = pow(X[i], 3) + 2 * pow(X[i], 2) + X[i] - 3 + distribution2(generator);
    //    //TEST cout << "\n X[" << i << "]=" << X[i] << ", Y[" << i << "]=" << Y[i] << endl;
    //}



    //MatrixXd Xp;
    //makeXp(Xp, X, n);
    //float h(1);
    //cout << "\n End Making data : __________________________" << endl;
    //cout << "\n Fitting to data : __________________________" << endl;


    //LPR_cubic o1 = LPR_cubic(Y, Xp, n, h);

    //cout << "\nEnd fitting to data : __________________________" << endl;

    ///*TEST double Yp[n];
    //double sum(0);

    //for (int i = 0; i < 200; i++) {
    //    Yp[i] = o1.value(X[i]);
    //    sum += pow(Y[i] - Yp[i], 2);
    //}
    //cout << "\n RS : " << sum << endl;
    //int i1(0);
    //cout << "\n i1 should be "<< pow(i1,3)+2*pow(i1,2)+i1-3 <<" but is really " << o1.value(i1) <<endl;
    //cout << "\n B: " << o1.B << endl;*/

    // COMPUTE k ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    const int n(200);
    double X[n];
    double Y[n];

    cout << "\n Making data : __________________________" << endl;
    default_random_engine generator;
    uniform_real_distribution<double> distribution1(0,10.0);
    normal_distribution<double> distribution2(0, 0.1);
    for (int i = 0; i < n; i++) {
        X[i] = distribution1(generator);
        Y[i] = sin(X[i]);
    }

    KDE* kde = new KDE();
    kde->set_kernel_type(1);
    kde->set_bandwidth_opt_type(1);

    for (int i = 0; i < n; i++) {
        kde->add_data(X[i]);
    }

    //const string file_name("C:/Users/andrew-pc/Desktop/data.txt");
    ////const string file_name2("C:/Users/andrew-pc/Desktop/data2.txt");
    //ofstream file(file_name.c_str());
    ////ofstream file2(file_name2.c_str());
    //for (int i = 0; i < 200; i++) {
    //    if (file) {
    //        file << i * 0.05 << "," << kde->pdf(i * 0.05) << endl;
    //    }
    //}
    //file.close();

    // CALCUL B
    double B(0);

    double U[n];
    for (int i = 0; i < n; i++) {
        U[i] = kde->cdf(X[i]);
    }
    // SORTING U AND Y
    vector<vector<double>> map;
    for (int i = 0; i < n; i++) {
        map.push_back(vector<double>(0));
        map[i].push_back(U[i]);
        map[i].push_back(Y[i]);
    };
    sort(map.begin(), map.end(), sort_map);   
    for (int i = 0; i < n; i++) {
        U[i] = map[i][0];
        Y[i] = map[i][1];
    }

    // END SORTING

    for (int i = 0; i < n; i++) {
        cout << " U[i] " << U[i] << " Y[i] " << Y[i] << endl;

    };
    

    MatrixXd Up;
    makeXp(Up, U, n);
    LPR_cubic approx = LPR_cubic(Y,Up,n,0.01);
    

    


    double* max = new double;
    *max = 0;
    double* tmp = new double;
    *tmp = 0;
    for (int i = 0; i < 30; i++) {
        approx.value((double)i * 0.03);
        *tmp += abs(approx.B(3) * 6 * (double)i*0.03 + 2*approx.B(2));
        /*if (*tmp > *max){
            *max = *tmp;
        }*/
    }
    //B = *tmp / 30;
    B = 10;

    cout << "\n B : " << B << endl;

    delete max;
    delete tmp;

    // FIN CALCUL B

    // CALCUL SIGMA
    double sigma2(0);
    for (int i = 0; i < n - 2; i++) {
        sigma2 += pow((0.809 * Y[i] - 0.5 * Y[i + 1] - 0.309 * Y[i + 2]),2);
    }
    sigma2 *= (1 / ((double)n - 2));

    cout << "sigma2 : " << sigma2 << endl;

    // FIN CALCUL SIGMA

    // CALCUL k
    int k;
    double* min = new double;
    *min = 10000;
    int* argmin = new int;
    *argmin = n;
    tmp = new double;

    for (int i = 1; i < (int)(((double)n - 1) / 2); i++) {
        *tmp = pow(B * 3 * i * ((double)i + 1) / 4 / ((double)n + 1) / (2 * (double)i + 1), 2);
        *tmp += 3 * sigma2 * ((double)n + 1) * ((double)n + 1) / i / ((double)i + 1) / (2 * (double)i + 1);
        if (*tmp < *min){
            *min = *tmp;
            *argmin = i;
        }
    }
    k = *argmin;
    cout << "\k : " << k << endl;

    delete argmin;
    delete min;
    delete tmp;
    delete kde;

    // FIN CALCUL K

    double weights[n][n];
    double* sum = new double;
    int* ki = new int;
    *ki = 0;

    *sum = 0;
    for (int i = 1; i < n; i++) {
        
        for (int j = 0; j < n; j++) {
            if (i <= n - k && i >= k + 1) {
                
                *sum = 0;
                for (int l = 1; l <= k; l++) {
                    *sum += pow(U[i + l - 1] - U[i - l -1],2);
                }
                weights[i][j] = pow((U[i + j - 1] - U[i - j - 1]), 2) / *sum;
            }
            else if (i<k + 1) {
                *ki = i - 1;
                *sum = 0;
                if (j <= *ki) {
                    for (int l = 1; l <= *ki; l++) {
                        *sum += pow(U[i + l - 1] - U[i - l - 1], 2);
                    }
                    for (int l = *ki+1; l <= k; l++) {
                        *sum += pow(U[i + l - 1] - U[i - 1], 2);
                    }
                    weights[i][j] = pow(U[i + j - 1] - U[i - j - 1], 2) / *sum;
                }else if (j > *ki) {
                    for (int l = 1; l <= *ki; l++) {
                        *sum += pow(U[i + l - 1] - U[i - l - 1], 2);
                    }
                    for (int l = *ki+1; l <= k; l++) {
                        *sum += pow(U[i + l - 1] - U[i - 1], 2);
                    }
                    weights[i][j] = pow(U[i + j - 1] - U[i - 1], 2) / *sum;
                }
            }else if (i > n - k) {
                *ki = n - i;
                *sum = 0;
                //cout << "\n i : " << i << " j : " << j << endl;
                if (j <= *ki) {
                    for (int l = 1; l <= *ki; l++) {
                        //cout << "\n i : " << i << ", l : " << l << endl;
                        *sum += pow(U[i + l - 1] - U[i - l - 1], 2);
                    }
                    for (int l = *ki+1; l <= k; l++) {
                        *sum += pow(U[i - 1] - U[i-l - 1], 2);
                    }
                    weights[i][j] = pow(U[i + j - 1] - U[i - j - 1], 2) / *sum;
                }
                else if (j > * ki) {
                    for (int l = 1; l <= *ki; l++) {
                        *sum += pow(U[i + l - 1] - U[i - l - 1], 2);
                    }
                    for (int l = *ki+1; l <= k; l++) {
                        *sum += pow(U[i - 1] - U[i-l - 1], 2);
                    }
                    weights[i][j] = pow(U[i - 1] - U[i-j - 1], 2) / *sum;
                }
            }
        }
    }
    
    double Yd[n];

    for (int i = 1; i < n; i++) {
        Yd[i-1] = 0;
        if (i <= n - k && i >= k + 1) {
            for (int j = 1; j <= k; j++) {
                Yd[i-1] += weights[i-1][j-1] * (Y[i + j - 1] - Y[i - j - 1]) / (U[i + j - 1] - U[i - j - 1]);
                
            }
        }
        else if (i < k + 1) {
            *ki = i - 1;
            if (i > 1) {
                for (int j = 1; j <= *ki; j++) {
                    if (i == 1) {
                    }
                    Yd[i - 1] += weights[i - 1][j - 1] * (Y[i + j - 1] - Y[i - j - 1]) / (U[i + j - 1] - U[i - j - 1]);

                }
            }
            for (int j = *ki+1; j <= k; j++) {
                Yd[i-1] += weights[i-1][j-1] * (Y[i + j - 1] - Y[i - 1]) / (U[i + j - 1] - U[i - 1]);
                
            }
            
        }
        else if (i > n - k) {
            *ki = n-i;
            for (int j = 1; j <= *ki; j++) {
                Yd[i-1] += weights[i-1][j-1] * (Y[i + j - 1] - Y[i - j - 1]) / (U[i + j - 1] - U[i - j - 1]);
            }
            for (int j = *ki+1; j <= k; j++) {
                Yd[i-1] += weights[i-1][j-1] * (Y[i - 1] - Y[i-j - 1]) / (U[i - 1] - U[i-j - 1]);
            }
        }
    }

    removeRow(Up, 49);
    removeRow(Up, 0);
    double Yn[n - 2];
    for (int i = 0; i < n - 2; i++) {
        Yn[i] = Yd[i + 1];
    }
    cout << Up;
    LPR_cubic final = LPR_cubic(Yn, Up, n-2, 1);
    cout << "\n val : " << final.value(U[47]);

    const string file_name("C:/Users/andrew-pc/Desktop/data.txt");
    ofstream file(file_name.c_str());
    const string file_name2("C:/Users/andrew-pc/Desktop/data2.txt");
    ofstream file2(file_name2.c_str());

    /*for (int i = 0; i < n - 2; i++) {
        cout << "weight"<<i<<" : " << final.waits[i] <<endl;
        file2 << i << "," << final.waits[i] << endl;
    }*/

    
    for (int i = 0; i < n-2; i++) {
        if (file) {
            file << U[i] << "," << Yn[i] << endl;
        }
    }

    for (int i = 0; i < 1000; i++) {
        if (file2) {
            file2 << 0.001*i << "," << final.value(i*0.001) << endl;
        }
    }


    delete sum;

    return 0;
}