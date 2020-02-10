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
    gamma_distribution<double> distribution1(2.0,2.0);
    normal_distribution<double> distribution2(0, 0.1);
    for (int i = 0; i < n; i++) {
        X[i] = distribution1(generator);
        Y[i] = pow(X[i], 3) + 2 * pow(X[i], 2) + X[i] - 3 + distribution2(generator);
    }

    KDE* kde = new KDE();
    kde->set_kernel_type(1);
    kde->set_bandwidth_opt_type(1);

    for (int i = 0; i < n; i++) {
        kde->add_data(X[i]);
    }

    const string file_name("C:/Users/andrew-pc/Desktop/data.txt");
    //const string file_name2("C:/Users/andrew-pc/Desktop/data2.txt");
    ofstream file(file_name.c_str());
    //ofstream file2(file_name2.c_str());
    for (int i = 0; i < 200; i++) {
        if (file) {
            file << i * 0.05 << "," << kde->pdf(i * 0.05) << endl;
        }
    }
    file.close();

    delete kde;


    return 0;
}