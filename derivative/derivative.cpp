// derivative.cpp : Defines the entry point for the application.
//

#include "derivative.h"


using Eigen::MatrixXd;
using namespace std;

int main()
{
    cout << "\n Making data : __________________________" << endl;
    default_random_engine generator;
    gamma_distribution<double> distribution1(2.0,2.0);
    normal_distribution<double> distribution2(0, 0.1);

    //
    //    for(int i = 0; i<10 ; i++){
    //        cout << distribution(generator) << endl;
    //    }


        //TEST LPRB
    const string file_name("C:/Users/andrew-pc/Desktop/data.txt");
    const string file_name2("C:/Users/andrew-pc/Desktop/data2.txt");
    ofstream file(file_name.c_str());
    ofstream file2(file_name2.c_str());
    const int n(200);
    double X[n];
    double Y[n];//TEST
    for (int i = 0; i < 200; i++) {
        X[i] = distribution1(generator);
        Y[i] = pow(X[i], 3) + 2 * pow(X[i], 2) + X[i] - 3 + distribution2(generator);
        cout << "\n X[" << i << "]=" << X[i] << ", Y[" << i << "]=" << Y[i] << endl;
        if (file) {
            file << X[i] << "," << Y[i] << endl;
        }
        else {
            cout << "file error" << endl;
        }
    }



    MatrixXd Xp;
    makeXp(Xp, X, n);
    float h(1);
    cout << "\n End Making data : __________________________" << endl;
    cout << "\n Fitting to data : __________________________" << endl;


    /*LPR_cubic o1 = LPR_cubic(Y, Xp, n, h);*/

    cout << "\nEnd fitting to data : __________________________" << endl;

    float htest(1);
    double Yp[n];
    LPRB o2 = LPRB(Y, Xp, n, htest);
    double sum(0);
    for (int i = 0; i < 200; i++) {
        Yp[i] = o2.value(X[i]);
        file2 << X[i] << "," << Y[i] << endl;
        cout << "\ndiff pour i="<<i<<" : " << Y[i]-Yp[i] << endl;
        sum += pow(Y[i] - Yp[i], 2);
    }
    cout << "\n sum : " << sum << endl;
    int i1(0);
    cout << "\n i1 should be "<< pow(i1,3)+2*pow(i1,2)+i1-3 <<" but is really " << o2.value(i1) <<endl;
    cout << "\n B: " << o2.B << endl;
    int i2(1);
    cout << "\n i2 should be " << pow(i2, 3) + 2 * pow(i2, 2) + i2 - 3 << " but is really " << o2.value(i2) << endl;
    int i3(2);
    cout << "\n i3 should be " << pow(i3, 3) + 2 * pow(i3, 2) + i3 - 3 << " but is really " << o2.value(i3) << endl;
    int i4(3);
    cout << "\n i4 should be " << pow(i4, 3) + 2 * pow(i4, 2) + i4 - 3 << " but is really " << o2.value(i4) << endl;
    int i5(4);
    cout << "\n i5 should be " << pow(i5, 3) + 2 * pow(i5, 2) + i5 - 3 << " but is really " << o2.value(i5) << endl;
    cout << "\n RSS is " << o2.RSS(htest) << endl;



    /*LPRB o = LPRB(Y, Xp, n, 1);
    for (int i = 1; i < 100; i++) {
        float j(i*0.1);
        o.setH(j);
        cout << "\n rss for " << i << " is " << o.RSS(j) << endl;
    }*/
    //cout << "\n Voila RSS pour h= " << h << " : " << o1.RSS(h) << " ; \n";
    //h = 2;
    //LPRB o2 = LPRB(Y, Xp, n, h);
    //cout << "\n Voila RSS pour h= " << h << " : " << o2.RSS(h) << " ; \n";
    //h = 1.5;
    //LPRB o3 = LPRB(Y, Xp, n, h);
    //cout << "\n Voila RSS pour h= " << h << " : " << o3.RSS(h) << " ; \n";
    //h = 0.1;
    //LPRB o4 = LPRB(Y, Xp, n, h);
    //cout << "\n Voila RSS pour h= " << h << " : " << o4.RSS(h) << " ; \n";
    ////    cout << o.value(6); TEST*
    //o4.setH(3);
    //cout << "\n Voila RSS pour h= " << h << " : " << o4.RSS(h) << " ; \n";



    return 0;
}