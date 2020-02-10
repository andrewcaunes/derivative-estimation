#include <iostream>
#include <C:\Users\andrew-pc\programmation\C\eigen-3.3.7\eigen-3.3.7\Eigen\Dense>
 using namespace Eigen;

void makeXp(MatrixXd& Xp, double X[], int n){
    Xp.resize(n,4);
    for(int i=0 ; i<n ; i++){
        for(int j=0; j<4 ; j++){
            Xp(i,j) = pow(X[i],j);
        }
    }


//    std::cout << X << std::endl;
//    std::cout << "\n voila Xp : \n" << Xp << " . " << std::endl;

}
