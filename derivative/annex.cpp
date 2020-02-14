#include <iostream>
#include <C:\Users\andrew-pc\programmation\C\eigen-3.3.7\eigen-3.3.7\Eigen\Dense>
#include<vector>
 using namespace Eigen;

 void makeXp(MatrixXd& Xp, double X[], int n) {
     Xp.resize(n, 4);
     for (int i = 0; i < n; i++) {
         for (int j = 0; j < 4; j++) {
             Xp(i, j) = pow(X[i], j);
         }
     }
 }

bool sort_map(const std::vector<double> &u,const std::vector<double> &y){
    return u[0] < y[0];
}

void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove)
{
    unsigned int numRows = matrix.rows() - 1;
    unsigned int numCols = matrix.cols();

    if (rowToRemove < numRows)
        matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) = matrix.block(rowToRemove + 1, 0, numRows - rowToRemove, numCols);

    matrix.conservativeResize(numRows, numCols);
}


//    std::cout << X << std::endl;
//    std::cout << "\n voila Xp : \n" << Xp << " . " << std::endl;