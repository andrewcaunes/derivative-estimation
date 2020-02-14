#ifndef LPR_CUBIC_H
#define LPR_CUBIC_H
#include <C:\Users\andrew-pc\programmation\C\eigen-3.3.7\eigen-3.3.7\Eigen\Dense>
#include <vector>


using namespace Eigen;
using namespace std;

class LPR_cubic
{
    public:
        LPR_cubic(double Yn[], MatrixXd Xn, int nn, float hn);
        virtual ~LPR_cubic();
        double value(double x);
        Vector4d B;
        vector<double> waits;

   


    protected:

    private:
        float h = 1;
        VectorXd Y;
        MatrixXd Xp;
        int n;
};

#endif // LPR_CUBIC_H
