// This file defines the LPR_cubic class, used as a local polynomial regressor, using gaussian kernel

#ifndef LPR_CUBIC_H
#define LPR_CUBIC_H
#include <eigen-3.3.7\Eigen\Dense>
#include <vector>


using namespace Eigen;
using namespace std;

class LPR_cubic
{
    public:
        LPR_cubic(double Yn[], MatrixXd Xn, int nn, float hn);
        virtual ~LPR_cubic();
        double value(double x);

    private:
        float h = 1;
        VectorXd Y;
        MatrixXd Xp;
        int n;
        Vector4d B;

};

#endif // LPR_CUBIC_H
