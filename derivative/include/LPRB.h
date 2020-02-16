// This file defines the LPRB class, used as a local polynomial regressor with specific bandwidth, using Binomial kernel.


#ifndef LPRB_H
#define LPRB_H
#include <C:\Users\andrew-pc\programmation\C\eigen-3.3.7\eigen-3.3.7\Eigen\Dense>


using namespace Eigen;
class LPRB
{
    public:
        LPRB(double Yn[], MatrixXd Xn, int nn, float hn);
        virtual ~LPRB();
        double value(double x);
        double RSS(float hn);

    private:
        float h = 1;
        VectorXd Y;
        MatrixXd Xp;
        Vector4d B;
        int n;
};

#endif // LPRB_H
