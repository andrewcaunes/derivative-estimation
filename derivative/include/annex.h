#ifndef ANNEX_H_INCLUDED
#define ANNEX_H_INCLUDED
#include <vector>

void makeXp(MatrixXd& Xp, double X[], int n);
void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove);
bool sort_map(const std::vector<double>& u, const std::vector<double>& y);

#endif // ANNEX_H_INCLUDED
