#ifndef WEIGHTEDMATRIX_H
#define WEIGHTEDMATRIX_H

#include <Eigen/Dense>
#include<Eigen/StdVector>
#include <string>
#include "scorematrix.h"

namespace ScoreMatrix {
// for Weighted Matrices
bool isWeightedMatrixFile(const std::string& fileName);
bool loadWeightedMatrix(const std::string& fileName);
}

#endif // WEIGHTEDMATRIX_H
