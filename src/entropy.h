#ifndef _ENTROPY_H_
#define _ENTROPY_H_

#include <map>
#include <vector>
#include <cmath>
#include "SafeVector.h"
#include "utility.h"

namespace Entropy {

typedef float T;
typedef SafeVector<T> VT;
typedef SafeVector<VT> VVT;
const T ZERO = 0.0;
const T ONE = 1.0f;

const int PROB_TABLE_SIZE = 0x7f + 1;

bool calcComplexProb(VVT& complexProb, const std::string& s1, const std::string& s2);
bool calcProb(VT& prob, const std::string s);
T getEntropy(const VT& prob);
T getMutualEntropy(const VT& prob1, const VT& prob2, const VVT& complexProb);
T getEER(const VT& prob1, const VT& prob2, const VVT& complexProb);
T getEER2(const VT& prob1, const VT& prob2, const VVT& complexProb);
T getECD(const VT& prob1, const VT& prob2, const VVT& complexProb);

}

#endif
