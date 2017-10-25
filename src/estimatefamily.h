#ifndef ESTIMATEFAMILY_H
#define ESTIMATEFAMILY_H

#include <string>
#include "utility.h"
#include "sequences.h"

namespace ScoreMatrix {

std::string findSCOPID(const std::string& s);
std::string estimateGroups(VS& group_list, const Sequences& querySeqs, bool findMostFamous);

}

#endif // ESTIMATEFAMILY_H
