#ifndef _PRIMARYLIBRARY_H_
#define _PRIMARYLIBRARY_H_

#include <string>
#include "distance.h"
#include "utility.h"
#include "mtrap.h"

void makePrimaryLibrary(const std::string& databaseFile, const std::string& outFile);
int getResiduePairsAndNumOfIdentity(std::vector<std::pair<int,int> >& indexPair, const pString& s1, const pString& s2);

#endif
