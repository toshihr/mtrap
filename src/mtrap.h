#ifndef _MTRAP_H_
#define _MTRAP_H_

#include <string>
#include "globaloption.h"
#include "utility.h"

void help_line(const std::string& option, const std::string& value, const std::string& explain);
void printHelp();

void setOptionDefaultValue();
int setOption(int argc, char ** argv, int priority = 1);
void setMatrixConfiguration();
void loadFilelist(const std::string& fileName);
void runMSA(const std::string& inFile, const std::string& outFile);
double run(const std::string& inFile, const std::string& outFile);
int main(int argc, char **argv);

#endif
