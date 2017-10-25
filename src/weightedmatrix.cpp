#include <cassert>
#include <Eigen/Dense>
#include<Eigen/StdVector>
#include <string>
#include <fstream>
#include <iostream>
#include "scorematrix.h"
#include "weightedmatrix.h"
#include "utility.h"

using namespace std;

namespace ScoreMatrix {
// for weighted matrices
    vector< AminoMatrix, Eigen::aligned_allocator<AminoMatrix> > matrices;
    vector<double> weights;
}


bool ScoreMatrix::isWeightedMatrixFile(const string& fileName)
{
    ifstream ifs(fileName.c_str());
    if(!ifs) return false;
    string buff;
    getline(ifs, buff);
    return (buff == "#WEIGHTED MATRIX#");
}


/*!
 \brief 重み置換行列ファイルを読み込み設定する.

 \param fileName
 \return bool
*/
bool ScoreMatrix::loadWeightedMatrix(const string& fileName)
{
    AminoMatrix locAM;
    double locWeight;

    amatrix.setZero();

    // File format: Comma separated file.
    // [weight],[subsutitution matrix file name]
    ifstream ifs(fileName.c_str());
    string buff;

    // omit the first comment line.
    assert(ifs && "[ScoreMatrix::loadWeightedMatrix] FILE IS NOT FOUND.");
    getline(ifs, buff);

    while(ifs && getline(ifs, buff))
    {
        vector<string> list;
        vector<string> vs;
        mySplit(buff, ", ", vs);
        for(vector<string>::iterator i=vs.begin(); i!=vs.end(); ++i)list.push_back(*i);

        if(vs.size() > 0)
        {
            locWeight = string2binary<double>(vs.at(0));
            const string matrixName = vs.at(1);
            loadMatrix(fixFilePath(matrixName), true, locAM);
            cout << "w=" << locWeight << ":" << matrixName << endl;
            //amatrix += locWeight * locAM;
            amatrix = locAM;
        }
    }

    ifs.close();

    // Calc. amatrix
    // Mode 00: normalized (distance), average
    // Mode 01: normalized (distance), square root
    // Mode 10: non-normalized, average
    // Mode 11: non-normalized, square root
    // 現在値がバグる
    return true;
}
