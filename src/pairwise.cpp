#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <algorithm>
#include <string>
#include <fstream>
#include <map>
#include <set>
#include <cstdlib>
#include <cassert>
#include "scorematrix.h"
#include "tqmatrix.h"
#include "pairwise.h"
#include "tree.h"
#include "globaloption.h"
#include "utility.h"
#include "SparseMatrix.h"
#include "sequences.h"
#include "entropy.h"
#include "distance.h"

using namespace std;
using namespace ScoreMatrix;

AlignerMTRAP::AlignerMTRAP()
    :pName1(NULL),pName2(NULL),pSeq1(NULL),pSeq2(NULL),pGlobalOption(NULL),initialized(false),passDone(false)
{
//    regionEnds1.push_back(3);
//    regionEnds2.push_back(2);
    regionEnds1.push_back(pSeq1->size());
    regionEnds2.push_back(pSeq2->size());
}

AlignerMTRAP::AlignerMTRAP(std::string* pN1, std::string* pS1, std::string* pN2, std::string* pS2, GlobalOption* pG)
    :pName1(pN1),pName2(pN2),pSeq1(pS1),pSeq2(pS2),pGlobalOption(pG),initialized(true),passDone(false)
{
//    regionEnds1.push_back(3);
//    regionEnds2.push_back(2);
    regionEnds1.push_back(pSeq1->size());
    regionEnds2.push_back(pSeq2->size());
}

PROFILE_TYPE AlignerMTRAP::calcPathMatrix() {
    assert(this->initialized);

    const string& A = *pSeq1;
    const string& B = *pSeq2;
    const PATH_TYPE numRow = A.size(); // vertical  <-> sequence 1
    const PATH_TYPE numCol = B.size(); // horizonal <-> sequence 2
    const int D = 0, U = 1, L = 2, N = 3;

    const ScoreMatrix::EncodedResidue PRE = ScoreMatrix::encode('@');
    const ScoreMatrix::EncodedResidue GAP = ScoreMatrix::encode(ScoreMatrix::GAP);
    const PATH_TYPE ZERO = PATH_ZERO;
    const PATH_TYPE INF = PATH_INF;

    // --- initialize ---
    assert(A[0] == PRE && B[0] == PRE && "calcPathMatrix: sequences must have the prefix '@'");
    assign3Matrix(PM, 3, numRow, numCol, ZERO);
    assign3Matrix(DM, 3, numRow, numCol, N);

    int RB = 0, CB = 0; // Begins of the regions Row, Column
    for(int region_i=0; region_i<regionEnds1.size(); region_i++)
    {
        const int RE = regionEnds1[region_i]; // End of the region Row
        const int CE = regionEnds2[region_i]; // End of the region Column

        // ====== FIRST ROW ======
        // --- [0,0] ---
        if(RB == 0 && CB == 0)
        {
            set3Matrix(PM,RB,CB,ZERO,INF,INF);
            set3Matrix(DM,RB,CB, N, N, N);
        }else{
            set3Matrix(PM,RB,CB,PM[D][RB-1][CB-1],INF,INF);
            set3Matrix(DM,RB,CB, D, N, N);
        }

        // --- [0,1]: gap open ---
        {
            const int c=CB+1;
            const PATH_TYPE l = PM[D][RB][CB] + ScoreMatrix::distTQplusDiff_INT[PRE][B[c-1]][GAP][B[c]];
            set3Matrix(PM,RB,c,INF,INF, l );
            set3Matrix(DM,RB,c, N , N , L );
        }

        // --- [0,c] for c=2,...,numCol-1: gap extension ---
        for(int c=CB+2; c<CE; c++)
        {
            const PATH_TYPE l = PM[L][RB][c-1] + ScoreMatrix::distTQplusDiff_INT[GAP][B[c-1]][GAP][B[c]];
            set3Matrix(PM,RB,c,INF,INF, l );
            set3Matrix(DM,RB,c, N , N , L );
        }

        // ====== SECOND ROW ======
        {
            // --- [1,0]: FIRST COLUMN ---
            const int r=RB+1;
            const PATH_TYPE u = PM[D][RB][CB] + ScoreMatrix::distTQplusDiff_INT[A[r-1]][PRE][A[r]][GAP];
            set3Matrix(PM,r,CB,INF, u ,INF);
            set3Matrix(DM,r,CB, N , U , N );

            // --- [1,1]: SECOND COLUMN ---
            {
                const int c=CB+1;
                // PM[0]
                const PATH_TYPE dd = PM[D][r-1][c-1] + ScoreMatrix::distTQplusDiff_INT[A[r-1]][B[c-1]][A[r]][B[c]];
                const pair<PATH_TYPE,int> minD(dd, D );
                // PM[1]
                const PATH_TYPE lu = PM[L][r-1][c  ] + ScoreMatrix::distTQplusDiff_INT[GAP   ][B[c  ]][A[r]][GAP ];
                const pair<PATH_TYPE,int> minU(lu, L );
                // PM[2]
                const PATH_TYPE ul = PM[U][r  ][c-1] + ScoreMatrix::distTQplusDiff_INT[A[r  ]][GAP   ][GAP ][B[c]];
                const pair<PATH_TYPE,int> minL(ul, U );
                // set PM and DM [0,1,2][r][c]
                set3Matrix(PM,r,c,minD.first ,minU.first ,minL.first );
                set3Matrix(DM,r,c,minD.second,minU.second,minL.second);
            }

            // --- [1,c] for c=2,...,numCol-1: rest columns ---
            for(int c=CB+2; c<CE; c++)
            {
                // PM[0]
                const PATH_TYPE ld = PM[L][r-1][c-1] + ScoreMatrix::distTQplusDiff_INT[GAP   ][B[c-1]][A[r]][B[c]];
                const pair<PATH_TYPE,int> minD(ld, L );
                // PM[1]
                const PATH_TYPE lu = PM[L][r-1][c  ] + ScoreMatrix::distTQplusDiff_INT[GAP   ][B[c  ]][A[r]][GAP ];
                const pair<PATH_TYPE,int> minU(lu, L );
                // PM[2]
                const PATH_TYPE dl = PM[D][r  ][c-1] + ScoreMatrix::distTQplusDiff_INT[A[r  ]][B[c-1]][GAP ][B[c]];
#ifndef OMIT_GAP_CROSS_PATH
                const PATH_TYPE ul = PM[U][r  ][c-1] + ScoreMatrix::distTQplusDiff_INT[A[r  ]][GAP   ][GAP ][B[c]];
#endif
                const PATH_TYPE ll = PM[L][r  ][c-1] + ScoreMatrix::distTQplusDiff_INT[GAP   ][B[c-1]][GAP ][B[c]];
#ifndef OMIT_GAP_CROSS_PATH
                const pair<PATH_TYPE,int> minL = choiceMin(dl,ul,ll, D , U , L );
#else
                const pair<PATH_TYPE,int> minL = choiceMin(dl,ll, D , L );
#endif
                // set PM and DM [0,1,2][r][c]
                set3Matrix(PM,r,c,minD.first ,minU.first ,minL.first );
                set3Matrix(DM,r,c,minD.second,minU.second,minL.second);
            }
        }

        // ====== REST ROWS ======
        for(int r=RB+2; r<RE; r++)
        {
            // --- [r,0] for r=2,...,numRow-1: FIRST COLUMN ---
            const PATH_TYPE u = PM[U][r-1][CB] + ScoreMatrix::distTQplusDiff_INT[A[r-1]][GAP][A[r]][GAP];
            set3Matrix(PM,r,CB,INF, u ,INF);
            set3Matrix(DM,r,CB, N , U , N );

            // --- [r,1] for r=2,...,numRow-1: SECOND COLUMN ---
            {
                const int c=CB+1;
                // PM[0]
                const PATH_TYPE ud = PM[U][r-1][c-1] + ScoreMatrix::distTQplusDiff_INT[A[r-1]][GAP   ][A[r]][B[c]];
                const pair<PATH_TYPE,int> minD(ud, U );
                // PM[1]
                const PATH_TYPE du = PM[D][r-1][c  ] + ScoreMatrix::distTQplusDiff_INT[A[r-1]][B[c  ]][A[r]][GAP ];
                const PATH_TYPE uu = PM[U][r-1][c  ] + ScoreMatrix::distTQplusDiff_INT[A[r-1]][GAP   ][A[r]][GAP ];
#ifndef OMIT_GAP_CROSS_PATH
                const PATH_TYPE lu = PM[L][r-1][c  ] + ScoreMatrix::distTQplusDiff_INT[GAP   ][B[c  ]][A[r]][GAP ];
                const pair<PATH_TYPE,int> minU = choiceMin(du,uu,lu, D , U , L );
#else
                const pair<PATH_TYPE,int> minU = choiceMin(du,uu, D , U );
#endif
                // PM[2]
                const PATH_TYPE ul = PM[U][r  ][c-1] + ScoreMatrix::distTQplusDiff_INT[A[r  ]][GAP   ][GAP ][B[c]];
                const pair<PATH_TYPE,int> minL(ul, U );
                // set PM and DM [0,1,2][r][c]
                set3Matrix(PM,r,c,minD.first ,minU.first ,minL.first );
                set3Matrix(DM,r,c,minD.second,minU.second,minL.second);
            }

            // --- [r,c] for c=2,...,numCol, r=2,...,numRow-1: rest columns ---
            for(int c=CB+2; c<CE; c++)
            {
#ifndef OMIT_GAP_CROSS_PATH
                typedef Eigen::Matrix<PATH_TYPE,9,1> Vector9i;
                const Vector9i V =
                        (Vector9i() << PM[D][r-1][c-1], /* dd */
                        PM[U][r-1][c-1], /* ud */
                        PM[L][r-1][c-1], /* ld */
                        PM[D][r-1][c  ], /* du */
                        PM[U][r-1][c  ], /* uu */
                        PM[L][r-1][c  ], /* lu */
                        PM[D][r  ][c-1], /* dl*/
                        PM[U][r  ][c-1], /* ul */
                        PM[L][r  ][c-1]  /* ll */).finished() +
                        (Vector9i() << ScoreMatrix::distTQplusDiff_INT[A[r-1]][B[c-1]][A[r]][B[c]], /* dd */
                        ScoreMatrix::distTQplusDiff_INT[A[r-1]][GAP   ][A[r]][B[c]], /* ud */
                        ScoreMatrix::distTQplusDiff_INT[GAP   ][B[c-1]][A[r]][B[c]], /* ld */
                        ScoreMatrix::distTQplusDiff_INT[A[r-1]][B[c  ]][A[r]][GAP ], /* du */
                        ScoreMatrix::distTQplusDiff_INT[A[r-1]][GAP   ][A[r]][GAP ], /* uu */
                        ScoreMatrix::distTQplusDiff_INT[GAP   ][B[c  ]][A[r]][GAP ], /* lu */
                        ScoreMatrix::distTQplusDiff_INT[A[r  ]][B[c-1]][GAP ][B[c]], /* dl */
                        ScoreMatrix::distTQplusDiff_INT[A[r  ]][GAP   ][GAP ][B[c]], /* ul */
                        ScoreMatrix::distTQplusDiff_INT[GAP   ][B[c-1]][GAP ][B[c]]  /* ll */).finished();

                const pair<PATH_TYPE,int> minD = choiceMin(V[0],V[1],V[2], D , U , L );
                const pair<PATH_TYPE,int> minU = choiceMin(V[3],V[4],V[5], D , U , L );
                const pair<PATH_TYPE,int> minL = choiceMin(V[6],V[7],V[8], D , U , L );
#else
                // if the size is a multiple of 128bits (16bytes), the fixed-size Eigen objects can be vectorized.
                // 32bit x 4 = 128bit, 64bit x 2 = 128bit
                typedef Eigen::Matrix<PATH_TYPE,8,1> Vector8i;
                const Vector8i V (
                            (Vector8i() << PM[D][r-1][c-1], /* dd */
                        PM[U][r-1][c-1], /* ud */
                        PM[L][r-1][c-1], /* ld */
                        PM[D][r-1][c  ], /* du */
                        PM[U][r-1][c  ], /* uu */
                        PM[D][r  ][c-1], /* dl */
                        PM[L][r  ][c-1], /* ll */ 0 ).finished() +
                        (Vector8i() << ScoreMatrix::distTQplusDiff_INT[A[r-1]][B[c-1]][A[r]][B[c]], /* dd */
                        ScoreMatrix::distTQplusDiff_INT[A[r-1]][GAP   ][A[r]][B[c]], /* ud */
                        ScoreMatrix::distTQplusDiff_INT[GAP   ][B[c-1]][A[r]][B[c]], /* ld */
                        ScoreMatrix::distTQplusDiff_INT[A[r-1]][B[c  ]][A[r]][GAP ], /* du */
                        ScoreMatrix::distTQplusDiff_INT[A[r-1]][GAP   ][A[r]][GAP ], /* uu */
                        ScoreMatrix::distTQplusDiff_INT[A[r  ]][B[c-1]][GAP ][B[c]], /* dl */
                        ScoreMatrix::distTQplusDiff_INT[GAP   ][B[c-1]][GAP ][B[c]], /* ll */ 0 ).finished() );

                const pair<PATH_TYPE,int> minD = choiceMin(V[0],V[1],V[2], D , U , L );
                const pair<PATH_TYPE,int> minU = choiceMin(V[3],V[4]     , D , U     );
                const pair<PATH_TYPE,int> minL = choiceMin(V[5]     ,V[6], D     , L );
#endif
                // set PM and DM [0,1,2][r][c]
                set3Matrix(PM,r,c,minD.first ,minU.first ,minL.first );
                set3Matrix(DM,r,c,minD.second,minU.second,minL.second);
            }
        }

        RB = RE;
        CB = CE;
    } // end of region loop

    this->passDone = true;
    return choiceMin(PM[D][numRow-1][numCol-1],PM[U][numRow-1][numCol-1],PM[L][numRow-1][numCol-1]) / ScoreMatrix::FIXED_SCALE;
}

VPROFILE* AlignerMTRAP::genProfile() {
    assert(this->passDone);

    const string& A = *pSeq1;
    const string& B = *pSeq2;
    const int numRow = A.size(); // vertical  <-> sequence 1
    const int numCol = B.size(); // horizonal <-> sequence 2
    const int D = 0, U = 1, L = 2;

    VPROFILE* profilePtr = new VPROFILE(numRow*numCol, PROFILE_ZERO);
    VPROFILE& profile = *profilePtr;

    // --- distance ---
//    string* alignment = this->genAlignment();
//    string* seq1 = genGappedSeq(*(this->pSeq1),*alignment,'X');
//    string* seq2 = genGappedSeq(*(this->pSeq2),*alignment,'Y');
//    const PROFILE_TYPE dist_gid = Distance::calcDistance<PROFILE_TYPE>(*seq1, *seq2, Distance::GID);
//    const PROFILE_TYPE dist_lid = Distance::calcDistance<PROFILE_TYPE>(*seq1, *seq2, Distance::LID);
//    const PROFILE_TYPE dist_eer2 = Distance::calcDistance<PROFILE_TYPE>(*seq1, *seq2, Distance::EER2);
//    delete alignment;
//    delete seq1;
//    delete seq2;

    // --- traceback ---
    // the first row and first column are all 0
    // matchにしか値はないため端に行った時点で終了で問題ない
    int r = numRow - 1;
    int c = numCol - 1;
    pair<PATH_TYPE,int> MM = choiceMin(PM[D][r][c],PM[U][r][c],PM[L][r][c], D , U , L );
    while( r != 0 && c != 0 )
    {
        switch(MM.second) {
        case D:
            // profile は row0 row1 ... rown と格納されている
            profile[r*numCol+c] = 1.0;
            MM.second = DM[D][r][c];
            r--;
            c--;
            break;
        case U:
            MM.second = DM[U][r][c];
            r--;
            break;
        case L:
            MM.second = DM[L][r][c];
            c--;
            break;
        default:
            assert(false && "genProfile: error!!");
        }
    }
    return profilePtr;
}

std::string* AlignerMTRAP::genAlignment() {
    assert(this->passDone);

    const string& A = *pSeq1;
    const string& B = *pSeq2;
    const int numRow = A.size(); // vertical  <-> sequence 1
    const int numCol = B.size(); // horizonal <-> sequence 2
    const int D = 0, U = 1, L = 2;

    string* alignment = new string(); assert(alignment);

    int r = numRow - 1;
    int c = numCol - 1;

    pair<PATH_TYPE,int> MM = choiceMin(PM[D][r][c],PM[U][r][c],PM[L][r][c], D , U , L );

    // --- traceback ---
    while( r != 0 || c != 0 )
    {
        switch(MM.second) {
        case D:
            alignment->push_back('B');
            MM.second = DM[D][r][c];
            r--;
            c--;
            break;
        case U:
            alignment->push_back('X');
            MM.second = DM[U][r][c];
            r--;
            break;
        case L:
            alignment->push_back('Y');
            MM.second = DM[L][r][c];
            c--;
            break;
        default:
            assert(false && "genProfile: error!!");
        }
    }

    reverse(alignment->begin(), alignment->end());
    return alignment;
}


void AlignerMTRAP::showPathMatrix()
{
    assert(this->passDone);
    PATH_MATRIX* pM = &PM;

    PATH_MATRIX& M = *pM;

    const string& A = *pSeq1;
    const string& B = *pSeq2;
    const int numRow = A.size(); // vertical  <-> sequence 1
    const int numCol = B.size(); // horizonal <-> sequence 2

    putLog("path matrix DUL:");
    for(int r=0; r<numRow; r++)
    {
        stringstream line;

        for(int c=0; c<numCol; c++)
        {
            line << "(";
            for(int i=0; i<3; i++)
            {
                if(PATH_INF == M[i][r][c])
                {
                    line << left << setw(10) << "INF";
                }else if(PATH_ZERO == M[i][r][c])
                {
                    line << left << setw(10) << "---";
                }else{
//                    line << scientific << left << setw(10) << setprecision(2) << (M[i][r][c] / ScoreMatrix::FIXED_SCALE);
                    line << scientific << left << setw(10) << setprecision(2) << M[i][r][c];
                }
            }
            line << ") ";
        }
        if(M == PM)
        {
            char MARK[4] = {'D','U','L','-'};
            line << "|";
            for(int c=0; c<numCol; c++)
            {
                for(int i=0; i<3; i++)
                {
                    line << left << setw(1) << MARK[DM[i][r][c]];
                }
                line << ' ';
            }
        }
        putLog(line.str());
    }
    putLog("");

}

void AlignerMTRAP::showPartMatrix(ScoreMatrix::VVVPART* pM)
{
    assert(pM != NULL);

    VVVPART& M = *pM;

    const string& A = *pSeq1;
    const string& B = *pSeq2;
    const int numRow = A.size(); // vertical  <-> sequence 1
    const int numCol = B.size(); // horizonal <-> sequence 2

    const PART_TYPE ZERO = PART_ZERO;

    putLog("partition function matrix DUL:");
    for(int r=0; r<numRow; r++)
    {
        stringstream line;

        for(int c=0; c<numCol; c++)
        {
            line << "(";
            for(int i=0; i<3; i++)
            {
                if(PART_INF - M[i][r][c] < PART_ONE/1000000000)
                {
                    line << left << setw(10) << "INF";
                }else if(fabs(ZERO - M[i][r][c]) < PART_ONE/1000000000)
                {
                    line << left << setw(10) << "---";
                }else{
                    line << scientific << left << setw(10) << setprecision(2) << M[i][r][c];
                }
            }
            line << ") ";
        }
        putLog(line.str());
    }
    putLog("");

}


/* The correspondence among AlignerMTRAP, Probalign and Miyazawa's paper
 * Z[D][r,c] <-> Zm <-> Zij
 * Z[U][r,c] <-> Zf <-> ZBij
 * Z[L][r,c] <-> Ze <-> ZAij
 */
ScoreMatrix::PART_TYPE AlignerMTRAP::calcPartForward(ScoreMatrix::VVVPART& Z) {
    assert(this->initialized);

    const string& A = *pSeq1;
    const string& B = *pSeq2;
    const int numRow = A.size(); // vertical  <-> sequence 1
    const int numCol = B.size(); // horizonal <-> sequence 2
    const int D = 0, U = 1, L = 2;

    const ScoreMatrix::EncodedResidue PRE = ScoreMatrix::encode('@');
    const ScoreMatrix::EncodedResidue GAP = ScoreMatrix::encode(ScoreMatrix::GAP);
    const PART_TYPE ZERO = PART_ZERO;
    const PART_TYPE ONE = PART_ONE;

    // --- initialize ---
    assert(A[0] == PRE && B[0] == PRE && "calcPartForward: sequences must have the prefix '@'");
    assign3Matrix(Z, 3, numRow, numCol, ZERO);

    // ====== FIRST ROW ======
    // --- [0,0] ---
    set3Matrix(Z,0,0,ONE,ZERO,ZERO);

    // --- [0,1]: gap open ---
    {
        const int c=1;
        const PART_TYPE l = Z[D][0][c-1] * ScoreMatrix::distTQplusDiffExt[PRE][B[c-1]][GAP][B[c]];
#ifdef PAIRWISE_PART_PROBALIGN_STYLE
        set3Matrix(Z,0,c, ZERO, ZERO, l );
#else   // miyazawa p.1009 (32) line 1 with i = 0
        set3Matrix(Z,0,c, l , ZERO, l );
#endif
    }

    // --- [0,c] for c=2,...,numCol-1: gap extension ---
    for(int c=2; c<numCol; c++)
    {
#ifdef PAIRWISE_PART_PROBALIGN_STYLE
        const PART_TYPE l = Z[L][0][c-1] * ScoreMatrix::distTQplusDiffExt[GAP][B[c-1]][GAP][B[c]];
        set3Matrix(Z,0,c, ZERO, ZERO, l );
#else   // miyazawa p.1009 (32) line 2 with i = 0
        const PART_TYPE l = Z[D][0][c-1] * ScoreMatrix::distTQplusDiffExt[PRE][B[c-1]][GAP][B[c]]
                          + Z[L][0][c-1] * ScoreMatrix::distTQplusDiffExt[GAP][B[c-1]][GAP][B[c]];
        set3Matrix(Z,0,c, l , ZERO, l );
#endif
    }

    // ====== SECOND ROW ======
    // these code is same as REST ROW if the PAIRWISE_PART_PROBALIGN_STYLE is defined
    {
        // --- [1,0]: FIRST COLUMN ---
        const int r=1;
        const PART_TYPE u = Z[D][r-1][0] * ScoreMatrix::distTQplusDiffExt[A[r-1]][PRE][A[r]][GAP];
#ifdef PAIRWISE_PART_PROBALIGN_STYLE
        set3Matrix(Z,r,0, ZERO, u , ZERO);
#else   // miyazawa p.1009 (32) line 3 with j = 0
        set3Matrix(Z,r,0,  u , u , ZERO);
#endif

        // --- [1,1]: HERE IS DIFFERENT FROM MTRAP ---
#ifndef PAIRWISE_PART_PROBALIGN_STYLE
        assert(false && "calcPartForward: [1,1] these code should be written!");
        // but the result becomes same one.
        // because, for example, Z[U,L][0][0] = 0, so the result becomes same one.
        // multiply 0 is the key.""
#endif

        // --- [1,c] for c=1,...,numCol-1: rest columns ---
        for(int c=1; c<numCol; c++)
        {
            // PM[0]
            const PART_TYPE dd = Z[D][r-1][c-1] * ScoreMatrix::distTQplusDiffExt[A[r-1]][B[c-1]][A[r]][B[c]];
            const PART_TYPE ud = Z[U][r-1][c-1] * ScoreMatrix::distTQplusDiffExt[A[r-1]][GAP   ][A[r]][B[c]];
            const PART_TYPE ld = Z[L][r-1][c-1] * ScoreMatrix::distTQplusDiffExt[GAP   ][B[c-1]][A[r]][B[c]];
            const PART_TYPE sumD = dd + ud + ld;
            // PM[1]
            const PART_TYPE du = Z[D][r-1][c  ] * ScoreMatrix::distTQplusDiffExt[A[r-1]][B[c  ]][A[r]][GAP ];
            const PART_TYPE uu = Z[U][r-1][c  ] * ScoreMatrix::distTQplusDiffExt[A[r-1]][GAP   ][A[r]][GAP ];
#ifndef OMIT_GAP_CROSS_PATH
            const PART_TYPE lu = Z[L][r-1][c  ] * ScoreMatrix::distTQplusDiffExt[A[GAP]][B[c-1]][A[r]][GAP ];
            const PART_TYPE sumU = du + uu + lu;
#else
            const PART_TYPE sumU = du + uu;
#endif
            // PM[2]
            const PART_TYPE dl = Z[D][r  ][c-1] * ScoreMatrix::distTQplusDiffExt[A[r  ]][B[c-1]][GAP ][B[c]];
#ifndef OMIT_GAP_CROSS_PATH
            const PART_TYPE ul = Z[U][r  ][c-1] * ScoreMatrix::distTQplusDiffExt[A[r  ]][GAP   ][GAP ][B[c]];
#endif
            const PART_TYPE ll = Z[L][r  ][c-1] * ScoreMatrix::distTQplusDiffExt[GAP   ][B[c-1]][GAP ][B[c]];
#ifndef OMIT_GAP_CROSS_PATH
            const PART_TYPE sumL = dl + ul + ll;
#else
            const PART_TYPE sumL = dl + ll;
#endif
            // set PM and DM [0,1,2][r][c]
            set3Matrix(Z,r,c, sumD , sumU , sumL);
        }
    }

    // ====== REST ROWS ======
     for(int r=2; r<numRow; r++)
    {
         // --- [r,0] for r=2,...,numRow-1: FIRST COLUMN ---
#ifdef PAIRWISE_PART_PROBALIGN_STYLE
        const PART_TYPE u = Z[U][r-1][0] * ScoreMatrix::distTQplusDiffExt[A[r-1]][GAP][A[r]][GAP];
        set3Matrix(Z,r,0, ZERO, u , ZERO);
#else   // miyazawa style
        const PART_TYPE u = Z[D][r-1][0] * ScoreMatrix::distTQplusDiffExt[A[r-1]][PRE][A[r]][GAP]
                          + Z[U][r-1][0] * ScoreMatrix::distTQplusDiffExt[A[r-1]][GAP][A[r]][GAP];
        set3Matrix(Z,r,0,  u , u , ZERO);
#endif

        // --- [r,1] for r=2,...,numRow-1: HERE IS DIFFERENT FROM MTRAP ---

        // --- [r,c] for c=1,...,numCol, r=2,...,numRow-1: rest columns ---
        for(int c=1; c<numCol; c++)
        {
            typedef Eigen::Array<PART_TYPE,8,1> Array8f;
            const Array8f V (
                        (Array8f() << Z[D][r-1][c-1], /* dd */
                                      Z[U][r-1][c-1], /* ud */
                                      Z[L][r-1][c-1], /* ld */
                                      Z[D][r-1][c  ], /* du */
                                      Z[U][r-1][c  ], /* uu */
                                      Z[D][r  ][c-1], /* dl */
                                      Z[L][r  ][c-1], /* ll */ PART_ZERO ).finished() *
                        (Array8f() << ScoreMatrix::distTQplusDiffExt[A[r-1]][B[c-1]][A[r]][B[c]], /* dd */
                                      ScoreMatrix::distTQplusDiffExt[A[r-1]][GAP   ][A[r]][B[c]], /* ud */
                                      ScoreMatrix::distTQplusDiffExt[GAP   ][B[c-1]][A[r]][B[c]], /* ld */
                                      ScoreMatrix::distTQplusDiffExt[A[r-1]][B[c  ]][A[r]][GAP ], /* du */
                                      ScoreMatrix::distTQplusDiffExt[A[r-1]][GAP   ][A[r]][GAP ], /* uu */
                                      ScoreMatrix::distTQplusDiffExt[A[r  ]][B[c-1]][GAP ][B[c]], /* dl */
                                      ScoreMatrix::distTQplusDiffExt[GAP   ][B[c-1]][GAP ][B[c]], /* ll */ PART_ZERO ).finished() );

            const PART_TYPE sumD = V.matrix().block<3,1>(0,0).sum(); // du + uu + lu
            const PART_TYPE sumU = V.matrix().block<2,1>(3,0).sum(); // du + uu       (du + uu + lu)
            const PART_TYPE sumL = V.matrix().block<2,1>(5,0).sum(); // dl + ll       (dl + ul + ll)

//            // PM[0]
//            const PART_TYPE dd = Z[D][r-1][c-1] * ScoreMatrix::distTQplusDiffExt[A[r-1]][B[c-1]][A[r]][B[c]];
//            const PART_TYPE ud = Z[U][r-1][c-1] * ScoreMatrix::distTQplusDiffExt[A[r-1]][GAP   ][A[r]][B[c]];
//            const PART_TYPE ld = Z[L][r-1][c-1] * ScoreMatrix::distTQplusDiffExt[GAP   ][B[c-1]][A[r]][B[c]];
//            const PART_TYPE sumD = dd + ud + ld;
//            // PM[1]
//            const PART_TYPE du = Z[D][r-1][c  ] * ScoreMatrix::distTQplusDiffExt[A[r-1]][B[c  ]][A[r]][GAP ];
//            const PART_TYPE uu = Z[U][r-1][c  ] * ScoreMatrix::distTQplusDiffExt[A[r-1]][GAP   ][A[r]][GAP ];
//#ifndef OMIT_GAP_CROSS_PATH
//            const PART_TYPE lu = Z[L][r-1][c  ] * ScoreMatrix::distTQplusDiffExt[A[GAP]][B[c-1]][A[r]][GAP ];
//            const PART_TYPE sumU = du + uu + lu;
//#else
//            const PART_TYPE sumU = du + uu;
//#endif
//            // PM[2]
//            const PART_TYPE dl = Z[D][r  ][c-1] * ScoreMatrix::distTQplusDiffExt[A[r  ]][B[c-1]][GAP ][B[c]];
//#ifndef OMIT_GAP_CROSS_PATH
//            const PART_TYPE ul = Z[U][r  ][c-1] * ScoreMatrix::distTQplusDiffExt[A[r  ]][GAP   ][GAP ][B[c]];
//#endif
//            const PART_TYPE ll = Z[L][r  ][c-1] * ScoreMatrix::distTQplusDiffExt[GAP   ][B[c-1]][GAP ][B[c]];
//#ifndef OMIT_GAP_CROSS_PATH
//            const PART_TYPE sumL = dl + ul + ll;
//#else
//            const PART_TYPE sumL = dl + ll;
//#endif
            // set PM and DM [0,1,2][r][c]
            set3Matrix(Z,r,c, sumD, sumU , sumL);
        }
    }

    this->passDone = true;
    return Z[D][numRow-1][numCol-1] + Z[U][numRow-1][numCol-1] + Z[L][numRow-1][numCol-1];
}

VPROFILE* AlignerMTRAP::genPartProfile() {
    assert(this->initialized);

    VVVPART forwardZ;
    const PART_TYPE denoZ = calcPartForward(forwardZ);

    // appending '@' right side is necessary for using ScoreMatrix::distTQplusDiffExt
    const string A = pSeq1->substr(1, pSeq1->size()-1) + ScoreMatrix::encode('@');
    const string B = pSeq2->substr(1, pSeq2->size()-1) + ScoreMatrix::encode('@');
    const int numRow = A.size(); // vertical  <-> sequence 1
    const int numCol = B.size(); // horizonal <-> sequence 2
    const int D = 0, U = 1, L = 2;

    const ScoreMatrix::EncodedResidue PRE = ScoreMatrix::encode('@');
    const ScoreMatrix::EncodedResidue GAP = ScoreMatrix::encode(ScoreMatrix::GAP);
    const PART_TYPE ZERO = PART_ZERO;
    const PART_TYPE ONE = PART_ONE;

    // --- initialize ---
    assert(A[numRow-1] == PRE && B[numCol-1] == PRE && "calcPartForward: sequences must have the prefix '@'");
    VVVPART Z;
    assign3Matrix(Z, 3, numRow, numCol, ZERO);

    VPROFILE *profilePtr = new VPROFILE(numRow * numCol);
    VPROFILE & profile = *profilePtr;
    VPROFILE::iterator ptr = profile.begin();

    // ====== FIRST ROW ======
    // --- [0,0] ---
    set3Matrix(Z,numRow-1,numCol-1,ONE,ZERO,ZERO);

    // --- [0,1]: gap open ---
    {
        const int c=numCol-2, cc = c;
        const PART_TYPE l = Z[D][numRow-1][c+1] * ScoreMatrix::distTQplusDiffExtBackward[PRE][B[cc+1]][GAP][B[cc]];
#ifdef PAIRWISE_PART_PROBALIGN_STYLE
        set3Matrix(Z,numRow-1,c, ZERO, ZERO, l );
#else   // miyazawa p.1009 (32) line 1 with i = 0
        set3Matrix(Z,numRow-1,c, l , ZERO, l );
#endif
    }

    // --- [0,c] for c=2,...,numCol-1: gap extension ---
    for(int c=numCol-3; c>=0; c--)
    {
        const int cc = c;
#ifdef PAIRWISE_PART_PROBALIGN_STYLE
        const PART_TYPE l = Z[L][numRow-1][c+1] * ScoreMatrix::distTQplusDiffExtBackward[GAP][B[cc+1]][GAP][B[cc]];
        set3Matrix(Z,numRow-1,c, ZERO, ZERO, l );
#else   // miyazawa p.1009 (32) line 2 with i = 0
        const PART_TYPE l = Z[D][numRow-1][c+1] * ScoreMatrix::distTQplusDiffExtBackward[PRE][B[cc+1]][GAP][B[cc]]
                          + Z[L][numRow-1][c+1] * ScoreMatrix::distTQplusDiffExtBackward[GAP][B[cc+1]][GAP][B[cc]];
        set3Matrix(Z,numRow-1,c, l , ZERO, l );
#endif
    }

    // ====== SECOND ROW ======
    {
        // --- [1,0]: FIRST COLUMN ---
        const int r=numRow-2, rr = r;
        const PART_TYPE u = Z[D][r+1][numCol-1] * ScoreMatrix::distTQplusDiffExtBackward[A[rr+1]][PRE][A[rr]][GAP];
#ifdef PAIRWISE_PART_PROBALIGN_STYLE
        set3Matrix(Z,r,numCol-1, ZERO, u , ZERO);
#else   // miyazawa p.1009 (32) line 3 with j = 0
        set3Matrix(Z,r,numCol-1,  u , u , ZERO);
#endif

        // --- [1,1]: HERE IS DIFFERENT FROM MTRAP ---

        // --- [1,c] for c=1,...,numCol-1: rest columns ---
        for(int c=numCol-2; c>=0; c--)
        {
            const int cc = c;
            // PM[0]
            const PART_TYPE dd = Z[D][r+1][c+1] * ScoreMatrix::distTQplusDiffExtBackward[A[rr+1]][B[cc+1]][A[rr]][B[cc]];
            const PART_TYPE ud = Z[U][r+1][c+1] * ScoreMatrix::distTQplusDiffExtBackward[A[rr+1]][GAP    ][A[rr]][B[cc]];
            const PART_TYPE ld = Z[L][r+1][c+1] * ScoreMatrix::distTQplusDiffExtBackward[GAP    ][B[cc+1]][A[rr]][B[cc]];
            const PART_TYPE sumD = dd + ud + ld;
            // PM[1]
            const PART_TYPE du = Z[D][r+1][c  ] * ScoreMatrix::distTQplusDiffExtBackward[A[rr+1]][B[cc  ]][A[rr]][GAP  ];
            const PART_TYPE uu = Z[U][r+1][c  ] * ScoreMatrix::distTQplusDiffExtBackward[A[rr+1]][GAP    ][A[rr]][GAP  ];
#ifndef OMIT_GAP_CROSS_PATH
            const PART_TYPE lu = Z[L][r+1][c  ] * ScoreMatrix::distTQplusDiffExtBackward[GAP    ][B[cc  ]][A[rr]][GAP  ];
            const PART_TYPE sumU = du + uu + lu;
#else
            const PART_TYPE sumU = du + uu;
#endif
            // PM[2]
            const PART_TYPE dl = Z[D][r  ][c+1] * ScoreMatrix::distTQplusDiffExtBackward[A[rr  ]][B[cc+1]][GAP  ][B[cc]];
#ifndef OMIT_GAP_CROSS_PATH
            const PART_TYPE ul = Z[U][r  ][c+1] * ScoreMatrix::distTQplusDiffExtBackward[A[rr  ]][GAP    ][GAP  ][B[cc]];
#endif
            const PART_TYPE ll = Z[L][r  ][c+1] * ScoreMatrix::distTQplusDiffExtBackward[GAP    ][B[cc+1]][GAP  ][B[cc]];
#ifndef OMIT_GAP_CROSS_PATH
            const PART_TYPE sumL = dl + ul + ll;
#else
            const PART_TYPE sumL = dl + ll;
#endif
            // set PM and DM [0,1,2][r][c]
            set3Matrix(Z,r,c, sumD , sumU , sumL);

            // --- calculate P(Ar ~ Bc) ---
            /*{ // [1] --- probalign (code) style ---
                const PART_TYPE s = forwardZ[D][r+1][c+1] * Z[D][r][c];
                const PART_TYPE scorez = ScoreMatrix::distTQplusDiffExtBackward[A[rr+1]][B[cc+1]][A[rr]][B[cc]];
                const PART_TYPE prob = s / (scorez * denoZ);
                //cerr << r << "," << c << ":" << s << "/(" << scorez << " x " << denoZ << ")=" << prob << endl;
                if(0.001 <= prob) ptr[rr*numCol+cc] = min(1.0, prob);
            }*/
            /*{ // [2] --- this code generates same result as probalign if epsilon == 0 ---
                const PART_TYPE s = forwardZ[D][r+1][c+1] * (Z[D][r+1][c+1]+Z[U][r+1][c+1]+Z[L][r+1][c+1]);
                const PART_TYPE prob = s / denoZ;
                if(0.001 <= prob) ptr[(r+1)*numCol+(c+1)] = min(1.0, prob);
            }*/
            { // [2] --- this code generates same result as probalign if epsilon == 0 ---
                const PART_TYPE s = forwardZ[D][r+1][c+1] * Z[D][r][c] / ScoreMatrix::distTQplusDiffExtBackward[PRE][PRE][A[rr]][B[cc]];
                const PART_TYPE prob = s / denoZ;
                if(0.001 <= prob) ptr[(r+1)*numCol+(c+1)] = min(static_cast<PART_TYPE>(1.0), prob);
            }
            /*{ // [2] --- this code generates same result as probalign if epsilon == 0 ---
                const PART_TYPE s = (forwardZ[D][r][c] + forwardZ[U][r][c] + forwardZ[L][r][c]) * (Z[D][r+1][c+1]+Z[U][r+1][c+1]+Z[L][r+1][c+1]);
                const PART_TYPE scorez = ScoreMatrix::distTQplusDiffExtBackward[PRE][PRE][A[rr]][B[cc]];
                const PART_TYPE prob = s * scorez / denoZ;
                if(0.001 <= prob) ptr[(r+1)*numCol+(c+1)] = min(1.0, prob);
            }*/
        }
    }

    // ====== REST ROWS ======
    for(int r=numRow-3; r>=0; r--)
    {
        // --- [r,0] for r=2,...,numRow-1: FIRST COLUMN ---
        const int rr = r;
#ifdef PAIRWISE_PART_PROBALIGN_STYLE
        const PART_TYPE u = Z[U][r+1][numCol-1] * ScoreMatrix::distTQplusDiffExtBackward[A[rr+1]][GAP][A[rr]][GAP];
        set3Matrix(Z,r,numCol-1, ZERO, u , ZERO);
#else   // miyazawa style
        const PART_TYPE u = Z[D][r+1][numCol-1] * ScoreMatrix::distTQplusDiffExtBackward[A[rr+1]][PRE][A[rr]][GAP]
                          + Z[U][r+1][numCol-1] * ScoreMatrix::distTQplusDiffExtBackward[A[rr+1]][GAP][A[rr]][GAP];
        set3Matrix(Z,r,numCol-1,  u , u , ZERO);
#endif

        // --- [r,1] for r=2,...,numRow-1: HERE IS DIFFERENT FROM MTRAP ---

        // --- [r,c] for c=1,...,numCol, r=2,...,numRow-1: rest columns ---
        for(int c=numCol-2; c>=0; c--)
        {
            const int cc = c;
            // PM[0]
            const PART_TYPE dd = Z[D][r+1][c+1] * ScoreMatrix::distTQplusDiffExtBackward[A[rr+1]][B[cc+1]][A[rr]][B[cc]];
            const PART_TYPE ud = Z[U][r+1][c+1] * ScoreMatrix::distTQplusDiffExtBackward[A[rr+1]][GAP    ][A[rr]][B[cc]];
            const PART_TYPE ld = Z[L][r+1][c+1] * ScoreMatrix::distTQplusDiffExtBackward[GAP    ][B[cc+1]][A[rr]][B[cc]];
            const PART_TYPE sumD = dd + ud + ld;
            // PM[1]
            const PART_TYPE du = Z[D][r+1][c  ] * ScoreMatrix::distTQplusDiffExtBackward[A[rr+1]][B[cc  ]][A[rr]][GAP  ];
            const PART_TYPE uu = Z[U][r+1][c  ] * ScoreMatrix::distTQplusDiffExtBackward[A[rr+1]][GAP    ][A[rr]][GAP  ];
#ifndef OMIT_GAP_CROSS_PATH
            const PART_TYPE lu = Z[L][r+1][c  ] * ScoreMatrix::distTQplusDiffExtBackward[GAP    ][B[cc  ]][A[rr]][GAP  ];
            const PART_TYPE sumU = du + uu + lu;
#else
            const PART_TYPE sumU = du + uu;
#endif
            // PM[2]
            const PART_TYPE dl = Z[D][r  ][c+1] * ScoreMatrix::distTQplusDiffExtBackward[A[rr  ]][B[cc+1]][GAP  ][B[cc]];
#ifndef OMIT_GAP_CROSS_PATH
            const PART_TYPE ul = Z[U][r  ][c+1] * ScoreMatrix::distTQplusDiffExtBackward[A[rr  ]][GAP    ][GAP  ][B[cc]];
#endif
            const PART_TYPE ll = Z[L][r  ][c+1] * ScoreMatrix::distTQplusDiffExtBackward[GAP    ][B[cc+1]][GAP  ][B[cc]];
#ifndef OMIT_GAP_CROSS_PATH
            const PART_TYPE sumL = dl + ul + ll;
#else
            const PART_TYPE sumL = dl + ll;
#endif
            // set PM and DM [0,1,2][r][c]
            set3Matrix(Z,r,c, sumD, sumU, sumL);

            // --- calculate P(Ar ~ Bc) ---
            /*{ // [1] --- probalign (code) style ---
                const PART_TYPE s = forwardZ[D][r+1][c+1] * Z[D][r][c];
                const PART_TYPE scorez = ScoreMatrix::distTQplusDiffExtBackward[A[rr+1]][B[cc+1]][A[rr]][B[cc]];
                const PART_TYPE prob = s / (scorez * denoZ);
                //cerr << r << "," << c << ":" << s << "/(" << scorez << " x " << denoZ << ")=" << prob << endl;
                if(0.001 <= prob) ptr[rr*numCol+cc] = min(1.0, prob);
            }*/
            /*{ // [2] --- this code generates same result as probalign if epsilon == 0 ---
                const PART_TYPE s = forwardZ[D][r+1][c+1] * (Z[D][r+1][c+1]+Z[U][r+1][c+1]+Z[L][r+1][c+1]);
                const PART_TYPE prob = s / denoZ;
                if(0.001 <= prob) ptr[(r+1)*numCol+(c+1)] = min(1.0, prob);
            }*/
            { // [2] --- this code generates same result as probalign if epsilon == 0 ---
                const PART_TYPE s = forwardZ[D][r+1][c+1] * Z[D][r][c] / ScoreMatrix::distTQplusDiffExtBackward[PRE][PRE][A[rr]][B[cc]];
                const PROFILE_TYPE prob = static_cast<PROFILE_TYPE>(s / denoZ);
                if(PROFILE_ONE/1000 <= prob) ptr[(r+1)*numCol+(c+1)] = min(PROFILE_ONE, prob);
            }
            /*{ // [2] --- this code generates same result as probalign if epsilon == 0 ---
                const PART_TYPE s = (forwardZ[D][r][c] + forwardZ[U][r][c] + forwardZ[L][r][c]) * (Z[D][r+1][c+1]+Z[U][r+1][c+1]+Z[L][r+1][c+1]);
                const PART_TYPE scorez = ScoreMatrix::distTQplusDiffExtBackward[PRE][PRE][A[rr]][B[cc]];
                const PART_TYPE prob = s * scorez / denoZ;
                if(0.001 <= prob) ptr[(r+1)*numCol+(c+1)] = min(1.0, prob);
            }*/
        }
    }
#ifdef DEBUG
    putLog("forward  Z=" + DS(denoZ));
    putLog("backward Z=" + DS(Z[D][0][0] + Z[U][0][0] + Z[L][0][0]));

    putLog("forward Z:");
    showPartMatrix(&forwardZ);
    putLog("backward Z:");
    showPartMatrix(&Z);
#endif
    return profilePtr;
}

// MSAProbs::MSA::DoIterativeRefinement()相当
void doIterativeRef(const VVpSM& sparseProfiles, Sequences* & alignment, VPROFILE* pSeqsWeights, const GlobalOption& gOption, int seedid) {
    const int numSeqs = alignment->seqs.size();
    VI indexList(numSeqs), g1(numSeqs/2), g2(numSeqs - g1.size());

    srand(seedid);

    for(int i=0; i<numSeqs; i++)indexList[i] = i;
    random_shuffle(indexList.begin(), indexList.end());

    const int g1size = g1.size();
    int i = 0;
    for(; i<g1size; i++) g1[i] = indexList[i];
    for(; i<numSeqs;i++) g2[i - g1size] = indexList[i];

    assert(!g1.empty());
    assert(!g2.empty());

    // split the alignment to two alignments
    Sequences* align1 = alignment->extract(g1);
    Sequences* align2 = alignment->extract(g2);
    delete alignment;

    // realign
    alignment = mergeAlignments(align1, align2, sparseProfiles, pSeqsWeights);

    delete align1;
    delete align2;
}


// MSAProbs::MSA::ComputeFinalAlignment()相当
Sequences* genFinalAlignment(tr_node* pGuideTree, const Sequences* pSeqs, const VVpSM& sparseProfiles, VPROFILE* pSeqsWeights, const GlobalOption& gOption) {
    Sequences* alignment = doProgressiveAlignment(pGuideTree, pSeqs, sparseProfiles, pSeqsWeights);
    // --- treat an ordering ---
    // TODO:

    // --- iterative refinements ---
    if(pSeqs->seqs.size() > 2)
    {
        for(int i=0; i<gOption.numIterativeRefReps; i++)
        {
            putLog("iterative " + IS(i+1) + ".");
            doIterativeRef(sparseProfiles, alignment, pSeqsWeights, gOption, i);
        }
    }
    // --- sort ---
    if(gOption.sortByInputOrder) alignment->sort();

    return alignment;
}


// MSAProbs::MSA::ProcessTree()相当
// 再帰するごとにpTreeは変化するが，pSeqs等は固定
Sequences* doProgressiveAlignment(tr_node* pTree, const Sequences* pSeqs, const VVpSM& sparseProfiles, VPROFILE* pSeqsWeights){
    Sequences* result;
//putLog("progressive alignment: node id=" + binary2string<int>(pTree->id));

    if(pTree->id < 0)
    {
        // node
        Sequences* alignLeft = doProgressiveAlignment(pTree->pLeft, pSeqs, sparseProfiles, pSeqsWeights);
        Sequences* alignRight = doProgressiveAlignment(pTree->pRight, pSeqs, sparseProfiles, pSeqsWeights);
        assert(alignLeft);
        assert(alignRight);
//putLog("merge the two alignments");
        result = mergeAlignments(alignLeft, alignRight, sparseProfiles, pSeqsWeights);
//putLog("done");
        assert(result);

        delete alignLeft;
        delete alignRight;
    }else{
        // leaf
        result = new Sequences();
        assert(result);
        result->addSequence(pSeqs->names[pTree->id], pSeqs->seqs[pTree->id], pTree->id);
    }

    return result;
}


// probcons::ProbabilisticModel.h::BuildPosterior()相当
// M[i,j] = sum{s in align1} sum{s in align2} f(s,t,i,j)
// where
// f(s,t,i,j) = P(s[i'] <-> t[j']) when s[i'] is a symbol in the ith site of align1 and
//                                      t[j'] is a symbol in the jth site of align2
//              0                  otherwise
VPROFILE* buildProfile(VPROFILE* pSeqsWeights, Sequences* align1, Sequences* align2, const VVpSM& sparseProfiles){
//    cout << "build profile:" << endl;
//    align1->output(cout);
//    align2->output(cout);

    const PROFILE_TYPE cutoff = 0.0;
    // size of the alignment 1 and 2 includes the size of prefix '@'
    // this is the same size of the new profile
    const int numNewRow = align1->seqs[0]->size();
    const int numNewCol = align2->seqs[0]->size();

    VPROFILE* profilePtr = new VPROFILE(numNewRow * numNewCol, PROFILE_ZERO); assert(profilePtr);
    VPROFILE& profile = *profilePtr;

    // compute the total sum of all weights
//putLog("compute the total sum of all weights");
    PROFILE_TYPE totalWeights = 0.0;
    for(size_t i=0; i<align1->seqs.size(); i++)
    {
        const int first = align1->index[i];
        const PROFILE_TYPE w1 = (*pSeqsWeights)[first];
        for(size_t j=0; j<align2->seqs.size(); j++)
        {
            const int second = align2->index[j];
            const PROFILE_TYPE w2 = (*pSeqsWeights)[second];

            totalWeights += w1 * w2;
        }
    }

    // generate the weighted profile
//putLog("generate the weighted profile");
    // for each s in align1
    for(size_t i=0; i<align1->seqs.size(); i++)
    {
        const int first = align1->index[i];
        const PROFILE_TYPE w1 = (*pSeqsWeights)[first];
        VI* mapping1 = genNonGapIndexList(*(align1->seqs[i]));
        for(size_t j=0; j<align2->seqs.size(); j++)
        {
            const int second = align2->index[j];
            const PROFILE_TYPE w2 = (*pSeqsWeights)[second];
            VI* mapping2 = genNonGapIndexList(*(align2->seqs[j]));

            PROFILE_TYPE w = w1 * w2 / totalWeights;

            pair<int,int> ij = myMinmax(first,second);
//cout << "ij=" << ij.first << "," << ij.second << endl;
            SparseMatrix* matrix = sparseProfiles[ij.first][ij.second];
            assert(matrix && "ij is wrong!!");

            int seq1Length = matrix->GetSeq1Length(); // numRow - 1
            int seq2Length = matrix->GetSeq2Length(); // numCol - 1
            // here, not use the trans matrix because of the calculation costs
            if(ij.first == first)
            {
                for(int ii=1; ii<=seq1Length; ii++)
                {
                    SafeVector<PID>::iterator row = matrix->GetRowPtr(ii);
                    const int base = (*mapping1)[ii] * numNewCol;
                    const int sparseRow = matrix->GetRowSize(ii);
                    // add in all relevant values
                    for(int jj=0; jj<sparseRow; jj++)
                    {
                        profile[base + (*mapping2)[row[jj].first]] += w * row[jj].second;
                    }
                    // subtract cutoff
                    for(int jj=0; jj<seq2Length; jj++)
                    {
                        profile[base + (*mapping2)[jj]] -= w * cutoff;
                    }
                }
            }else{
                for(int jj=1; jj<=seq1Length; jj++)
                {
                    SafeVector<PID>::iterator row = matrix->GetRowPtr(jj);
                    const int base = (*mapping2)[jj];
                    const int sparseRow = matrix->GetRowSize(jj);
                    // add in all relevant values
                    for(int ii=0; ii<sparseRow; ii++)
                    {
                        profile[base + (*mapping1)[row[ii].first] * numNewCol] += w * row[ii].second;
                    }
                    // subtract cutoff
                    for(int ii=0; ii<seq2Length; ii++)
                    {
                        profile[base + (*mapping1)[ii] * numNewCol] -= w * cutoff;
                    }
                }
            }
            delete mapping2;
        }
        delete mapping1;
    }


    return profilePtr;
}

// probcons::MSA::AlignAlignments()相当
Sequences* mergeAlignments(Sequences* align1, Sequences* align2, const VVpSM& sparseProfiles, SafeVector<PROFILE_TYPE>* pSeqsWeights){
    // make a profile
//putLog("make a profile from alignment");
//putLog("group1:");
//align1->output(cout);
//putLog("group2:");
//align2->output(cout);
    VPROFILE* profile = buildProfile(pSeqsWeights, align1, align2, sparseProfiles);
//SparseMatrix sm = SparseMatrix(align1->seqs[0]->size(), align2->seqs[0]->size(), *profile);
//sm.Print(cout);
    // do profile alignment
//putLog("do profile alignment");
    pair<string*, PROFILE_TYPE> alignment = doProfileAlignment(align1->seqs[0]->size(), align2->seqs[0]->size(), *profile);
    // finalize
    delete profile;
    // generate the final alignment
//putLog("generate the alignment");
    Sequences* result = new Sequences();
    for(size_t i=0; i<align1->seqs.size(); i++) result->addSequence(align1->names[i], pString(genGappedSeq(*align1->seqs[i], *alignment.first, 'X')), align1->index[i]);
    for(size_t i=0; i<align2->seqs.size(); i++) result->addSequence(align2->names[i], pString(genGappedSeq(*align2->seqs[i], *alignment.first, 'Y')), align2->index[i]);

    delete alignment.first;

//putLog("merged:");
//result->output(cout);
    return result;
}

// probcons::ProbabilisticModel.h::ComputeAlignment()相当
// return the best alignment string and the best alignment score
pair<std::string*, PROFILE_TYPE> doProfileAlignment(int seq1Length, int seq2Length, const VPROFILE& profile) {
    const int numRow = seq1Length; // with the length of prefix '@'
    const int numCol = seq2Length; // with the length of prefix '@'
    PROFILE_TYPE* twoRows = new PROFILE_TYPE[numCol*2]; assert(twoRows);
    PROFILE_TYPE* oldRow = twoRows;
    PROFILE_TYPE* newRow = twoRows + numCol;

    char* tracebackMatrix = new char[numRow * numCol]; assert(tracebackMatrix);
    char* tracebackPtr = tracebackMatrix;

    assert(profile.size() == numRow * numCol && "the size of profile is wrong!!");

    // initialization
    for(int c=0; c<numCol; c++)
    {
        oldRow[c] = 0.0;
        *(tracebackPtr++) = 'L'; // Left
    }

    // omit the first row
    VPROFILE::const_iterator profilePtr = profile.begin() + numCol;

    // fill the matrix (dynamic programming)
    for(int r=1; r<numRow; r++)
    {
        // initialize left column
        newRow[0] = 0.0;
        profilePtr++;
        *(tracebackPtr++) = 'U'; // Upper

        // fill in rest of row
        for(int c=1; c<numCol; c++)
        {
            // current node (r,c) is newRow[c] so
            // Diag  (r-1,c-1) = oldRow[c-1]
            // Left  (r,c-1) = newRow[c-1]
            // Upper (r-1,c) = oldRow[c]
            // TODO: dの計算はベクトル化し高速化できる(Eigen3の出番)
            const PROFILE_TYPE d = *(profilePtr++) + oldRow[c-1];
            const PROFILE_TYPE l = newRow[c-1];
            const PROFILE_TYPE u = oldRow[c];
            pair<PROFILE_TYPE,char> best = choiceMax(d,u,l,'D','U','L');
            newRow[c] = best.first;
            *(tracebackPtr++) = best.second;
        }

        // swap rows
        PROFILE_TYPE* temp = oldRow;
        oldRow = newRow;
        newRow = temp;
    }

    // store alignment score
    PROFILE_TYPE total = oldRow[seq2Length-1];
    delete [] twoRows;

    // debug: output the tracebackMatrix
    /*cout << "numRow,numCol=" << numRow << "," << numCol << endl;
    for(int i=0; i<numRow; i++)
    {
        for(int j=0; j<numCol; j++)
        {
            cout << tracebackMatrix[i*numCol+j];
        }
        cout << endl;
    }*/

    // compute traceback
    string* alignment = new string(); assert(alignment);

    int r = numRow - 1;
    int c = numCol - 1;
    while( r != 0 || c != 0)
    {
        const char ch = tracebackMatrix[r * numCol + c];
        switch(ch){
        case 'L': c--; alignment->push_back('Y'); break; // seq1にギャップ
        case 'U': r--; alignment->push_back('X'); break; // seq2にギャップ
        case 'D': c--; r--; alignment->push_back('B'); break; // seq1にギャップ
        default: assert(false);
        }
    }

    // finalize
    delete [] tracebackMatrix;

    reverse(alignment->begin(), alignment->end());
    return make_pair(alignment, total);
}

// MSA::DoRelaxation()相当
// TODO: アルゴリズムの理解
void doConsistencyTrans(const VPROFILE& seqsWeights, const Sequences* pSeqs, VVpSM& sparseProfiles) {
    const int numSeqs = pSeqs->seqs.size();

    VVpSM newSparseProfiles(numSeqs, VpSM(numSeqs, static_cast<pSparseMatrix>(NULL)));

#ifdef _OPENMP
    int pairIndex;
#pragma omp parallel for private(pairIndex) default(shared) schedule(dynamic)
    for(pairIndex=0; pairIndex<globalOption.numProfilePairs; pairIndex++)
    {
        int i = globalOption.profilePairs[pairIndex].first;
        int j = globalOption.profilePairs[pairIndex].second;
        const PROFILE_TYPE wi = seqsWeights[i];
        const PROFILE_TYPE wj = seqsWeights[j];
#else
    for(int i=0; i<numSeqs; i++)
    {
        const PROFILE_TYPE wi = seqsWeights[i];
        for(int j=i+1; j<numSeqs; j++)
        {
            const PROFILE_TYPE wj = seqsWeights[j];
#endif
            string& seq1 = *(pSeqs->seqs[i]);
            string& seq2 = *(pSeqs->seqs[j]);

            VPROFILE* profilePtr = sparseProfiles[i][j]->GetPosterior(); assert(profilePtr);
            VPROFILE& profile = *profilePtr;

            const int numRow = seq1.size();
            const int numCol = seq2.size();

            // --- consistency calculation ---
            // STEP1: only calculate the summation where z = x and z = y
            const PROFILE_TYPE w1 = wi*wi*wj + wi*wj*wj;
            PROFILE_TYPE sumW = w1;
            for(int k=0; k<numRow*numCol; k++) profile[k] *= w1;
            // STEP2: calculate the other parts
            for(int k=0; k<numSeqs; k++) if(k != i && k != j)
            {
                const PROFILE_TYPE wk = seqsWeights[k];
                const PROFILE_TYPE w2 = wi*wj*wk;
                sumW += w2;
                if(k < i)
                {
                    consistencyZXZY(w2, sparseProfiles[k][i], sparseProfiles[k][j], profile);
                }else if(k > i && k < j){
                    consistencyXZZY(w2, sparseProfiles[i][k], sparseProfiles[k][j], profile);
                }else{
                    // TODO: 転置するコストを考えると特化consistency関数を作成すべき
                    SparseMatrix* t = sparseProfiles[j][k]->ComputeTranspose();
                    consistencyXZZY(w2, sparseProfiles[i][k], t, profile);
                    delete t;
                }
            } // for k
            // STEP3:
            for(int k=0; k<numRow*numCol; k++) profile[k] /= sumW;
            // STEP4: mask out
            SparseMatrix* matXY = sparseProfiles[i][j];
            for(int y=0; y<numCol; y++) profile[y] = 0.0;
            for(int x=1; x<numRow; x++)
            {
                SafeVector<PID>::iterator XYptr = matXY->GetRowPtr(x);
                SafeVector<PID>::iterator XYend = XYptr + matXY->GetRowSize(x);
                VPROFILE::iterator base = profile.begin() + x * numCol;
                int curr = 0;
                for( ;XYptr != XYend; curr++, ++XYptr)
                {
                    // set 0 before the first filled column
                    for( ;curr < XYptr->first; curr++) base[curr] = 0.0;
                }
                // set 0 after last column
                for( ; curr<numCol; curr++) base[curr] = 0.0;
            }
            // --- save the new profile ---
            newSparseProfiles[i][j] = new SparseMatrix(seq1.size()-1, seq2.size()-1, profile);

            // --- finalize ---
            delete profilePtr;
#ifndef _OPENMP
        } // for j
#endif
    } // for i


    // copy
    for (int i=0; i<numSeqs-1; i++){
        for (int j=i+1; j<numSeqs; j++){
            delete sparseProfiles[i][j];
            sparseProfiles[i][j] = newSparseProfiles[i][j];
        }
    }
}

// MSA::Relax()相当
void consistencyXZZY(PROFILE_TYPE weight, SparseMatrix* matXZ, SparseMatrix* matZY, VPROFILE& profile)
{
    assert(matXZ);
    assert(matZY);

    const int numRow = matXZ->GetSeq1Length() + 1; // X
    const int numCol = matZY->GetSeq2Length() + 1; // Y
    assert(matXZ->GetSeq2Length() == matZY->GetSeq1Length());

    for(int i=1; i<numRow; i++)
    { // every x[i]
        SafeVector<PID>::iterator XZptr = matXZ->GetRowPtr(i);
        SafeVector<PID>::iterator XZend = XZptr + matXZ->GetRowSize(i);

        VPROFILE::iterator base = profile.begin() + i * numCol;

        while (XZptr != XZend)
        { // x[i] <-> z[k]
            SafeVector<PID>::iterator ZYptr = matZY->GetRowPtr(XZptr->first);
            SafeVector<PID>::iterator ZYend = ZYptr + matZY->GetRowSize(XZptr->first);
            const PROFILE_TYPE valXZ = XZptr->second;

            while(ZYptr != ZYend)
            { // z[k] <-> y[j]
                base[ZYptr->first] += weight * valXZ * ZYptr->second;
                ZYptr++;
            }
            XZptr++;
        }
    }
}

// MSA::Relax1()相当
void consistencyZXZY(PROFILE_TYPE weight, SparseMatrix* matZX, SparseMatrix* matZY, VPROFILE& profile)
{
    assert(matZX);
    assert(matZY);

    const int numRow = matZX->GetSeq1Length() + 1; // Z
    const int numCol = matZY->GetSeq2Length() + 1; // Y

    for(int k=1; k<numRow; k++)
    {
        SafeVector<PID>::iterator ZXptr = matZX->GetRowPtr(k);
        SafeVector<PID>::iterator ZXend = ZXptr + matZX->GetRowSize(k);

        while (ZXptr != ZXend)
        { // x[i] <-> z[k]
            SafeVector<PID>::iterator ZYptr = matZY->GetRowPtr(k);
            SafeVector<PID>::iterator ZYend = ZYptr + matZY->GetRowSize(k);
            const PROFILE_TYPE valZX = ZXptr->second;
            VPROFILE::iterator base = profile.begin() + ZXptr->first * numCol;

            while(ZYptr != ZYend)
            { // z[k] <-> y[j]
                base[ZYptr->first] += weight * valZX * ZYptr->second;
                ZYptr++;
            }
            ZXptr++;
        }
    }
}

void showProfile(VPROFILE* profile, int numRow, int numCol)
{
    for(int r=0; r<numRow; r++)
    {
        stringstream line;
        for(int c=0; c<numCol; c++)
        {
            line << left << setw(8) << setprecision(3) << (*profile)[r*numCol+c];
        }
        putLog(line.str());
    }
}

