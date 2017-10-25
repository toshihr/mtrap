#pragma once
#include <stdint.h> // C99
#include <string>
//#include <vector>
#include <map>
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include "SafeVector.h"
#include "scorematrix.h"
#include "globaloption.h"
#include "utility.h"

#define PART_TQMODE (0)
//#define PART_EPS_ORIGINAL

namespace ScoreMatrix {

// used for path matrix
typedef unsigned int PATH_TYPE;
typedef VVVUI PATH_MATRIX;
extern const PATH_TYPE PATH_INF;
extern const PATH_TYPE PATH_ZERO;

// used for partition function matrix
// floatにすると精度がガタ落ちする
typedef double PART_TYPE;
typedef SafeVector<PART_TYPE> VPART;
typedef SafeVector<VPART> VVPART;
typedef SafeVector<VVPART> VVVPART;
const PART_TYPE PART_INF = DBL_MAX;
const PART_TYPE PART_ZERO = 0.0;
const PART_TYPE PART_ONE = 1.0;

// used for fixed-point arthmetic calculation
extern const double FIXED_SCALE;

extern bool countModeSkipGapGapSite;
extern bool countModeSetZeroGapSite;
extern bool countModeWithoutLowerCase;
extern bool countModeSymmetrizeFreq;
extern bool countModeNoRound;
extern int countModeCluster;

extern bool countModeDirectProb;
extern bool countModeGlobalTQ;

extern double tmatrix[RANK+1][RANK+1][RANK+1][RANK+1]; // 推移距離行列 RANK + GAP

// Transition quantity + Symbol difference
// td((AA11,AA21) -> (AA12,AA22)) = e x TQ() + (1-e) x S()     when (AA11,AA21) -> (AA12,AA22) = (AA*,AA*) -> (AA,AA)
//                                  e x TQ() + (1-e) x Open    when (AA11,AA21) -> (AA12,AA22) = (AA ,*  ) -> (*, AA)
//                                                                                            or (*  ,AA ) -> (AA,* )
//                                  e x TQ() + (1-e) x Ext     when (AA11,AA21) -> (AA12,AA22) = (*  ,AA ) -> (*, AA)
//                                                                                            or (AA ,*  ) -> (AA,* )
extern PATH_TYPE distTQplusDiff_INT[RANK+2][RANK+2][RANK+1][RANK+1];            // Ar-1,Bc-1: RANK+GAP+PRE,  Ar,Bc: RANK+GAP
extern double distTQplusDiff[RANK+2][RANK+2][RANK+1][RANK+1]; // Ar-1,Bc-1: RANK+GAP+PRE,  Ar,Bc: RANK+GAP
extern PART_TYPE distTQplusDiffExt[RANK+2][RANK+2][RANK+1][RANK+1]; // for partition function
extern PART_TYPE distTQplusDiffExtBackward[RANK+2][RANK+2][RANK+1][RANK+1]; // for partition function

typedef Eigen::Matrix<double, RANK, RANK> AminoMatrix;
extern double epsilon;								// T-scoreの重み. (1-eps)*Score + eps*T-Score
extern bool gapWithoutTscore;				// Gapへ対するεの重みをなくす. Substitution Matrixのみにεがかかる

inline const double& getTQ(ScoreMatrix::EncodedResidue Ar_1, ScoreMatrix::EncodedResidue Bc_1, ScoreMatrix::EncodedResidue Ar, ScoreMatrix::EncodedResidue Bc){
    return tmatrix[Ar_1][Bc_1][Ar][Bc];
}

void initTQ();
bool loadTMatrix(const std::string& fileName, bool ignorePostprocess = false); // load transition quantity matrix
void calcDistTQplusDiffTable(const GlobalOption& g);

// --- for the counting mode ---
inline int BLOCKS_INDEX(int n, int col, int row){
    return ( col*n - (col*(col+3))/2 - 1 + row );
}

void countToProb(const VVVVD& countTable, VVVVD& probTable);
void calcTQuantity(const VVVVD& probTable, VVVVD& quantityTable);
void symmetrize(VVD& srcTable);
void symmetrize(VVVVD& srcTable);
int getNextNonGapSitePosition(const std::string& s1, const std::string& s2, int startPos);
void BLOCKS_cluster_seqs(VS& block_seqs, VI& block_cluster_id, VI& block_ncluster);
void loadFrequecyFromBLOCKS(const std::string& fileName, VD& count1Table, VVD& count2Table, VVVVD& count4Table, uint64_t& totAas, uint64_t& totPairs, uint64_t& totTuples);
void loadFrequecy(const std::string& listFile, const std::map<std::string,int>& seqFreqTable, VD& count1Table, VVD& count2Table, VVVVD& count4Table, uint64_t& totAas, uint64_t& totPairs, uint64_t& totTuples);
void generateMatrix(const std::string& listFile, const std::string& seqFreqFile, const std::string& outFile);
void outputBinaryTQ(const std::string& outFile);
}
