#pragma once

#include <string>
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include "tree.h"
#include "globaloption.h"
#include "sequences.h"
#include "scorematrix.h"
#include "tqmatrix.h"
#include "SparseMatrix.h"

//#define DEBUG_LOG

// if not defined, all three path directions are used.
#define OMIT_GAP_CROSS_PATH

// if not defined, using miyazawa paper's style will be used.
// This is same as MTRAP code, so should be defined
#define PAIRWISE_PART_PROBALIGN_STYLE

//TODO: random refinementを一様分布以外にしたがって行うとよくなる？
//TODO: EERなUPGMAに対応する

// AlignerMTRAP オブジェクトが開放されるまで配列は生きている必要がある
class AlignerMTRAP
{
private:
    // --- data ---
    // sequences must have prefix '@' for easy access by using the same index as PM,DM
    std::string* pName1;
    std::string* pName2;
    std::string* pSeq1; // Ar: @,A1,A2,...,AnumRow
    std::string* pSeq2; // Bc: @,B1,B2,...,BnumCol
    GlobalOption* pGlobalOption;
    // region information: all regions are connected with each neighbors
    VI regionEnds1; // the set of right edges (STL style). last value = AnumRow
    VI regionEnds2; // the set of right edges (STL style). last value = BnumCol

    // --- private variables ---
    // structure of array (SOA) style
    ScoreMatrix::PATH_MATRIX PM;
    VVVI DM;

    bool initialized;
    bool passDone;

    // --- partition function ---
    ScoreMatrix::PART_TYPE calcPartForward(ScoreMatrix::VVVPART& Z);

public:
    AlignerMTRAP();
    AlignerMTRAP(std::string* pN1, std::string* pS1, std::string* pN2, std::string* pS2, GlobalOption* pG);

    // --- STEP 1 ---
    PROFILE_TYPE calcPathMatrix();
    // --- STEP 2 ---
    VPROFILE* genProfile();
    std::string* genAlignment();

    // --- partition function ---
    VPROFILE* genPartProfile();

    // --- DEBUG ---
    void showPathMatrix();
    void showPartMatrix(ScoreMatrix::VVVPART* pM);
};


template<typename T>
inline void set3Matrix(SafeVector<SafeVector<SafeVector< T > > > & M, int i, int j, const T& v1, const T& v2, const T& v3) {
    M[0][i][j] = v1;
    M[1][i][j] = v2;
    M[2][i][j] = v3;
}

template<typename T>
inline void assign3Matrix(SafeVector<SafeVector<SafeVector< T > > > & M, int i1, int i2, int i3, const T& v) {
    M.assign(i1, SafeVector<SafeVector< T > >(i2, SafeVector< T >(i3, v) ) );
}


// --- profile alignments ---
// sequence length は prefix の長さを含む物
void doIterativeRef(const VVpSM& sparseProfiles, Sequences* & alignment, VPROFILE* pSeqsWeights, const GlobalOption& gOption, int seedid);
Sequences* genFinalAlignment(tr_node* pGuideTree, const Sequences* pSeqs, const VVpSM& sparseProfiles, VPROFILE* pSeqsWeights, const GlobalOption& gOption);
Sequences* doProgressiveAlignment(tr_node* pTree, const Sequences* pSeqs, const VVpSM& sparseProfiles, VPROFILE* pSeqsWeights);
VPROFILE* buildProfile(SafeVector<PROFILE_TYPE>* pSeqsWeights, Sequences* align1, Sequences* align2, const VVpSM& sparseProfiles);
Sequences* mergeAlignments(Sequences* align1, Sequences* align2, const VVpSM& sparseProfiles, VPROFILE* pSeqsWeights);
pair<std::string*, PROFILE_TYPE> doProfileAlignment(int seq1Length, int seq2Length, const VPROFILE& profile);
void doConsistencyTrans(const VPROFILE& seqsWeights, const Sequences* pSeqs, VVpSM& sparseProfiles);
void consistencyXZZY(PROFILE_TYPE weight, SparseMatrix* matXZ, SparseMatrix* matZY, VPROFILE& profile);
void consistencyZXZY(PROFILE_TYPE weight, SparseMatrix* matZX, SparseMatrix* matZY, VPROFILE& profile);
// --- DEBUG ---
void showProfile(VPROFILE* profile, int numRow, int numCol);
