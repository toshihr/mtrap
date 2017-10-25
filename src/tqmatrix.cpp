#include <vector>
#include <string>
#include <fstream>
#include <cfloat>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdint.h>
#include "scorematrix.h"
#include "tqmatrix.h"

#ifdef COMPRESS
  #include "snappy.h"
#endif

using namespace std;

namespace ScoreMatrix {

const PATH_TYPE PATH_INF = INT_MAX;
const PATH_TYPE PATH_ZERO = 0U;

const double FIXED_SCALE = 1000.0;

bool countModeSkipGapGapSite;
bool countModeSetZeroGapSite;
bool countModeDirectProb;
bool countModeWithoutLowerCase;
bool countModeSymmetrizeFreq;
bool countModeGlobalTQ;
bool countModeNoRound;
int countModeCluster;

double tmatrix[RANK+1][RANK+1][RANK+1][RANK+1];
PATH_TYPE distTQplusDiff_INT[RANK+2][RANK+2][RANK+1][RANK+1];
double distTQplusDiff[RANK+2][RANK+2][RANK+1][RANK+1];
PART_TYPE distTQplusDiffExt[RANK+2][RANK+2][RANK+1][RANK+1];
PART_TYPE distTQplusDiffExtBackward[RANK+2][RANK+2][RANK+1][RANK+1];
double epsilon;
bool gapWithoutTscore;

}

void ScoreMatrix::initTQ()
{
    for(int i1=0; i1<RANK+1; i1++)
    for(int i2=0; i2<RANK+1; i2++)
    for(int i3=0; i3<RANK+1; i3++)
    for(int i4=0; i4<RANK+1; i4++)
        tmatrix[i1][i2][i3][i4] = 0.0;
}

void ScoreMatrix::calcDistTQplusDiffTable(const GlobalOption& g) {
    const double e = ScoreMatrix::epsilon;
    const ScoreMatrix::EncodedResidue GAP = ScoreMatrix::encode(ScoreMatrix::GAP);
    const ScoreMatrix::EncodedResidue PRE = ScoreMatrix::encode('@');

    for(ScoreMatrix::EncodedResidue Ar_1=0; Ar_1<RANK+2; Ar_1++) // residues + gap + pre
    {
        for(ScoreMatrix::EncodedResidue Bc_1=0; Bc_1<RANK+2; Bc_1++) // residues + gap + pre
        {
            for(ScoreMatrix::EncodedResidue Ar=0; Ar<RANK+1; Ar++) // residues + gap
            {
                for(ScoreMatrix::EncodedResidue Bc=0; Bc<RANK+1; Bc++) // residues + gap
                {
                    assert(Ar != PRE && Bc != PRE);

                    const double TQ = (Ar_1 != PRE && Bc_1 != PRE) ? ScoreMatrix::getTQ(Ar_1, Bc_1, Ar, Bc) : 0.0;
                    const double nextTQ = (Ar_1 != PRE && Bc_1 != PRE) ? ScoreMatrix::getTQ(Ar, Bc, Ar_1, Bc_1) : 0.0;
                    double diffArBc, scoreArBc_part;

                    if(Ar != GAP && Bc != GAP)
                    {
                        // difference between amino acid Ar and Ac
                        diffArBc = ScoreMatrix::amatrix(Ar, Bc);
                        scoreArBc_part = ScoreMatrix::pmatrix(Ar, Bc);
                    }else if((Ar != GAP && Bc == GAP && Bc_1 != GAP) || (Ar == GAP && Bc != GAP && Ar_1 != GAP)){ // (Ar_1,Bc_1) -> (Ar,-) or (Ar_1,Bc_1) -> (-,Bc)
                        // gap open (inside & wallside)
                        diffArBc = g.gap_open_metric;
                        scoreArBc_part = g.gap_open_part;
                    }else if((Ar != GAP && Bc == GAP && Bc_1 == GAP) || (Ar == GAP && Bc != GAP && Ar_1 == GAP)){  // (Ar_1,-) -> (Ar,-) or (-,Bc_1) -> (-,Bc)
                        // gap ext (inside & wallside)
                        diffArBc = g.gap_ext_metric;
                        scoreArBc_part = g.gap_ext_part;
                    }else if(Ar == GAP && Bc == GAP){ // (Ar_1,Bc_1) -> (-,-)
                        // this part may not be used
                        diffArBc = 0.0;
                        scoreArBc_part = 0.0;
                    }else{
                        // error
                        cerr << ScoreMatrix::RESIDUE[Ar_1] << "," << ScoreMatrix::RESIDUE[Bc_1] << "," << ScoreMatrix::RESIDUE[Ar] << "," << ScoreMatrix::RESIDUE[Bc] << endl;
                        assert(false && "calcDistTQplusDiffTable() error");
                    }

                    // --- calc (1 - e) * Diff(Ar,Bc) + e * TQ(Ar-1,Bc-1,Ar,Bc) ---
                    // Table[Ar-1,Bc-1,Ar,Bc] = (1 - e) * TQ(Ar-1,Bc-1,Ar,Bc) + e * S(Ar,Bc)
                    distTQplusDiff[Ar_1][Bc_1][Ar][Bc] = (1.0 - e) * diffArBc + e * TQ;
                    distTQplusDiff_INT[Ar_1][Bc_1][Ar][Bc] = static_cast<PATH_TYPE>(distTQplusDiff[Ar_1][Bc_1][Ar][Bc] * ScoreMatrix::FIXED_SCALE);

                    // --- for partition function ---
#ifdef PART_EPS_ORIGINAL
                    // epsilonの値によって元々のProbalignのほうが敏感に値が変化してしまい、扱いにくい
                    // Table[Ar-1,Bc-1,Ar,Bc] = exp[ -beta * (1 - e) * TQ(Ar-1,Bc-1,Ar,Bc) ] * exp[ beta * e * S(Ar,Bc) ]
                    distTQplusDiffExt[Ar_1][Bc_1][Ar][Bc] = exp(g.beta * scoreArBc_part * (1.0 - e));
                    distTQplusDiffExt[Ar_1][Bc_1][Ar][Bc] *= exp(-g.betaTQ * TQ * e);

                    // Table[Ar+1,Bc+1,Ar,Bc] = exp[ -beta * (1 - e) * TQ(Ar,Bc,Ar+1,Bc+1) ] * exp[ beta * e * S(Ar,Bc) ]
                    distTQplusDiffExtBackward[Ar_1][Bc_1][Ar][Bc] = exp(g.beta * scoreArBc_part * (1.0 - e));
                    distTQplusDiffExtBackward[Ar_1][Bc_1][Ar][Bc] *= exp(-g.betaTQ * nextTQ * e);
#else
                    // partition functionではepsilonは使わず, TQの重みはbetaTQで代用
                    // Table[Ar-1,Bc-1,Ar,Bc] = exp[ -beta * (1 - e) * TQ(Ar-1,Bc-1,Ar,Bc) ] * exp[ beta * e * S(Ar,Bc) ]
                    distTQplusDiffExt[Ar_1][Bc_1][Ar][Bc] = exp(g.beta * scoreArBc_part);
                    distTQplusDiffExt[Ar_1][Bc_1][Ar][Bc] *= exp(-g.betaTQ * TQ);

                    // Table[Ar+1,Bc+1,Ar,Bc] = exp[ -beta * (1 - e) * TQ(Ar,Bc,Ar+1,Bc+1) ] * exp[ beta * e * S(Ar,Bc) ]
                    distTQplusDiffExtBackward[Ar_1][Bc_1][Ar][Bc] = exp(g.beta * scoreArBc_part);
                    distTQplusDiffExtBackward[Ar_1][Bc_1][Ar][Bc] *= exp(-g.betaTQ * nextTQ);
#endif

/*#ifdef DEBUG
                    cerr << ScoreMatrix::RESIDUE[Ar_1];
                    cerr << "," << ScoreMatrix::RESIDUE[Bc_1];
                    cerr << "->" << ScoreMatrix::RESIDUE[Ar];
                    cerr << "," << ScoreMatrix::RESIDUE[Bc];
                    cerr << ":" << distTQplusDiffBeta[Ar_1][Bc_1][Ar][Bc] << endl;

#endif*/
                }
            }
        }
    }


#ifdef DEBUG
    putLog("ScoreMatrix(tqmatrix.cpp)::calcDistTQplusDiffTable():");
    putLog("exp(-beta * score metric) matrix (eps=" + DS(e) + "):");

    stringstream header;
    header << "  ";
    for(int c=0; c<RANK; c++)
    {
        header << left << setw(5) << ScoreMatrix::RESIDUE[c];
    }
    putLog(header.str());

    for(int r=0; r<RANK; r++)
    {
        stringstream line;
        line << ScoreMatrix::RESIDUE[r] << ' ';
        for(int c=0; c<RANK; c++)
        {
            line << left << setw(5) << setprecision(2) << distTQplusDiffExt[PRE][PRE][r][c];
        }
        putLog(line.str());
    }
#endif

}


void ScoreMatrix::countToProb(const VVVVD& countTable, VVVVD& probTable){
    const int rmax = ScoreMatrix::RANK+1;

    VVD allCountForZTable;
    init2Vec(allCountForZTable, rmax);
    init4Vec(probTable, rmax);

    for(int i1=0; i1<rmax; ++i1){
        for(int i2=0; i2<rmax; ++i2){
            double allCountForZ = 0;
            for(int i3=0; i3<rmax; ++i3){
                for(int i4=0; i4<rmax; ++i4){
                    allCountForZ += countTable[i1][i2][i3][i4];
                }
            }
            allCountForZTable[i1][i2] = allCountForZ;
        }
    }

    for(int i1=0; i1<rmax; ++i1){
        for(int i2=0; i2<rmax; ++i2){
            for(int i3=0; i3<rmax; ++i3){
                for(int i4=0; i4<rmax; ++i4){
                    if(allCountForZTable[i1][i2] > 0){
                        probTable[i1][i2][i3][i4] = countTable[i1][i2][i3][i4]/allCountForZTable[i1][i2];
                    }else{
                        probTable[i1][i2][i3][i4] = 0.0;
                    }
                }
            }
        }
    }

 }


void ScoreMatrix::calcTQuantity(const VVVVD& probTable, VVVVD& quantityTable){
    const int rmax = ScoreMatrix::RANK+1;
    init4Vec(quantityTable, rmax);

    if(countModeDirectProb)
    {
        VVD maxProb;
        init2Vec(maxProb,rmax);

        for(int i1=0;i1<rmax; ++i1){
            for(int i2=0;i2<rmax; ++i2){
                for(int i3=0;i3<rmax; ++i3){
                    for(int i4=0;i4<rmax; ++i4){
                        const double prob = probTable[i1][i2][i3][i4];
                        if(prob > maxProb[i1][i2])
                        {
                            maxProb[i1][i2] = prob;
                        }
                    }
                }
            }
        }

        for(int i1=0;i1<rmax; ++i1){
            for(int i2=0;i2<rmax; ++i2){
                for(int i3=0;i3<rmax; ++i3){
                    for(int i4=0;i4<rmax; ++i4){
                        if(abs(maxProb[i1][i2]) > EPS)
                        {
                            quantityTable[i1][i2][i3][i4] = 1.0 - probTable[i1][i2][i3][i4] / maxProb[i1][i2];
                        }else{
                            quantityTable[i1][i2][i3][i4] = 1.0;
                        }
                    }
                }
            }
        }
    }else if(countModeGlobalTQ){
        double maxLog = 0.0;

        for(int i1=0;i1<rmax; ++i1){
            for(int i2=0;i2<rmax; ++i2){
                for(int i3=0;i3<rmax; ++i3){
                    for(int i4=0;i4<rmax; ++i4){
                        double logexy = ((double)(-1))*log(probTable[i1][i2][i3][i4]);
                        if(!(Isinf(logexy))){
                            if(logexy > maxLog)maxLog = logexy;
                        }
                    }
                }
            }
        }

        for(int i1=0; i1<rmax; ++i1){
            for(int i2=0; i2<rmax; ++i2){
                for(int i3=0; i3<rmax; ++i3){
                    for(int i4=0; i4<rmax; ++i4){
                        if(probTable[i1][i2][i3][i4] > 1.0e-15){
                            double logexy = ((double)(-1))*log(probTable[i1][i2][i3][i4]);
                            if(!(Isinf(logexy))){
                                quantityTable[i1][i2][i3][i4] = logexy/maxLog;
                            }else{
                                quantityTable[i1][i2][i3][i4] = 1.0;
                            }
                        }else{
                            quantityTable[i1][i2][i3][i4] = 1.0;
                        }
                    }
                }
            }
        }
    }else{
        // --- original transition quantity ---
        VVD maxLogExy;
        init2Vec(maxLogExy,rmax);

        for(int i1=0;i1<rmax; ++i1){
            for(int i2=0;i2<rmax; ++i2){
                for(int i3=0;i3<rmax; ++i3){
                    for(int i4=0;i4<rmax; ++i4){
                        double logexy = ((double)(-1))*log(probTable[i1][i2][i3][i4]);
                        if(!(Isinf(logexy))){
                            if(logexy > maxLogExy[i1][i2]){
                                maxLogExy[i1][i2] = logexy;
                                //cout << maxLogExy[i1][i2] << endl;
                            }
                        }
                    }
                }
            }
        }

        for(int i1=0; i1<rmax; ++i1){
            for(int i2=0; i2<rmax; ++i2){
                for(int i3=0; i3<rmax; ++i3){
                    for(int i4=0; i4<rmax; ++i4){
                        if(probTable[i1][i2][i3][i4] > 1.0e-15){
                            double logexy = ((double)(-1))*log(probTable[i1][i2][i3][i4]);
                            if(!(Isinf(logexy))){
                                quantityTable[i1][i2][i3][i4] = logexy/maxLogExy[i1][i2];
                            }else{
                                quantityTable[i1][i2][i3][i4] = 1.0;
                            }
                        }else{
                            quantityTable[i1][i2][i3][i4] = 1.0;
                        }
                    }
                }
            }
        }
    } // end of if countModeDirectProb
}

void ScoreMatrix::symmetrize(VVD& srcTable){
    const int rmax = static_cast<int>(srcTable.size());
    VVD destTable;
    init2Vec(destTable, rmax);

    for(int i1=0; i1<rmax; ++i1){
        for(int i2=0; i2<rmax; ++i2){
            destTable[i1][i2] = (srcTable[i1][i2]+srcTable[i2][i1])/2;
        }
    }

    for(int i1=0; i1<rmax; ++i1){
        for(int i2=0; i2<rmax; ++i2){
            srcTable[i1][i2] = destTable[i1][i2];
        }
    }
}

void ScoreMatrix::symmetrize(VVVVD& srcTable){
    const int rmax = static_cast<int>(srcTable.size());
    VVVVD destTable;
    init4Vec(destTable, rmax);

    for(int i1=0; i1<rmax; ++i1){
        for(int i2=0; i2<rmax; ++i2){
            for(int i3=0; i3<rmax; ++i3){
                for(int i4=0; i4<rmax; ++i4){
                    const double TQ = srcTable[i1][i2][i3][i4];
                    // AA -> BC and AA -> CB
                    if(i1 == i2){
                        const double rTQ = srcTable[i1][i2][i4][i3];
                        destTable[i1][i2][i3][i4] = (TQ+rTQ)/2;
                        destTable[i1][i2][i4][i3] = (TQ+rTQ)/2;
                    // AB -> CC and BA -> CC
                    }else if(i3 == i4){
                        const double rTQ = srcTable[i2][i1][i3][i4];
                        destTable[i1][i2][i3][i4] = (TQ+rTQ)/2;
                        destTable[i2][i1][i3][i4] = (TQ+rTQ)/2;
                    // AB -> CD and BA -> DC
                    }else{
                        double rTQ = srcTable[i2][i1][i4][i3];;
                        destTable[i1][i2][i3][i4] = (TQ+rTQ)/2;
                        destTable[i2][i1][i4][i3] = (TQ+rTQ)/2;
                    }
                }
            }
        }
    }

    for(int i1=0; i1<rmax; ++i1){
        for(int i2=0; i2<rmax; ++i2){
            for(int i3=0; i3<rmax; ++i3){
                for(int i4=0; i4<rmax; ++i4){
                    srcTable[i1][i2][i3][i4] = destTable[i1][i2][i3][i4];
                }
            }
        }
    }
}

// return the next non gap site position, if not exist then return the length of the sequences
int ScoreMatrix::getNextNonGapSitePosition(const string& s1, const string& s2, int startPos)
{
    assert(s1.size() == s2.size() && "the length of s1 and s2 is wrong");
    const ScoreMatrix::EncodedResidue EGAP = ScoreMatrix::encode(ScoreMatrix::GAP);
    const int l = s1.size();

    int retPos = startPos;
    for(; retPos<l; ++retPos)
    { // Non gap site search
        const ScoreMatrix::EncodedResidue R1 = s1[retPos];
        const ScoreMatrix::EncodedResidue R2 = s2[retPos];
        // Mode 1: Cancel gap-gap site
        if(ScoreMatrix::countModeSkipGapGapSite && R1 == EGAP && R2 == EGAP)
        {
            continue;
        }
        break;
    }
    return retPos;
}


void ScoreMatrix::loadFrequecyFromBLOCKS(const std::string& fileName, VD& count1Table, VVD& count2Table, VVVVD& count4Table, uint64_t& totAas, uint64_t& totPairs, uint64_t& totTuples)
{
    const int rmax_pure = 20;
    const int rmax = ScoreMatrix::RANK+1;
    const ScoreMatrix::EncodedResidue EGAP = ScoreMatrix::encode(ScoreMatrix::GAP);
    init1Vec(count1Table,rmax_pure);
    init2Vec(count2Table,rmax_pure);
    init4Vec(count4Table,rmax);
    totAas = 0;
    totPairs = 0;
    totTuples = 0;

    //int totBlocks = 0;

    string ac;
    string line;
    ifstream ifs(fileName.c_str());
    while(ifs && getline(ifs, line))
    {
        // STEP1: search AC line
        while(ifs)
        {
            if(line.size() >= 5 && line.substr(0,5) == "AC   "){
                ac = line.substr(5,8);
                break;
            }
            getline(ifs, line);
        }
        // STEP2: search data line
        while(ifs)
        {
            if(line.size() >= 13 && line.substr(2,3) != "   " && line[12] == '(') break;
            getline(ifs, line);
        }
        // STEP3: read a block data (same algorithm as fill_block in blocum.c)
        VS block_seqs;
        VI block_cluster_id; // cluster # for seq
        VI block_ncluster; // #seqs in same cluster
        //double block_totDiag = 0.0;
        //double block_totOffd = 0.0;
        int now_cluster_id = 0;
        int now_ncluster = 0;
        while(ifs)
        {
            if(line.size() == 0)
            { // new cluster
                // set #seqs in cluster to seqs in previous cluster
                if(now_ncluster > 0)
                {
                    for(int i=0; i<static_cast<int>(block_seqs.size()); i++)
                    {
                        if(block_cluster_id[i] == now_cluster_id) block_ncluster[i] = now_ncluster;
                    }
                }
                now_cluster_id++;
                now_ncluster = 0;
            }else if(line.size() > 0){
                if(line.size() >= 2 && line.substr(0,2) == "//") break;

                int loc = line.find(')',0);
                while(line[loc] == ' ') loc++;

                // add to block
                block_seqs.push_back(line.substr(loc));
                block_cluster_id.push_back(now_cluster_id);
                // add to cluster
                now_ncluster++;

                // count frequencies
                string& lastSeq = block_seqs.back();
                ScoreMatrix::encode(lastSeq);
                for(string::iterator i=lastSeq.begin(); i!=lastSeq.end(); ++i)
                {
                    if(*i != EGAP)
                    {
                        count1Table[*i]++;
                        totAas++;
                    }
                }
            }
            getline(ifs, line);
        }
        // compute weights for the last cluster
        if(now_ncluster > 0)
        {
            for(int i=0; i<static_cast<int>(block_seqs.size()); i++)
            {
                if(block_cluster_id[i] == now_cluster_id) block_ncluster[i] = now_ncluster;
            }
        }
        //totSeqs += static_cast<int>(block_seqs.size());
        //totWidth += static_cast<int>(block_seqs.at(0).size());

        // STEP4: re-clustering (same algorithm as cluster_seqs in blocum.c)
        if(ScoreMatrix::countModeCluster >= 0)
        {

        }
    }

}


void ScoreMatrix::BLOCKS_cluster_seqs(VS& block_seqs, VI& block_cluster_id, VI& block_ncluster)
{
    const int nseq = static_cast<int>(block_seqs.size());
    const int width = static_cast<int>(block_seqs.at(0).size());
    const int npair = nseq*(nseq-1)/2;
    //const int threshold = static_cast<int>(ScoreMatrix::countModeCluster*width/100);
    //int clus;
    VI nclus;
    //int minclus, oldclus;
    // first: %id, second: cluster #
    vector<std::pair<int,int> > pairs(npair);

    // compute local %ID for all possible pairs of sequences
    for(int s1=0; s1<nseq-1; s1++)
    {
        for(int s2=s1+1; s2<nseq; s2++)
        {
            int px = ScoreMatrix::BLOCKS_INDEX(nseq, s1, s2);
            pairs[px].first = 0;
            pairs[px].second = -1;
            for(int i=0; i<width; i++)
            {
                if(block_seqs[s1][i] == block_seqs[s2][i]) pairs[px].first++;
            }
        }
    }



}

void ScoreMatrix::loadFrequecy(const string& listFile, const map<string,int>& seqFreqTable, VD& count1Table, VVD& count2Table, VVVVD& count4Table, uint64_t& totAas, uint64_t& totPairs, uint64_t& totTuples)
{
    const int rmax_pure = 20;
    const int rmax = ScoreMatrix::RANK+1;
    const char EGAP = ScoreMatrix::encode(ScoreMatrix::GAP);
    init1Vec(count1Table,rmax_pure);
    init2Vec(count2Table,rmax_pure);
    init4Vec(count4Table,rmax);
    totAas = 0;
    totPairs = 0;
    totTuples = 0;

    ifstream ifs(listFile.c_str());
    string fastaFileName;
    while(ifs && getline(ifs, fastaFileName))
    {
        VS seqs;
        VS names;
        getMultiFASTA(names,seqs,fastaFileName,false,countModeWithoutLowerCase);
        for(VS::iterator i=seqs.begin(); i!=seqs.end(); ++i) ScoreMatrix::encode(*i);

        // n: the number of sequences in the alignment
        const int n = static_cast<int>(seqs.size());
        // l: the length of the alignment
        const int l = static_cast<int>(seqs[0].size());

        //cout << "COUNT: " << fastaFileName;
        //cout << "[" << n << "," << l << "]" << endl;

        // --- major loop ---
        for(int i=0; i<n-1; ++i)
        {
            for(int j=i+1; j<n; ++j)
            {
                assert(seqs[i].size() == seqs[j].size() && "size is wrong");
                // --- filter: local %ID > 0 ---
                {
                    int pw_n = 0;
                    int pw_len = 0;
                    for(int k=0; k<l; k++)
                    {
                        if(seqs[i][k] != EGAP && seqs[j][k] != EGAP)
                        {
                            if(seqs[i][k] == seqs[j][k])pw_n++;
                            pw_len++;
                        }
                    }
                    const double id = (pw_len > 0) ? static_cast<double>(pw_n) / static_cast<double>(pw_len) : 0.0;
                    if(abs(id) < EPS)continue;
                }

                int k = ScoreMatrix::getNextNonGapSitePosition(seqs[i], seqs[j], 0);
                while(k<l)
                {
                    // --- sequence weight ---
                    double weight = 1.0;
                    {
                        map<string, int>::const_iterator freq_i = seqFreqTable.find(names[i]);
                        map<string, int>::const_iterator freq_j = seqFreqTable.find(names[j]);
                        if(freq_i != seqFreqTable.end() && freq_j != seqFreqTable.end())
                        {
                            weight = weight / sqrt(static_cast<double>(freq_i->second * freq_j->second));
                        }
                    }

                    int next_k = ScoreMatrix::getNextNonGapSitePosition(seqs[i], seqs[j], k+1);

                    // --- pair frequency ---
                    const ScoreMatrix::EncodedResidue pre_Ri = seqs[i][k];
                    const ScoreMatrix::EncodedResidue pre_Rj = seqs[j][k];
                    if(pre_Ri < rmax_pure && pre_Rj < rmax_pure)
                    { // ignore including gap sites
                        count1Table[pre_Ri]++; // TODO: weight
                        totAas++;
                        count1Table[pre_Rj]++;
                        totAas++;
                        count2Table[pre_Ri][pre_Rj] += weight;
                        totPairs++;
                    }

                    // --- pair transition frequency ---
                    if(next_k<l)
                    {
                        const ScoreMatrix::EncodedResidue Ri = seqs[i][next_k];
                        const ScoreMatrix::EncodedResidue Rj = seqs[j][next_k];
                        count4Table[pre_Ri][pre_Rj][Ri][Rj] += weight;
                        totTuples++;
                    }

                    k = next_k;
                }
            }
        }

    }
    ifs.close();
}


/* Counting consequtive pairs of residues for the transition-quantity
    listFile: input database which contains filenames
    seqFreqFile: frequency information for weighting method
    outFile: output filename

    TODO: treating a gap. why blosum62.bla containing -4 as a gap value
*/
void ScoreMatrix::generateMatrix(const string& listFile, const string& seqFreqFile, const string& outFile)
{
    const int rmax_pure = 20;
    //const int rmax_extend = ScoreMatrix::RANK;
    const int rmax = ScoreMatrix::RANK+1;
    const char EGAP = ScoreMatrix::encode(ScoreMatrix::GAP);

    // --- STEP1: load frequencies of each sequence in the database for the weighting ---
    // input:  seqFreqFile
    // otuput: seqFreqTable
    map<string,int> seqFreqTable;
    cout << "load sequence frequencies..." <<  endl;
    if(seqFreqFile != "")
    {
        string buff;
        ifstream ifsFreq(seqFreqFile.c_str());
        while(ifsFreq && getline(ifsFreq, buff))
        {
            VS vs;
            mySplit(buff, ", ", vs);
            if(vs.size() > 0)
            {
                seqFreqTable[vs[0]] = SB<int>(vs[1]);
                cout << "add:" << vs[0] << "->" << vs[1] << endl;
            }
        }
        ifsFreq.close();
    }

    // --- STEP2: calculate frequencies ---
    // input:  listFile, seqFreqTable
    // output: count1Table, count2Table, count4Table: frequencies of symbols, pairs(Symbol,Symbol) / tuples(Symbol,Symbol,Symbol,Symbol) in the database
    uint64_t totAas, totPairs, totTuples;
    VD count1Table; // count1Table[i] is equal to AaFreq[i] in blosum
    VVD count2Table;
    VVVVD count4Table;

    cout << "calculate frequencies..." << endl;

    ScoreMatrix::loadFrequecy(listFile, seqFreqTable, count1Table, count2Table, count4Table, totAas, totPairs, totTuples);

    // --- STEP3: post process ---
    // input:  count4Table
    // output: count4Table
    cout << "post process..." << endl;
    if(ScoreMatrix::countModeSetZeroGapSite)
    {   // omit the transition including gap
        for(int i1=0; i1<rmax; ++i1)
        for(int i2=0; i2<rmax; ++i2)
        for(int i3=0; i3<rmax; ++i3)
        for(int i4=0; i4<rmax; ++i4)
        {
            if(i1 == EGAP || i2 == EGAP || i3 == EGAP || i4 == EGAP)
            {
                count4Table[i1][i2][i3][i4] = 0;
            }
        }
    }

    // --- STEP4 symmetrize ---
    ScoreMatrix::symmetrize(count2Table);
    if(countModeSymmetrizeFreq)
    {
        cout << "symmetrize frequencies..." << endl;
        ScoreMatrix::symmetrize(count4Table);
    }

    // --- STEP5: probability ---
    cout << "calculate probability..." << endl;
    VD prob1Table; // prob1Table[i] = Pi = AaPairs[i]/totPairs (not equal to AaFreq[i]) in blosum
    VVD prob2Table; // equal to qij in blosum
    init1Vec(prob1Table, rmax_pure);
    init2Vec(prob2Table, rmax_pure);
    double totPairProbs = 0.0;
    for(int i1=0; i1<rmax_pure; ++i1)
    {
        for(int i2=0; i2<rmax_pure; ++i2)
        {
            prob2Table[i1][i2] = count2Table[i1][i2]/totPairs;
            totPairProbs += prob2Table[i1][i2];
            prob1Table[i1] += prob2Table[i1][i2];
        }
        cout << "P(" << ScoreMatrix::decode(i1) << ")=" << prob1Table[i1] << endl;
    }
    cout << "totPairProbs=" << totPairProbs << endl;

    VVVVD prob4Table;
    ScoreMatrix::countToProb(count4Table,prob4Table);

    // --- STEP6 Substitution Matrix ---
    cout << "calculate substitution matrix..." << endl;
    double mutualEntropy = 0.0;
    VVD oddsTable; // equal to oij in blosum
    VVD logOddsTable; // equal to sij in blosum
    ScoreMatrix::AminoMatrix substitutionMatrix; // equal to iijMatrix in blosum
    init2Vec(oddsTable, rmax_pure);
    init2Vec(logOddsTable, rmax_pure);
    substitutionMatrix.setZero();
    for(int i1=0; i1<rmax_pure; ++i1)
    {
        for(int i2=0; i2<rmax_pure; ++i2)
        {
            if(prob1Table[i1] * prob1Table[i2] > EPS)
            {
                oddsTable[i1][i2] = prob2Table[i1][i2] / prob1Table[i1] / prob1Table[i2];
            }else{
                oddsTable[i1][i2] = 0;
            }
            //cout << i1 << "," << i2 << "oij:" << oddsTable[i1][i2] << endl;
            if(oddsTable[i1][i2] > EPS)
            {
                logOddsTable[i1][i2] = Log2(oddsTable[i1][i2]);
            }else{
                logOddsTable[i1][i2] = -20.0; // -infinity
            }
            //cout << i1 << "," << i2 << "sij:" << logOddsTable[i1][i2] << endl;
            mutualEntropy += prob2Table[i1][i2] * Log2(oddsTable[i1][i2]);
        }
    }
    const double iscale = max(2.0,Round(2.0/sqrt(mutualEntropy)));
    for(int i1=0; i1<rmax_pure; ++i1)
    {
        for(int i2=0; i2<rmax_pure; ++i2)
        {
            substitutionMatrix(i1,i2) = logOddsTable[i1][i2] * iscale;
        }
    }

    // --- treatment of BZX (blosum style) ---
    cout << "treat BZX..." << endl;
    VD& FF = count1Table;
    VVD& SS = logOddsTable;
    const ScoreMatrix::EncodedResidue B = ScoreMatrix::encode('B');
    const ScoreMatrix::EncodedResidue N = ScoreMatrix::encode('N');
    const ScoreMatrix::EncodedResidue D = ScoreMatrix::encode('D');
    const ScoreMatrix::EncodedResidue Z = ScoreMatrix::encode('Z');
    const ScoreMatrix::EncodedResidue Q = ScoreMatrix::encode('Q');
    const ScoreMatrix::EncodedResidue E = ScoreMatrix::encode('E');
    const ScoreMatrix::EncodedResidue X = ScoreMatrix::encode('X');
    { // (B,i) = (freq(N)*(N,i) + freq(D)*(D,i))*iscale / (freq(N)+freq(D))
        for(int i1=0; i1<rmax_pure; ++i1)
        {
            substitutionMatrix(i1,B) = (FF[N] * SS[i1][N] + FF[D] * SS[i1][D]) / (FF[N]+FF[D]) * iscale;
            substitutionMatrix(B,i1) = substitutionMatrix(i1,B);
        }
        // (B,B)
        substitutionMatrix(B,B) = (FF[N]*FF[N]*SS[N][N]*iscale
                                        + FF[N]*FF[D]*SS[N][D]*iscale
                                        + FF[D]*FF[N]*SS[D][N]*iscale
                                        + FF[N]*FF[N]*SS[D][D]*iscale)
                / ((FF[N]+FF[D])*(FF[N]+FF[D]));
    }
    { // (Z,i) = (freq(Q)*(Q,i) + freq(E)*(E,i))*iscale / (freq(Q)+freq(E))
        for(int i1=0; i1<rmax_pure; ++i1)
        {
            substitutionMatrix(i1,Z) = (FF[Q] * SS[i1][Q] + FF[E] * SS[i1][E]) / (FF[Q]+FF[E]) * iscale;
            substitutionMatrix(Z,i1) = substitutionMatrix(i1,Z);
        }
        // (Z,B)
        const double dtemp = (FF[N]*FF[Q]*SS[N][Q]*iscale
                            + FF[N]*FF[E]*SS[N][E]*iscale
                            + FF[D]*FF[Q]*SS[D][Q]*iscale
                            + FF[D]*FF[E]*SS[D][E]*iscale) / ((FF[N]+FF[D])*(FF[N]+FF[D]));
        const double dtemp1= (FF[Q]*FF[N]*SS[Q][N]*iscale
                            + FF[Q]*FF[D]*SS[Q][D]*iscale
                            + FF[E]*FF[N]*SS[E][N]*iscale
                            + FF[E]*FF[D]*SS[E][D]*iscale) / ((FF[Q]+FF[E])*(FF[Q]+FF[E]));
        substitutionMatrix(Z,B) = dtemp*dtemp1;
        substitutionMatrix(B,Z) = substitutionMatrix(Z,B);
        // (Z,Z)
        substitutionMatrix(Z,Z) = (FF[Q]*FF[Q]*SS[Q][Q]*iscale
                                  + FF[Q]*FF[E]*SS[Q][E]*iscale
                                  + FF[E]*FF[Q]*SS[E][Q]*iscale
                                  + FF[E]*FF[E]*SS[E][E]*iscale)
                / ((FF[Q]+FF[E])*(FF[Q]+FF[E]));
    }
    { // X
        { // (X,X)
            double xx = 0.0, dtemp = 0.0;
            for(int i1=0; i1<rmax_pure; ++i1)
            {
                double x = 0.0;
                for(int i2=0; i2<rmax_pure; ++i2)
                {
                    x += FF[i2]*SS[i1][i2]*iscale;
                    xx += FF[i1]*FF[i2]*SS[i1][i2]*iscale;
                    dtemp += FF[i1]*FF[i2];
                }
                x /= static_cast<double>(totAas);
                substitutionMatrix(i1,X) = x;
                substitutionMatrix(X,i1) = substitutionMatrix(i1,X);
            }
            substitutionMatrix(X,X) = xx/dtemp;
        }
        { // (X,i)
            double x = 0.0, xx = 0.0;
            for(int i1=0; i1<rmax_pure; ++i1)
            {
                x += FF[i1]*FF[N]*SS[i1][N]*iscale
                    +FF[i1]*FF[D]*SS[i1][D]*iscale;
                xx += FF[i1]*FF[Q]*SS[i1][Q]*iscale
                     +FF[i1]*FF[E]*SS[i1][E]*iscale;
            }
            x = x / (FF[N]+FF[D]) / static_cast<double>(totAas);
            xx = xx / (FF[Q]+FF[E]) / static_cast<double>(totAas);
            substitutionMatrix(X,B) = x;
            substitutionMatrix(B,X) = substitutionMatrix(X,B);
            substitutionMatrix(X,Z) = xx;
            substitutionMatrix(Z,X) = substitutionMatrix(X,Z);
        }
    }

    // --- STEP7 Transition Quantity ---
    cout << "calculate transition quantity" << endl;
    VVVVD quantityTable;
    ScoreMatrix::calcTQuantity(prob4Table,quantityTable);
    if(countModeSymmetrizeFreq == false)ScoreMatrix::symmetrize(quantityTable);

    // --- STEP8: generate text file ---
    cout << "generate outputs..." << endl;
    ScoreMatrix::saveEMBOSSformat(outFile + ".mat", substitutionMatrix, iscale, listFile, mutualEntropy);

    { // transition quantity
        ofstream ofs((outFile + ".tq").c_str());
        ofs << "format_ITO=" << totTuples << endl;
        for(int i1=0; i1<rmax; ++i1)
        {
            const char c1 = ScoreMatrix::decode(i1);
            for(int i2=0; i2<rmax; ++i2)
            {
                const char c2 = ScoreMatrix::decode(i2);
                for(int i3=0; i3<rmax; ++i3)
                {
                    const char c3 = ScoreMatrix::decode(i3);
                    for(int i4=0; i4<rmax; ++i4)
                    {
                        const char c4 = ScoreMatrix::decode(i4);
                        // AAAA,TQ,prob,symmetrized freq
                        ofs << c1 << c2 << c3 << c4;
                        ofs << "," << quantityTable[i1][i2][i3][i4];
                        ofs << "," << prob4Table[i1][i2][i3][i4];
                        ofs << "," << count4Table[i1][i2][i3][i4];
                        ofs << endl;
                    }
                }
            }
        }
        ofs.close();
    }
    // --- STEP9: generate binary file ---
    ScoreMatrix::loadTMatrix((outFile + ".tq").c_str(), false);
    ScoreMatrix::outputBinaryTQ((outFile + ".btq").c_str());
}


void ScoreMatrix::outputBinaryTQ(const string& outFile)
{
    const int rmax = ScoreMatrix::RANK+1;
    const int size_tmatrix = sizeof(double)*rmax*rmax*rmax*rmax;
#ifdef COMPRESS
    string buff;
    const uint32_t output_size = static_cast<uint32_t>(snappy::Compress((char *)&ScoreMatrix::tmatrix,size_tmatrix,&buff));

    ofstream ofs(outFile.c_str(), ios::out | ios::trunc | ios::binary);
    ofs << "format_ZIP=" << output_size << '@';
    ofs << 'b'; // flag for the check
    ofs.write(&*buff.begin(), output_size );
    ofs << 'e'; // flag for the check
    ofs.close();
#else
    string buff;
    ofstream ofs(outFile.c_str(), ios::out | ios::trunc | ios::binary);
    ofs << "format_BIN2=" << size_tmatrix << '@';
    ofs << 'b'; // flag for the check
    ofs.write(( char * ) &ScoreMatrix::tmatrix, size_tmatrix );
    ofs << 'e'; // flag for the check
    ofs.close();
#endif
}

bool ScoreMatrix::loadTMatrix(const string& fileName, bool ignorePostprocess)
{
    ifstream ifs(fileName.c_str(), ios::binary);
    if(!ifs) putError("Can't open the transition-quantity matrix.");

    string buff;
    int iAA[4];

    initTQ();

    // read header
    getline(ifs,buff,'=');
    if(buff == "format_ITO"){
        // omit the right side of =
        getline(ifs,buff);
        // read data
        while(ifs && getline(ifs, buff)){
            VS vs;
            mySplit(buff, ", ", vs);
            VS::iterator it = vs.begin();
            iAA[0] = RESIDUE_INDEX[(*it).at(0)];
            iAA[1] = RESIDUE_INDEX[(*it).at(1)];
            iAA[2] = RESIDUE_INDEX[(*it).at(2)];
            iAA[3] = RESIDUE_INDEX[(*it).at(3)];
            it++;
            const double p = string2binary<double>(*it);
            tmatrix[iAA[0]][iAA[1]][iAA[2]][iAA[3]] = p;
        }
    }else if(buff == "format_ZIP"){
#ifndef COMPRESS
        putError("this binary does not support compressed file. should be compiled with snappy library.", true);
#endif
        getline(ifs,buff,'@');
        const uint32_t size_compressed_tmatrix = SB<uint32_t>(buff);
#ifdef DEBUG
        cout << "size of compressed tmatrix = " << size_compressed_tmatrix << endl;
#endif
        char c;
        ifs >> c;										// 確認フラグ
        assert(c == 'b');
        buff.resize(size_compressed_tmatrix);
        ifs.read((char *)&*buff.begin(),size_compressed_tmatrix);
        ifs >> c;										// 確認フラグ
        assert(c == 'e' && "THE_SIZE_OF_TMATRIX_IS_WRONG");
        snappy::RawUncompress((char *)&*buff.begin(),size_compressed_tmatrix,(char *)&ScoreMatrix::tmatrix);
    }else if(buff == "format_BIN2"){
        getline(ifs,buff,'@');
        const unsigned int size_tmatrix = static_cast<int>(sizeof(double) * pow(static_cast<double>(ScoreMatrix::RANK+1),4));
#ifdef DEBUG
        cout << "size of tmatrix = " << size_tmatrix << endl;
#endif
        assert(SB<uint32_t>(buff) == size_tmatrix && "size is wrong!!");

        char c;
        ifs >> c;										// 確認フラグ
        assert(c == 'b');
        ifs.read(( char * ) &ScoreMatrix::tmatrix, size_tmatrix );
        ifs >> c;										// 確認フラグ
        assert(c == 'e' && "THE_SIZE_OF_TMATRIX_IS_WRONG");
    }else if(buff == "format_BIN"){
        getline(ifs,buff);
        const unsigned int size_tmatrix = static_cast<int>(sizeof(double) * pow(static_cast<double>(ScoreMatrix::RANK+1),4));
#ifdef DEBUG
        cout << "size of tmatrix = " << size_tmatrix << endl;
#endif
        assert(SB<uint32_t>(buff) == size_tmatrix && "size is wrong!!");

        char c;
        ifs >> c;										// 確認フラグ
        assert(c == 'b');
        ifs.read(( char * ) &ScoreMatrix::tmatrix, size_tmatrix );
        ifs >> c;										// 確認フラグ
        assert(c == 'e' && "THE_SIZE_OF_TMATRIX_IS_WRONG");
    }else{
    //旧フォーマットを読み込むとき
        ifs.close();
        ifs.open(fileName.c_str());
        while(ifs && getline(ifs, buff))
        {
//			tokenizer tokens(buff, sep);
//			tokenizer::iterator it=tokens.begin();
            vector<string> vs;
            mySplit(buff, ", ", vs);
            vector<string>::iterator it = vs.begin();
            //cout << *it << endl;
            iAA[0] = RESIDUE_INDEX[(*it).at(0)];
            iAA[1] = RESIDUE_INDEX[(*it).at(1)];
            it++;
            iAA[2] = RESIDUE_INDEX[(*it).at(0)];
            iAA[3] = RESIDUE_INDEX[(*it).at(1)];
            it++;
            const double p = string2binary<double>(*it);
            tmatrix[iAA[0]][iAA[1]][iAA[2]][iAA[3]] = p;
        }
    }

    ifs.close();

    // 拡張アミノ酸表記B=D or N, Z=E or Q, Xへの対応．ギャップとの距離以外を平均化
    if(ignorePostprocess == false && averageExtendedAA == true)
    {
        // B,Z,Xの平均化の順番で結果が若干変わるが気にしない
        // === B ===
        // [Bj] = [Dj] + [Nj]
        // [ik]   [ik] + [ik]
        for(int i=0; i<RANK-1; i++)
            if( i != RESIDUE_INDEX['B'] ) for(int j=0; j<RANK-1; j++)
                for(int k=0; k<RANK-1; k++)
                {
                    tmatrix[RESIDUE_INDEX['B']][i][j][k] = (tmatrix[RESIDUE_INDEX['D']][i][j][k] + tmatrix[RESIDUE_INDEX['N']][i][j][k]) / 2;
                }
        // 対象行列であるこを仮定．同時にやると微妙にへんなことになる可能性がありそうなのでループを分けて行う．
        // [ik]   [Bj]
        // [Bj] = [ik]
        for(int i=0; i<RANK-1; i++)
            if( i != RESIDUE_INDEX['B'] ) for(int j=0; j<RANK-1; j++)
                for(int k=0; k<RANK-1; k++)
                {
                    tmatrix[i][RESIDUE_INDEX['B']][k][j] = tmatrix[RESIDUE_INDEX['B']][i][j][k];
                }
        // [jB] = [jD] + [jN]
        // [ki]   [ki] + [ki]
        for(int i=0; i<RANK-1; i++)
            if( i != RESIDUE_INDEX['B'] ) for(int j=0; j<RANK-1; j++)
                for(int k=0; k<RANK-1; k++)
                {
                    tmatrix[j][k][RESIDUE_INDEX['B']][i] = (tmatrix[j][k][RESIDUE_INDEX['D']][i] + tmatrix[j][k][RESIDUE_INDEX['N']][i]) / 2;
                }
        // 対象行列であるこを仮定．同時にやると微妙にへんなことになるのでループを分けて行う．
        // [ki]   [jB]
        // [jB] = [ki]
        for(int i=0; i<RANK-1; i++)
            if( i != RESIDUE_INDEX['B'] ) for(int j=0; j<RANK-1; j++)
                for(int k=0; k<RANK-1; k++)
                {
                    tmatrix[k][j][i][RESIDUE_INDEX['B']] = tmatrix[j][k][RESIDUE_INDEX['B']][i];
                }
        // === Z ===
        // [Zj] = [Ej] + [Qj]
        // [ik]   [ik] + [ik]
        for(int i=0; i<RANK-1; i++)
            if( i != RESIDUE_INDEX['Z'] ) for(int j=0; j<RANK-1; j++)
                if( j != RESIDUE_INDEX['Z'] ) for(int k=0; k<RANK-1; k++)
                {
                    if( k != RESIDUE_INDEX['Z'] ) tmatrix[RESIDUE_INDEX['Z']][i][j][k] = (tmatrix[RESIDUE_INDEX['E']][i][j][k] + tmatrix[RESIDUE_INDEX['Q']][i][j][k]) / 2;
                }
        // 対象行列であるこを仮定．同時にやると微妙にへんなことになる可能性がありそうなのでループを分けて行う．
        // [ik]   [Zj]
        // [Zj] = [ik]
        for(int i=0; i<RANK-1; i++)
            if( i != RESIDUE_INDEX['Z'] ) for(int j=0; j<RANK-1; j++)
                for(int k=0; k<RANK-1; k++)
                {
                    tmatrix[i][RESIDUE_INDEX['Z']][k][j] = tmatrix[RESIDUE_INDEX['Z']][i][j][k];
                }
        // [jZ] = [jE] + [jQ]
        // [ki]   [ki] + [ki]
        for(int i=0; i<RANK-1; i++)
            if( i != RESIDUE_INDEX['Z'] ) for(int j=0; j<RANK-1; j++)
                for(int k=0; k<RANK-1; k++)
                {
                    tmatrix[j][k][RESIDUE_INDEX['Z']][i] = (tmatrix[j][k][RESIDUE_INDEX['E']][i] + tmatrix[j][k][RESIDUE_INDEX['Q']][i]) / 2;
                }
        // 対象行列であるこを仮定．同時にやると微妙にへんなことになるのでループを分けて行う．
        // [ki]   [jZ]
        // [jZ] = [ki]
        for(int i=0; i<RANK-1; i++)
            if( i != RESIDUE_INDEX['Z'] ) for(int j=0; j<RANK-1; j++)
                for(int k=0; k<RANK-1; k++)
                {
                    tmatrix[k][j][i][RESIDUE_INDEX['Z']] = tmatrix[j][k][RESIDUE_INDEX['Z']][i];
                }
        // === X ===
        // [Xj] = Σ[lj] / (RANK-2)
        // [ik]    [ik]
        for(int i=0; i<RANK-1; i++)
            if( i != RESIDUE_INDEX['X'] ) for(int j=0; j<RANK-1; j++)
                for(int k=0; k<RANK-1; k++)
                {
                    double sum=0.0;
                    for(int l=0; l<RANK-1; l++) if(l != RESIDUE_INDEX['X']) sum += tmatrix[l][i][j][k];
                    tmatrix[RESIDUE_INDEX['X']][i][j][k] = sum / (RANK-2);
                }
        // 対象行列であるこを仮定．同時にやると微妙にへんなことになる可能性がありそうなのでループを分けて行う．
        // [ik]   [Xj]
        // [Xj] = [ik]
        for(int i=0; i<RANK-1; i++)
            if( i != RESIDUE_INDEX['X'] ) for(int j=0; j<RANK-1; j++)
                for(int k=0; k<RANK-1; k++)
                {
                    tmatrix[i][RESIDUE_INDEX['X']][k][j] = tmatrix[RESIDUE_INDEX['X']][i][j][k];
                }
        // [jX] = Σ[jl] / (RANK-2)
        // [ki]    [ki]
        for(int i=0; i<RANK-1; i++)
            if( i != RESIDUE_INDEX['X'] ) for(int j=0; j<RANK-1; j++)
                for(int k=0; k<RANK-1; k++)
                {
                    double sum=0.0;
                    for(int l=0; l<RANK-1; l++) if(l != RESIDUE_INDEX['X']) sum += tmatrix[j][k][l][i];
                    tmatrix[j][k][RESIDUE_INDEX['X']][i] = sum / (RANK-2);
                }
        // 対象行列であるこを仮定．同時にやると微妙にへんなことになる可能性がありそうなのでループを分けて行う．
        // [ki]   [jX]
        // [jX] = [ki]
        for(int i=0; i<RANK-1; i++)
            if( i != RESIDUE_INDEX['X'] ) for(int j=0; j<RANK-1; j++)
                for(int k=0; k<RANK-1; k++)
                {
                    tmatrix[k][j][i][RESIDUE_INDEX['X']] = tmatrix[j][k][RESIDUE_INDEX['X']][i];
                }

    } // end of	if(averageExtendedAA)

    return true;
}
