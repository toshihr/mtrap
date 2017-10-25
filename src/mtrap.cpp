/*
 @file
 @author Toshihide Hara, Tokyo university of science

 @brief MTRAP: (Pairwise) Sequence alignment algorithm by a new measure based on transition probability between two consecutive pairs of residues

 === HOW TO USE ===
 [1] Sequence Alignment
 mtrap [option] inName outName
 [2] Output Binary format Transition-quantity Matrix
 mtrap -tm inName -outbintm outName
 or
 bash> for i in `ls *.tq` ; do mtrap -tm $i -outbintm ${i%.*}.btq  ; done

=== TODO ===
20120316
= general =
・new が失敗したときNULLがかえってくるとはかぎらない
・staticは，初期化においてマルチスレッド+VCで問題を起こす
20120311
・OPEN MPをOnにすると，Profile作成は高速化するがIterativeRefinementsはたまに遅くなる．原因不明．HyperThreadやキャッシュが原因？
= genetic distance & guide tree =
・案内木構築法としてNJ,WSPS等を実装
・案内木のための遺伝距離行列としてEER以外にもp-distance(clustal互換)等を実装すべき

 $Revision: $
*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <cctype>
#include <map>
#include <string>
#include <fstream>
#include <cmath>
#include <ctime>
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <shand/format.hpp>
#include "distance.h"
#include "scorematrix.h"
#include "tqmatrix.h"
#include "mtrap.h"
#include "pairwise.h"
#include "probalign.h"
#include "utility.h"
#include "entropy.h"
#include "primarylibrary.h"
#include "search.h"
#include "sequences.h"
#include "globaloption.h"
#include "tree.h"
#include "estimatefamily.h"
#include "config.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

void help_line(const string& option, const string& value, const string& explain){
    cerr << left << setw(10) << option;
    if(value != "")
    {
        cerr << right << setw(17) << ("[" + value + "]  ");
    }else{
        cerr << right << setw(17) << " ";
    }
    cerr << explain;
    cerr << endl;
}

void printHelp(){
    cerr << "==================" << endl;
    cerr << "  MTRAP Ver. " << string(VERSION);
#ifdef _OPENMP
    cerr << " [USE OPENMP]";
#endif
#ifdef USE_PROBALIGN
    cerr << " [USE PROBALIGN]";
#endif
#ifdef DEBUG
    cerr << " [DEBUG VERSION]";
#endif
    cerr << endl;

    cerr << "=== HOW TO USE ===" << endl;
    cerr << "mtrap [OPTIONS] inputfile outputfile" << endl;
    cerr << endl;
    cerr << "===   OPTIONS  === " << endl;
    cerr << "= GENERAL =" << endl;
    help_line("-I","...","add search path e.g. -I ~/mymatrix");
    help_line("-filelist","...","");
    help_line("-i","...","input file (FASTA format)");
    help_line("-o","...","output file");
    help_line("-nosort","","the results are not sorted");
    help_line("-noestimation","","do not estimate a family");
    help_line("-cds","GENETICCODES (comma separate)","\"protein\" coding DNA sequence");
    cerr << "             GENETICCODE 0: The Standard Code (transl_table=1)" << endl;
    cerr << "                         1: The Vertebrate Mitochondrial Code (transl_table=2)" << endl;
    cerr << "                         2: The Invertebrate Mitochondrial Code (transl_table=5)" << endl;
    cerr << "                         3: The Bacterial, Archaeal and Plant Plastid Code (transl_table=11)" << endl;
    cerr << endl;
    cerr << "= GAP =" << endl;
    help_line("-go",DS(globalOption.gap_open),"gap open [score]");
    help_line("-ge",DS(globalOption.gap_ext),"gap extension [score]");
    cerr << endl;
    cerr << "= MULTIPLE ALIGNMENT =" << endl;
    help_line("-itr",IS(globalOption.numIterativeRefReps),"set the number of iterations");
    help_line("-c",IS(globalOption.numConsistencyReps),"set the number of consistency transformations");
    cerr << endl;
    cerr << "= SUBSTITUTION MATRIX =" << endl;
    help_line("-m",globalOption.smatrixName,"substitution matrix e.g. CGONNET250, EBLOSUM62");
    help_line("-rm",globalOption.rmatrixName,"ramachandran matrix e.g. MRAMA1");
    help_line("-averageExtendedAA","","average extended amino acids (B,Z and X) substitutions");
    cerr << endl;
    cerr << "= TRANSITION QUANTITY MATRIX =" << endl;
    help_line("-tm",globalOption.tmatrixName,"transition-quantity matrix");
    help_line("-e",DS(ScoreMatrix::epsilon),"epsilon, weight for transition quantity");
    cerr << endl;
    cerr << "= RAMACHANDRAN QUANTITY =" << endl;
    help_line("-gamma",DS(ScoreMatrix::gamma),"gamma, weight for ramachandran matrix");
    cerr << endl;
    cerr << "= PARTITION FUNCTION POSTERIOR PROBABILITY =" << endl;
    help_line("-pf",DS(1 - globalOption.mtrap_degree),"the degree of partition function");
    help_line("-pm",globalOption.pmatrixName,"substitution matrix used for partition function");
    help_line("-pgo",DS(globalOption.gap_open_part),"gap open [score] for partition function");
    help_line("-pge",DS(globalOption.gap_ext_part),"gap extension [score] for partition function");
    help_line("-beta",DS(globalOption.beta),"beta, weight for partition function");
#ifndef USE_PROBALIGN
    help_line("-beta2",DS(globalOption.betaTQ),"beta2, weight for a part of TQ of partition function");
#endif
    cerr << endl;
    cerr << "= GENERATING MODE =" << endl;
    help_line("-primarylibrary","ID","output T-Coffee primary library, ID: sequence identity (same as T-Coffee)");
    cerr << endl;
    help_line("-count","MODES FILELIST WEIGHTLIST","count the frequency of transition");
    cerr << "           MODES 0: normal 1:skip gap-gap site 2:set zero gap site" << endl;
    cerr << "                 3: direct prob 4: without lower case 5: global TQ" << endl;
    cerr << "                 f: symmetrize the frequencies not the tq" << endl;
    cerr << "                 n: do not round the substitution matrix" << endl;
    help_line("-outbintm","FILENAME","output binary format transition-quantity matrix");
    cerr << endl;
//    cerr << "= DATABASE SEARCH MODE (UNDER CONSTRUCTION) =" << endl;
//    cerr << "-search query database : query and database must be multiple fasta format. The option -o should be used together." << endl;
//    cerr << "-search_allpair database : database must be ASTRAL SCOP fasta format. The option -o should be used together." << endl;
//    cerr << "-searchlimit query database" << endl;
//    cerr << endl;
//  cerr << "= DEBUG =" << endl;
//	cerr << "-pass output pass matrix filename" << endl;
//	cerr << "-oldeer use traditional entropy evolution rate" << endl;
//	cerr << "-gmode value[0] 0: A,B 1: width=1(B=0)  2: mirror only(A=-1,B=0)" << endl;
//	cerr << "!!! Set -gmode before -m -tm !!!" << endl;
//	cerr << "-withoutT: Gap cost without transition value" << endl;
}


//! load the list for batch style aligment.
/*!

    \param fileName comma separated style, i.e. input filename, output filename[ret].
*/
void loadFilelist(const string& fileName)
{
    ifstream ifs(fileName.c_str());
    string buff;

    while(ifs && getline(ifs, buff))
    {
        vector<string> vs;
        mySplit(buff, ", ", vs);
        if(vs.size() == 0)break;
        globalOption.inputFilenameList.push_back(vs[0]);
        globalOption.outputFilenameList.push_back(vs[1]);
    }

    ifs.close();
}

void setOptionDefaultValue()
{
    // set primary default value
    globalOption.tmatrixName = "SABmark1.63_sup_weighted.btq";
    globalOption.smatrixName = "VTML200I"; // "CGONNET250";
    globalOption.pmatrixName = "VTML200I"; // CGONNET160";
    globalOption.rmatrixName = "MRAMA1";
    globalOption.fmatrixDir = "fsmatrix-sabmark1.65";

    globalOption.csvSCOP_SAB_ID_Name = "scop1.75-id-with-sab1.65-group.csv";
    globalOption.csvSAB_seq_group_Name = "fsmatrix-sabmark1.65/summary_sequence.csv";

    globalOption.familyEstimation = false;
    globalOption.familyDatabaseName = "fsmatrix-sabmark1.65/SABmark1.65.fasta";

    globalOption.CDSMode = false;
    globalOption.CDS_type.clear();

    globalOption.gap_open = -11;
    globalOption.gap_ext = -0.3;
    globalOption.gap_open_part = -22;
    globalOption.gap_ext_part = -1;
    globalOption.mtrap_degree = 0.3;
    globalOption.beta = 0.2;
    globalOption.betaTQ = 1.5;
    ScoreMatrix::gapWithoutTscore = false;
    ScoreMatrix::epsilon = 0.775;
    ScoreMatrix::gamma = 0.0;
    globalOption.primaryLibraryDistance = Distance::GID;
    globalOption.numConsistencyReps = 2;
    globalOption.numIterativeRefReps = 10;
    globalOption.sortByInputOrder = true;

    // set secondary default value


}

int setOption(int argc, char ** argv, int priority)
{
    vec_string vecARGV;
    for(int i=1; i<argc; i++)vecARGV.push_back(string(argv[i]));

    // First priority options
    if(priority==0)
    {
        // default value
        setOptionDefaultValue();

        for(vector<string>::iterator i=vecARGV.begin(); i!=vecARGV.end(); ++i)
        {
            if(i->size() > 0 && i->at(0) == '-')
            {
                if(*i == "-I"){
                    i++; fileSearchDirList.push_back(*i + "/");
                }
            }
        }
        return 0;
    }else if(priority==1){
        // Second priority options
        for(vector<string>::iterator i=vecARGV.begin(); i!=vecARGV.end(); ++i)
        {
            if(i->at(0) == '-')
            {
                if(*i == "-o"){i++; globalOption.outputFilenameList.push_back(*i);
                }else if(*i == "-i"){i++; globalOption.inputFilenameList.push_back(*i);
                }else if(*i == "-nosort"){ globalOption.sortByInputOrder = false;
                }else if(*i == "-noestimation"){ globalOption.familyEstimation = false;
                }else if(*i == "-cds"){
                    i++;
                    globalOption.CDSMode = true;
                    VS tmpVec;
                    mySplit(*i, ",", tmpVec);
                    globalOption.CDS_type.clear();
                    for(VS::const_iterator itr=tmpVec.begin(); itr!=tmpVec.end(); ++itr) globalOption.CDS_type.push_back(SI(*itr));
                }else if(*i == "-pass"){i++; globalOption.outputFilename_pass = *i;
                }else if(*i == "-go"){i++; globalOption.gap_open = string2binary<double>(*i);
                }else if(*i == "-ge"){i++; globalOption.gap_ext = string2binary<double>(*i);
                }else if(*i == "-pgo"){i++; globalOption.gap_open_part = string2binary<double>(*i);
                }else if(*i == "-pge"){i++; globalOption.gap_ext_part = string2binary<double>(*i);
                }else if(*i == "-itr"){i++; globalOption.numIterativeRefReps = string2binary<int>(*i);
                }else if(*i == "-c"){i++; globalOption.numConsistencyReps = string2binary<int>(*i);
                }else if(*i == "-gmode"){i++; ScoreMatrix::TRANS_MODE = string2binary<int>(*i);
                }else if(*i == "-oldeer"){globalOption.optOLDEER = true;
                }else if(*i == "-m"){i++; globalOption.smatrixName = *i;
                }else if(*i == "-tm"){i++; globalOption.tmatrixName = *i;
                }else if(*i == "-rm"){i++; globalOption.rmatrixName = *i;
                }else if(*i == "-pm"){i++; globalOption.pmatrixName = *i;
                }else if(*i == "-averageExtendedAA"){ScoreMatrix::averageExtendedAA = true;
                }else if(*i == "-e"){i++; ScoreMatrix::epsilon = string2binary<double>(*i);
                }else if(*i == "-beta"){i++; globalOption.beta = string2binary<double>(*i);
                }else if(*i == "-beta2"){i++; globalOption.betaTQ = string2binary<double>(*i);
                }else if(*i == "-gamma"){i++; ScoreMatrix::gamma = string2binary<double>(*i);
                }else if(*i == "-pf"){i++; globalOption.mtrap_degree = 1.0 - string2binary<double>(*i);
                }else if(*i == "-withoutT"){ScoreMatrix::gapWithoutTscore = true;
                }else if(*i == "-filelist"){i++; loadFilelist(*i);
                }else if(*i == "-outbintm"){i++; globalOption.outputFilenameList.push_back(*i); globalOption.binoutMode = true;
                }else if(*i == "-search"){i++; globalOption.inputFilenameList.push_back(*i); i++; globalOption.inputFilenameList.push_back(*i); globalOption.searchMode = true;
                }else if(*i == "-search_allpair"){i++; globalOption.inputFilenameList.push_back(*i); globalOption.search_allPairMode = true;
                }else if(*i == "-searchlimit"){i++; globalOption.limitQuery = string2binary<int>(*i); i++; globalOption.limitDatabase = string2binary<int>(*i);
                }else if(*i == "-primarylibrary"){
                    i++; globalOption.primaryLibraryMode = true;
                    if(*i == "gid"){globalOption.primaryLibraryDistance = Distance::GID;
                    }else if(*i == "lid"){globalOption.primaryLibraryDistance = Distance::LID;
                    }else if(*i == "pdist"){globalOption.primaryLibraryDistance = Distance::PDIST;
                    }else if(*i == "eer"){globalOption.primaryLibraryDistance = Distance::EER;
                    }else if(*i == "eer2"){globalOption.primaryLibraryDistance = Distance::EER2;
                    }else if(*i == "ecd"){globalOption.primaryLibraryDistance = Distance::ECD;
                    }else if(*i == "ecdr"){globalOption.primaryLibraryDistance = Distance::ECDR;
                    }else{
                        putError(*i + " mode is not implemented", false);
                        return 1;
                    }
                }else if(*i == "-count"){
                    i++; if(i==vecARGV.end()) { putError("command line may have some mistakes.", true); }
                    for(int j=0; j<(int)i->size(); j++)
                    {
                        if(i->at(j) == '0'){
                            ScoreMatrix::countModeSkipGapGapSite = false;
                            ScoreMatrix::countModeSetZeroGapSite = false;
                            ScoreMatrix::countModeDirectProb = false;
                        }else if(i->at(j) == '1'){
                            ScoreMatrix::countModeSkipGapGapSite = true;
                        }else if(i->at(j) == '2'){
                            ScoreMatrix::countModeSetZeroGapSite = true;
                        }else if(i->at(j) == '3'){
                            ScoreMatrix::countModeDirectProb = true;
                        }else if(i->at(j) == '4'){
                            ScoreMatrix::countModeWithoutLowerCase = true;
                        }else if(i->at(j) == '5'){
                            ScoreMatrix::countModeGlobalTQ = true;
                        }else if(i->at(j) == 'f'){
                            ScoreMatrix::countModeSymmetrizeFreq = true;
                        }else if(i->at(j) == 'n'){
                            ScoreMatrix::countModeNoRound = true;
                        }
                    }
                    // input: database
                    i++; if(i==vecARGV.end()) { putError("command line may have some mistakes.", true); }
                    globalOption.inputFilenameList.push_back(*i);
                    // input: sequence weight
                    i++; if(i==vecARGV.end()) { putError("command line may have some mistakes.", true); }
                    globalOption.inputFilenameList.push_back(*i);
                    // output
                    i++; if(i==vecARGV.end()) { putError("command line may have some mistakes.", true); }
                    globalOption.outputFilenameList.push_back(*i);
                    globalOption.countMode = true;
                }else if(*i == "-I"){
                    ++i;
                    // reservation for first priority.
                }else{
                    return 1;
                }
            }else{
                if(globalOption.inputFilenameList.size() == 0){
                    globalOption.inputFilenameList.push_back(*i);
                }else if(globalOption.outputFilenameList.size() == 0){
                    globalOption.outputFilenameList.push_back(*i);
                }else{
                    putError("too many filenames.", false);
                    return 1;
                }
            }
        }
    } // end of priority

    // === error check ===
    if(globalOption.searchMode)
    {
        assert(globalOption.inputFilenameList.size() == 2 && globalOption.outputFilenameList.size() == 1);
    }else if(globalOption.search_allPairMode){
        assert(globalOption.inputFilenameList.size() == 1 && globalOption.outputFilenameList.size() == 1);
    }else if(globalOption.primaryLibraryMode){
        assert(globalOption.inputFilenameList.size() == 1 && globalOption.outputFilenameList.size() == 1);
    }else if(globalOption.countMode){
        assert(globalOption.inputFilenameList.size() == 2 && globalOption.outputFilenameList.size() == 1);
    }else if(globalOption.binoutMode){
        assert(globalOption.inputFilenameList.size() == 0 && globalOption.outputFilenameList.size() == 1);
    }else{
        if(globalOption.inputFilenameList.size() == 0)
        {
            cerr << "no input filename." << endl;
            return 1;
        }
        if(globalOption.outputFilenameList.size() == 0)
        {
            cerr << "no output filename." << endl;
            return 1;
        }
        assert(globalOption.inputFilenameList.size() == globalOption.outputFilenameList.size());
    }


    if(globalOption.familyEstimation)
    {
        // --- load { SCOP sccs: SABmark group } information ---
        globalOption.setMapSccsGroup();
        // --- load { SABmark sequence: SABmark group } information ---
        globalOption.setMapSeqGroup();
    }

    // --- load matrices ---
    setMatrixConfiguration();

    return 0;
}


void setMatrixConfiguration()
{
    static string previous_smatrix = "";
    static string previous_rmatrix = "";
    static string previous_pmatrix = "";
    static string previous_tmatrix = "";

    bool shouldRecalcDist = false;

    // load ramachandran matrix
    if(previous_rmatrix != globalOption.rmatrixName && ScoreMatrix::gamma > 0.0)
    {
        ScoreMatrix::setRamachandran(fixFilePath(globalOption.rmatrixName));
        previous_rmatrix = globalOption.rmatrixName;
    }

    // load substitution matrix used for partition function
    // this matrix should NOT be normalized
    if(previous_pmatrix != globalOption.pmatrixName)
    {
        ScoreMatrix::loadMatrix(fixFilePath(globalOption.pmatrixName), false, ScoreMatrix::pmatrix);
        previous_pmatrix = globalOption.pmatrixName;
        shouldRecalcDist = true;
    }

    // load substitution matrix
    if(previous_smatrix != globalOption.smatrixName)
    {
        ScoreMatrix::loadMatrix(fixFilePath(globalOption.smatrixName), true, ScoreMatrix::amatrix);
        previous_smatrix = globalOption.smatrixName;
        shouldRecalcDist = true;

        // normalize gap scores
        if(ScoreMatrix::normalized)
        {
            globalOption.gap_open_metric = ScoreMatrix::score_to_dist(globalOption.gap_open, 1);
            globalOption.gap_ext_metric = ScoreMatrix::score_to_dist(globalOption.gap_ext, 1);
            shouldRecalcDist = true;
        }
    }

    // load transition quantity matrix
    if(globalOption.tmatrixName != "" && previous_tmatrix != globalOption.tmatrixName)
    {
        ScoreMatrix::loadTMatrix(fixFilePath(globalOption.tmatrixName));
        previous_tmatrix = globalOption.tmatrixName;
        shouldRecalcDist = true;
    }

    if(shouldRecalcDist)
    {
        ScoreMatrix::calcDistTQplusDiffTable(globalOption);
    }

    // --- DEBUG ---
#ifdef DEBUG
    putLog("::setOption():");
    putLog("score metric matrix (amatrix):");
    ScoreMatrix::showScoreMatrix(ScoreMatrix::amatrix);
    putLog("score metric matrix (pmatrix):");
    ScoreMatrix::showScoreMatrix(ScoreMatrix::pmatrix);
    putLog("gap open = " + DS(globalOption.gap_open));
    putLog("gap ext  = " + DS(globalOption.gap_ext));
    putLog("gap open metric = " + DS(globalOption.gap_open_metric));
    putLog("gap ext  metric = " + DS(globalOption.gap_ext_metric));
    putLog("beta = " + DS(globalOption.beta));
    putLog("betaTQ = " + DS(globalOption.betaTQ));
    putLog("exp(beta * gap open) = " + DS(exp(globalOption.beta * globalOption.gap_open)));
    putLog("exp(beta * gap ext)  = " + DS(exp(globalOption.beta * globalOption.gap_ext)));
    putLog("exp(-beta * gap open) = " + DS(exp(-globalOption.beta * globalOption.gap_open_metric)));
    putLog("exp(-beta * gap ext)  = " + DS(exp(-globalOption.beta * globalOption.gap_ext_metric)));
    putLog("TRANS_A=" + DS(ScoreMatrix::TRANS_A));
    putLog("TRANS_B=" + DS(ScoreMatrix::TRANS_B));
#endif
}

/* Profile version of Multiple Sequence Alignment
* outFile == "" の時は標準出力にアライメントを出力する
* profileは２次元配列を１次元で管理したもの．サイズは(seq1len+1)*(seq2len+1)．また，profile[i][0]=profile[0][j]=0である．
* sparseMatrixはprofileの２次元管理版．０に近い成分は切り捨てたもの
* TODO: ペアワイズアライメント時はProfileを生成せずに直接アライメントを出力し高速化をすべき
*/
void runMSA(const string& inFile, const string& outFile)
{
    clock_t t1, t2, t3, t4;
    t1 = clock();

    // --- read the input sequences ---
    const Sequences* pInputSeqs = new Sequences(inFile, true);
    const int numSeqs = pInputSeqs->seqs.size();

    // --- translate protein coding DNA sequence to AAs
    const Sequences* pCDSSeqs = ((globalOption.CDSMode) ? pInputSeqs->genTranslatedSeqs(globalOption.CDS_type) : NULL);
    if(globalOption.CDSMode) putLog("CDS translated.");

    const Sequences& refSeqs = ((pCDSSeqs != NULL) ? *pCDSSeqs : *pInputSeqs);

    // --- generate the sparse profiles and distances ---
    // 右上三角成分のみ利用(numSeqs*(numSeqs-1)/2だけ余分にメモリを確保している)
    //putLog("generate the sparse profiles and distances");
    VVPROFILE distances(numSeqs, VPROFILE(numSeqs, PROFILE_ZERO));
    VVpSM sparseProfiles(numSeqs, VpSM(numSeqs, (SparseMatrix*)NULL));

    // --- estimate family ---
    if(globalOption.familyEstimation)
    {
        string fmatrixName = globalOption.smatrixName;
        VS estimatedGroups;
        const string mostFamousGroup = ScoreMatrix::estimateGroups(estimatedGroups, refSeqs, true);
        if(mostFamousGroup != "")
        {
            fmatrixName = fixFilePath(globalOption.fmatrixDir + "/" + mostFamousGroup + ".mat");
            putLog("family is estimated. use family specific matrix [" + fmatrixName + "]");
        }else{
            putLog("family is not estimated. use default matrix [" + fmatrixName + "]");
        }

        // --- load family specific matrix ---
        const string old_smatrix = globalOption.smatrixName;
        globalOption.smatrixName = fmatrixName;
        setMatrixConfiguration();
        globalOption.smatrixName = old_smatrix;
    }

    // --- output the parameters ---
//    putLog("=== ALIGNMENT PARAMETERS ===");
//    putLog("-gap_open = " + DS(globalOption.gap_open) + "(" + DS(globalOption.gap_open_metric) + ")");
//    putLog("-gap_ext  = " + DS(globalOption.gap_ext) + "(" + DS(globalOption.gap_ext_metric) + ")");
//    putLog("");

    // --- Open MP ---
#ifdef _OPENMP
    globalOption.numProfilePairs = numSeqs * (numSeqs - 1) / 2;
    globalOption.profilePairs = new pair<int,int>[globalOption.numProfilePairs];
    int pairIndex = 0;
    for(int i=0; i<numSeqs-1; i++)
    {
        for(int j=i+1; j<numSeqs; j++)
        {
            globalOption.profilePairs[pairIndex].first = i;
            globalOption.profilePairs[pairIndex].second = j;
            pairIndex++;
        }
    }
#endif
#ifdef _OPENMP
#pragma omp parallel for private(pairIndex) default(shared) schedule(dynamic)
    for(pairIndex=0; pairIndex<globalOption.numProfilePairs; pairIndex++)
    {
        int i = globalOption.profilePairs[pairIndex].first;
        int j = globalOption.profilePairs[pairIndex].second;
#else
    for (int i=0; i<numSeqs-1; i++)
    {
        for (int j=i+1; j<numSeqs; j++)
        {
#endif
            VPROFILE* profile_mtrap(NULL);
            VPROFILE* profile_postprob(NULL);

            string& seq1 = *(refSeqs.seqs[i]);
            string& seq2 = *(refSeqs.seqs[j]);
            string& name1 = *(refSeqs.names[i]);
            string& name2 = *(refSeqs.names[j]);
            // --- compute profile and distance ---
            // MTRAP
            AlignerMTRAP aligner_mtrap(&name1, &seq1, &name2, &seq2, &globalOption);
            PROFILE_TYPE distance_mtrap = aligner_mtrap.calcPathMatrix();
            profile_mtrap = aligner_mtrap.genProfile();

            // MTRAP with partition function
            if(fabs(globalOption.mtrap_degree - 1.0) > 1.0e-10)
            {
#ifdef USE_PROBALIGN
                // Probalign
                Probalign aligner_probalign;
                profile_postprob = aligner_probalign.genProfile(&seq1, &seq2);
#else
                // Transition quantity with partition function posterior probability
                // To use beta2=0.0 is similar to Probalign
                profile_postprob = aligner_mtrap.genPartProfile();
#endif
            }

#ifdef DEBUG
            aligner_mtrap.showPathMatrix();
            string* dAlignStr = aligner_mtrap.genAlignment();
            putLog("alignment: " + *dAlignStr);
            delete dAlignStr;
            putLog("profile (mtrap):");
            showProfile(profile_mtrap, seq1.size(), seq2.size());
            if(fabs(globalOption.mtrap_degree - 1.0) > 1.0e-10)
            {
                putLog("profile (postprob):");
                showProfile(profile_postprob, seq1.size(), seq2.size());
            }
#endif
            // --- combine all profiles and generate the sparse profile ---
            if(fabs(globalOption.mtrap_degree - 1.0) > 1.0e-10)
            {
                VPROFILE::iterator ptr_mtrap = profile_mtrap->begin();
                VPROFILE::iterator ptr_postprob = profile_postprob->begin();
                // TODO: Eigenにマップして合成すれば高速化になるはず
                for(size_t ii=0; ii<seq1.size(); ii++)
                {
                    for(size_t jj=0; jj<seq2.size(); jj++)
                    {
                        *ptr_mtrap = globalOption.mtrap_degree * (*ptr_mtrap) + (1.0 - globalOption.mtrap_degree) * (*ptr_postprob);
                        ptr_mtrap++;
                        ptr_postprob++;
                    }
                }
            }

            // consider the size of prefix '@'
            sparseProfiles[i][j] = new SparseMatrix(seq1.size()-1, seq2.size()-1, *profile_mtrap);
            // --- set the distance ---
            distances[i][j] = distance_mtrap;
            // --- finalize ---
            delete profile_mtrap;
            if(fabs(globalOption.mtrap_degree - 1.0) > 1.0e-10)
            {
                delete profile_postprob;
            }
#ifndef _OPENMP
        }
#endif
    }

    // --- create the guide tree ---
    t2 = clock();
#ifdef DEBUG
    putLog("create the guide tree");
#endif
    tr_node* pGuideTree = genTreeUPGMA(distances);
#ifdef DEBUG
    putLog(pGuideTree->toString());
#endif
    // --- calculate the sequeces weights ---
#ifdef DEBUG
    putLog("calculate the sequences weights");
#endif
    VPROFILE seqsWeights(numSeqs);
    pGuideTree->setClustalWeight(seqsWeights);
#ifdef DEBUG
    for(size_t i=0; i<seqsWeights.size(); i++) putLog("weight of leaf" + binary2string<int>(i) + "=" + binary2string<double>(seqsWeights[i]));
#endif

    // --- perform the consistency transformation ---
#ifdef DEBUG
    putLog("perform the consistency tranformation");
#endif
    for (int r=0; r<globalOption.numConsistencyReps; r++)
    {
#ifdef DEBUG
        putLog("consistency transform " + IS(r+1) + ".");
#endif
        doConsistencyTrans(seqsWeights, &refSeqs, sparseProfiles);
    }

    // --- compute the final multiple alignment ---
    t3 = clock();
    //putLog("compute the final multiple alignment");
    Sequences* finalAlignment = genFinalAlignment(pGuideTree, &refSeqs, sparseProfiles, &seqsWeights, globalOption);

    // --- output ---
    // console
    //finalAlignment->output(std::cout, ((globalOption.CDSMode) ? const_cast<Sequences*>(pInputSeqs) : NULL));
    // file
    if(outFile != "")
    {
        ofstream ofs(outFile.c_str());
        finalAlignment->output(ofs, ((globalOption.CDSMode) ? const_cast<Sequences*>(pInputSeqs) : NULL));
        ofs.close();
    }

    // --- finalize ---
    for (int i=0; i<numSeqs-1; i++)
    {
        for (int j=i+1; j<numSeqs; j++)
        {
            delete sparseProfiles[i][j];
        }
    }

    delete pInputSeqs;
    delete pCDSSeqs;
    delete pGuideTree;

    delete finalAlignment;

#ifdef _OPENMP
    delete [] globalOption.profilePairs;
#endif

    t4 = clock();
    cerr << "generate the profiles  :" << (double)(t2 - t1) / CLOCKS_PER_SEC << endl;
    cerr << "generate the guide tree:" << (double)(t3 - t2) / CLOCKS_PER_SEC << endl;
    cerr << "iterative refinements  :" << (double)(t4 - t3) / CLOCKS_PER_SEC << endl;
    cerr << "all                    :" << (double)(t4 - t1) / CLOCKS_PER_SEC << endl;
}


int main(int argc, char **argv)
{
    // === INITIALIZE ===
#ifdef _OPENMP
    Eigen::initParallel();
    Eigen::setNbThreads(0);
#endif

    initFileSearchDir();
    ScoreMatrix::init();
    setOption(argc, argv, 0);

    // === OPTION ===
    if(setOption(argc, argv, 1) > 0)
    {
        printHelp();
        return 1;
    }

    // === OUTPUT THE INPUT PARAMETERS ===
    putLog("=== SEARCH DIRS ===");
    for(VS::iterator i=fileSearchDirList.begin(); i!=fileSearchDirList.end(); ++i) putLog(*i);
    putLog("");

    if(globalOption.searchMode)
    {
        cout << "query: " << globalOption.inputFilenameList.at(0) << endl;
        cout << "database: " << globalOption.inputFilenameList.at(1) << endl;
        cout << "output: " << globalOption.outputFilenameList.at(0) << endl;
        search(globalOption.inputFilenameList.at(0), globalOption.inputFilenameList.at(1), globalOption.outputFilenameList.at(0));
    }else if(globalOption.search_allPairMode){
        cout << "database: " << globalOption.inputFilenameList.at(0) << endl;
        cout << "output: " << globalOption.outputFilenameList.at(0) << endl;
        search_allPair(globalOption.inputFilenameList.at(0), globalOption.outputFilenameList.at(0));
    }else if(globalOption.primaryLibraryMode){
        cout << "database: " << globalOption.inputFilenameList.at(0) << endl;
        cout << "output: " << globalOption.outputFilenameList.at(0) << endl;
        makePrimaryLibrary(globalOption.inputFilenameList.at(0), globalOption.outputFilenameList.at(0));
    }else if(globalOption.countMode){
        cout << "database: " << globalOption.inputFilenameList.at(0) << endl;
        cout << "sequence weight: " << globalOption.inputFilenameList.at(1) << endl;
        cout << "output: " << globalOption.outputFilenameList.at(0) << endl;
        ScoreMatrix::generateMatrix(globalOption.inputFilenameList.at(0), globalOption.inputFilenameList.at(1), globalOption.outputFilenameList.at(0));
    }else if(globalOption.binoutMode){
        cout << "generate binary format Transition-quantity matrix." << endl;
        ScoreMatrix::outputBinaryTQ( globalOption.outputFilenameList.at(0) );
    }else{
        for(int i=0; i<(int)globalOption.inputFilenameList.size(); i++)
        {
#ifdef _OPENMP
            //set OpenMP to use dynamic number of threads which is equal to the number of processor cores on the host
            //omp_set_num_threads(2);
            cerr << "OpenMP : Enabled (Max # of threads = " << omp_get_max_threads() << ")" << endl;
#endif
            runMSA(globalOption.inputFilenameList.at(i), globalOption.outputFilenameList.at(i));
        }
    }

    return 0;
}
