#pragma once

#include <string>
#include <map>
#include <set>
#include "scorematrix.h"
#include "distance.h"
#include "utility.h"

// --- run mode ---
typedef enum { MODE_ALIGN, MODE_SEARCH, MODE_SEARCHALLPAIR } RUNMODE;

class GlobalOption {
public:
    Distance::MODE primaryLibraryDistance;

    VS inputFilenameList;
    VS outputFilenameList;
    std::string tmatrixName;
    std::string smatrixName;
    std::string rmatrixName;
    std::string pmatrixName; // used for partition function
    std::string fmatrixDir; // used for family specific matrix
    std::string outputFilename_pass;

    // --- SCOP & SABmark ID information ---
    // sid   -- SCOP identifier. e.g. d1danl2
    // sccs  -- SCOP concise classification strings.  e.g. b.1.2.1
    // sunid -- SCOP unique identifier for this domain
    // group -- SAMmark group
    std::string csvSCOP_SAB_ID_Name;
    std::string csvSAB_seq_group_Name;
    std::map<std::string, std::string> mapID_sccs_group; // SCOP's sccs -> SABmark's group
    std::map<std::string, std::string> mapID_seq_group; // SABmark's sequence name -> SABmark's group

    // --- family information ---
    bool familyEstimation;
    std::string familyDatabaseName;
    std::set<std::string> groupBlackList;
    std::set<std::string> groupWhiteList;

    // --- global options ---

    // --- output options ---
    bool sortByInputOrder;

    // --- mode flags ---
    // TODO: to be enumerate
    bool searchMode;
    bool search_allPairMode;
    bool primaryLibraryMode;
    bool countMode;
    bool binoutMode;
    bool CDSMode;

    // --- CDS mode ---
    VI CDS_type;

    // --- align mode ---
    // global options
    int alignment_size;		//アライメントの長さ(結果)
    // MTRAP
    double gap_open; // Score
    double gap_ext; // Score
    double gap_open_metric; // Metric
    double gap_ext_metric; // Metric
    double mtrap_degree; // degree of MTRAP against MSAProbs
    int numConsistencyReps;
    int numIterativeRefReps;
    // part MTRAP
    double beta;
    double betaTQ;
    double gap_open_part; // for partition function
    double gap_ext_part; // for partition function

    // --- database search mode ---
    int limitQuery;
    int limitDatabase;

    // --- legacy code ---
    bool optOLDEER;
#ifdef _OPENMP
    pair<int,int> * profilePairs;
    int numProfilePairs;
#endif

    // --- functions ---
    void setMapSccsGroup();
    void setMapSeqGroup();
};
extern GlobalOption globalOption;

