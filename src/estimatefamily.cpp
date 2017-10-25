#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <map>
#include <set>
#include "utility.h"
#include "globaloption.h"
#include "sequences.h"

using namespace std;

namespace ScoreMatrix {

std::string findSCOP_sccs(const std::string& s){
    size_t dotPos = -1;

    // --- search '.' ---
    for(int i=0;i<3;i++)
    {
        size_t dotPos = s.find('.',dotPos+1);
        if(dotPos == std::string::npos) break;

        // --- check format ---
        if(0 < dotPos && dotPos+5 < s.size() && s[dotPos+2] == '.' && s[dotPos+4] == '.')
        {
            return s.substr(dotPos-1,7);
        }
    }

    return "";
}

std::string estimateGroups(VS& group_list, const Sequences& querySeqs, bool findMostFamous){
    group_list.clear();

    // --- [FIRST PRIORITY] search SCOP sccs from the name header of a sequence ---
    {
        VS scop_sccs_list;

        putLog("search SCOP sccs...");

        // --- search SCOP sccs ID from the header ---
        for(vec_pString::const_iterator i=querySeqs.names.begin(); i!=querySeqs.names.end(); ++i)
        {
            // if sccs is not found then the "" is stored
            scop_sccs_list.push_back( findSCOP_sccs(**i) );
        }
        assert(querySeqs.names.size() == scop_sccs_list.size());

        // --- search the group name from SCOP sccs w.r.t. each sequence ---
        for(VS::iterator sccs=scop_sccs_list.begin(); sccs!=scop_sccs_list.end(); ++sccs)
        {
            // if sccs is not found or group name of the sccs is not found then the "" is stored
            if(*sccs != "")
            {
                if(globalOption.mapID_sccs_group.find(*sccs) != globalOption.mapID_sccs_group.end())
                {
                    group_list.push_back(globalOption.mapID_sccs_group[*sccs]);
                }else{
                    group_list.push_back("");
                }
            }else{
                group_list.push_back("");
            }
        }
    }

    // --- [SECOND PRIORITY] estimate group name by BLAST search ---
    // apply this step when the all sequences do not be linked at the first priority step
    if(count(group_list.begin(),group_list.end(),"") == group_list.size()){
        putLog("estimate group by BLAST...");

        const static int BUF = 256;
        const static string tmpFile = "_tmp";
        string errStr = "";
        typedef std::pair<int,string> pair_queryindex_subject;
        multimap<double,pair_queryindex_subject> mapEvalueSABSeq;

        // --- make a temporary query file based on the all sequences in the alignemnt ---
        // TODO: omit the queries that is already linked with the group name, for the speed up
        ofstream ofs(tmpFile.c_str());
        querySeqs.output(ofs);
        ofs.close();

        // --- search the most similar sequence from the database ---
        FILE	*fp;
        char	buf[BUF];
        const string cmdline = "blastp -evalue 1.0e-3 -db " + fixFilePath(globalOption.familyDatabaseName)
                + " -query " + tmpFile
                + " -outfmt \"6 sacc evalue\" -num_alignments 1";

        putLog("RUN: " + cmdline);

        if ( (fp=popen(cmdline.c_str(),"r")) == NULL)
        {
            errStr = "blast+ not found.";
        }else{
            int qidx = 0;
            while(fgets(buf, BUF, fp) != NULL) {
                VS line;
                mySplit(string(buf), "\t", line);
                // NOTE: subject name is stripped with length less than 10 by BLAST
                mapEvalueSABSeq.insert( std::pair<double,pair_queryindex_subject>(SD(line[1]),pair_queryindex_subject(qidx,line[0])) );
                qidx++;
            }
            pclose(fp);
        }

        // --- finalize ---
        remove(tmpFile.c_str());
        if(errStr != "")
        {
            putError(errStr, true);
        }

        // --- put estimated group name for the sequence that is not linked with the group name yet ---
        putLog("mearge estimated information...");

        for(multimap<double, pair_queryindex_subject>::iterator i=mapEvalueSABSeq.begin(); i!=mapEvalueSABSeq.end(); ++i)
        {
            const double evalue = i->first;
            const string &queryName = *querySeqs.names[i->second.first];
            const string &subjectName = i->second.second;
            const int index = getIndex(querySeqs.names,queryName);
            if(group_list[index] != "") group_list[index] = globalOption.mapID_seq_group[subjectName];
        }
    }

    // --- [FINALIZE STEP] black list filtering ---
    // if the group name is in a black list then the entry is replaced by ""
    for(VS::iterator i=group_list.begin(); i!=group_list.end(); ++i)
    {
        // black list filtering
        if(globalOption.groupBlackList.find(*i) != globalOption.groupBlackList.end()
                && globalOption.groupWhiteList.find(*i) == globalOption.groupWhiteList.end() )
        {
            cout << *i << " is omitted." << endl;
            *i = "";
        }
    }

    // --- [OPTION STEP] find the most famous ---
    if(findMostFamous)
    {
        putLog("find most famous group...");

        map<string,int> mapGroupFreq;
        for(VS::iterator i=group_list.begin(); i!=group_list.end(); ++i)
        {
            if(*i != "") mapGroupFreq[*i]++;
        }

        string maxGroup = "";
        int maxFreq = -1;
        for(map<string,int>::iterator i=mapGroupFreq.begin(); i!=mapGroupFreq.end(); ++i)
        {
            if(i->second > maxFreq)
            {
                maxGroup = i->first;
                maxFreq = i->second;
            }
        }

        return maxGroup;
    }

    return "";
}


} // end of namespace
