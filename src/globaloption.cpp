#include <vector>
#include <set>
#include <string>
#include <iostream>
#include <fstream>
#include "globaloption.h"
#include "utility.h"

using namespace std;

void GlobalOption::setMapSccsGroup(){
    string buff;

    // --- init ---
    this->mapID_sccs_group.clear();

    if(this->csvSCOP_SAB_ID_Name != "")
    {
        // --- open ---
        ifstream ifs(fixFilePath(this->csvSCOP_SAB_ID_Name).c_str());

        // --- read data [sunid,sccs,sid] ---
        while(ifs && getline(ifs, buff))
        {
            VS vs;
            mySplit(buff, ", ", vs);
            if(vs.size() == 4)
            {
                if(vs[2] == "-") vs[2] = "";
                if(vs[3] == "-") vs[3] = "";
                const string &sunid = vs[0];
                const string &sccs = vs[1];
                const string &sid = vs[2];
                const string &group = vs[3];
                this->mapID_sccs_group[sccs] = group;
            }
        }
    }
}


void GlobalOption::setMapSeqGroup(){
    string buff;

    // init
    this->mapID_seq_group.clear();

    if(this->csvSAB_seq_group_Name != "")
    {
        // open
        ifstream ifs(fixFilePath(this->csvSAB_seq_group_Name).c_str());

        // omit the header line
        getline(ifs, buff);

        // read data [Sequence,Group,SCOP sunid]
        while(ifs && getline(ifs, buff))
        {
            VS list;
            VS vs;
            mySplit(buff, ", ", vs);
            for(VS::iterator i=vs.begin(); i!=vs.end(); ++i)list.push_back(*i);
            if(list.size() == 3)
            {
                this->mapID_seq_group[list[0]] = list[1];
            }
        }
    }

    putLog("family information (" + IS(this->mapID_seq_group.size())  + " families) is loaded.");

    // --- set black list ---
    this->groupBlackList.clear();
    this->groupBlackList.insert("group169");
    this->groupBlackList.insert("group353");
    this->groupBlackList.insert("group357");

    // --- set white list ---
    this->groupWhiteList.clear();
    this->groupWhiteList.insert("group1");
    this->groupWhiteList.insert("group38");
    this->groupWhiteList.insert("group52");
    this->groupWhiteList.insert("group53");
    this->groupWhiteList.insert("group76");
    this->groupWhiteList.insert("group78");
    this->groupWhiteList.insert("group79");
    this->groupWhiteList.insert("group96");
    this->groupWhiteList.insert("group99");
    this->groupWhiteList.insert("group106");
    this->groupWhiteList.insert("group112");
    this->groupWhiteList.insert("group128");
    this->groupWhiteList.insert("group131");
    this->groupWhiteList.insert("group136");
    this->groupWhiteList.insert("group138");
    this->groupWhiteList.insert("group145");
    this->groupWhiteList.insert("group156");
    this->groupWhiteList.insert("group167");
    this->groupWhiteList.insert("group171");
    this->groupWhiteList.insert("group182");
    this->groupWhiteList.insert("group186");
    this->groupWhiteList.insert("group192");
    this->groupWhiteList.insert("group198");
    this->groupWhiteList.insert("group206");
    this->groupWhiteList.insert("group215");
    this->groupWhiteList.insert("group245");
    this->groupWhiteList.insert("group251");
    this->groupWhiteList.insert("group254");
    this->groupWhiteList.insert("group264");
    this->groupWhiteList.insert("group273");
    this->groupWhiteList.insert("group274");
    this->groupWhiteList.insert("group278");
    this->groupWhiteList.insert("group283");
    this->groupWhiteList.insert("group303");
    this->groupWhiteList.insert("group305");
    this->groupWhiteList.insert("group311");
    this->groupWhiteList.insert("group336");
    this->groupWhiteList.insert("group337");
    this->groupWhiteList.insert("group339");
    this->groupWhiteList.insert("group343");
    this->groupWhiteList.insert("group346");
    this->groupWhiteList.insert("group361");
    this->groupWhiteList.insert("group362");
    this->groupWhiteList.insert("group365");
    this->groupWhiteList.insert("group373");
    this->groupWhiteList.insert("group377");
    this->groupWhiteList.insert("group383");
    this->groupWhiteList.insert("group385");
    this->groupWhiteList.insert("group396");
    this->groupWhiteList.insert("group397");
}

// TODO: 将来的には廃止.mtrapにメンバとして持たせる
GlobalOption globalOption;

