/*!
    \file   primarylibrary.cpp
    \author Toshihide Hara
    \date   Sun Jul  3 17:39:33 2011

    \brief  generate T-Coffee style primary library.


*/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include "utility.h"
#include "scorematrix.h"
#include "pairwise.h"
#include "entropy.h"
#include "primarylibrary.h"
#include "sequences.h"

using namespace std;


//! make Primary Library of T-Coffee from input database.
/*!

    \param databaseFile input filename that contains Multiple FASTA
    \param outFile output filename
*/
// TODO: these code should be checked!!
void makePrimaryLibrary(const string& databaseFile, const string& outFile)
{
    // --- read database sequences ---
    Sequences* pSeqs = new Sequences(databaseFile, true);
    const int numSeqs = (globalOption.limitDatabase != 0) ? min(static_cast<int>(pSeqs->seqs.size()), globalOption.limitDatabase) : pSeqs->seqs.size();

    cout << endl;
    cout << "=== PROGRESS ===" << endl;

    // --- output the header of primary library ---
    ofstream ofs(outFile.c_str());
    // HEADER
    ofs << "! TC_LIB_FORMAT_01" << endl;
    // NUMBER
    ofs << numSeqs << endl;
    // SEQUENCES
    for(int i=0; i<numSeqs; i++)
    {
        ofs << *(pSeqs->names[i]);
        ofs << ' ';
        ofs << pSeqs->seqs[i]->size();
        ofs << ' ';
        ofs << *(pSeqs->seqs[i]);
        ofs << endl;
    }

    // --- output the body of primary library ---
    for (int i=0; i<numSeqs-1; i++)
    {
        for (int j=i+1; j<numSeqs; j++)
        {
            string& seq1 = *(pSeqs->seqs[i]);
            string& seq2 = *(pSeqs->seqs[j]);
            string& name1 = *(pSeqs->names[i]);
            string& name2 = *(pSeqs->names[j]);

            AlignerMTRAP aligner_mtrap(&name1, &seq1, &name2, &seq2, &globalOption);
            aligner_mtrap.calcPathMatrix();
            string* strAlignment = aligner_mtrap.genAlignment();

            vec_pString aligned;
            aligned.push_back(pString(genGappedSeq(seq1, *strAlignment, 'X')));
            aligned.push_back(pString(genGappedSeq(seq2, *strAlignment, 'Y')));

            delete strAlignment;

            vector<pair<int,int> > indexPair;
            getResiduePairsAndNumOfIdentity(indexPair, aligned.at(0), aligned.at(1));

            double distance = Distance::calcDistance<double>(*aligned.at(0), *aligned.at(1), globalOption.primaryLibraryDistance);

            // --- output the log ---
            cout << "pair: " << i << "," << j << endl;

            // --- output the primary library ---
            ofs << "#" << i+1 << " " << j+1 << endl;
            for(vector<pair<int,int> >::iterator i=indexPair.begin(); i!=indexPair.end(); ++i)
            {
                ofs << setw(5) << i->first;
                ofs << " ";
                ofs << setw(5) << i->second;
                ofs << " ";
                ofs << setw(5) << (int)(distance * 100.0 * 10.0); // for ver. >=5.31?
//				ofs << setw(5) << (int)(distance * 100.0); // for ver. <5.31?
                ofs << endl;

            }
        }
    }

    ofs << "! CPU 0" << endl;
    ofs << "! SEQ_1_TO_N" << endl;

    ofs.close();
}

// 指定されたアライメント中に存在する文字ペアをリストに登録．また，一致ペア数を返す(配列一致率での利用を想定)．
// 登録されるペアはアライメントにおけるカラム番号(origin 1)
int getResiduePairsAndNumOfIdentity(vector<pair<int,int> >& indexPair, const pString& s1, const pString& s2)
{
    const int length = s1->size();
    int index1 = 0;
    int index2 = 0;
    int same = 0;

    indexPair.reserve(length);
    indexPair.clear();

    for(int i=0; i<length; ++i)
    {
        const bool ngap1 = (s1->at(i) != ScoreMatrix::GAP);
        const bool ngap2 = (s2->at(i) != ScoreMatrix::GAP);

        if(ngap1)index1++;
        if(ngap2)index2++;
        if(ngap1 && ngap2)
        {
            indexPair.push_back(pair<int,int>(index1, index2));
            if(s1->at(i) == s2->at(i))same++;
        }
    }

    return same;
}


