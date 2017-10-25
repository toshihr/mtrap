#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include "scorematrix.h"
#include "utility.h"
#include "search.h"
#include "pairwise.h"
#include "sequences.h"
#include "entropy.h"

using namespace std;

// database search for all pairs of the database
void search_allPair(const string& databaseFile, const string& outFile)
{
    search(databaseFile, databaseFile, outFile);
}

// database search
void search(const string& queryFile, const string& databaseFile, const string& outFile)
{
    // --- read database sequences ---
    Sequences* querySeqs = new Sequences(queryFile, true);
    Sequences* databaseSeqs = new Sequences(databaseFile, true);
    const int numQuery = querySeqs->names.size();
    const int numDatabase = databaseSeqs->names.size();

    // --- initialize output file ---
    ofstream ofs(outFile.c_str());

    // --- output the log ---
    cout << endl;
    cout << "=== PROGRESS ===" << endl;

    // --- main ---
    for (int i=0; i<numQuery; i++)
    {
        // output the progress
        cout << "query:" << *(querySeqs->names[i]) << endl;

        for (int j=0; j<numDatabase; j++)
        {
            string& seq1 = *(querySeqs->seqs[i]);
            string& seq2 = *(databaseSeqs->seqs[j]);
            string& name1 = *(querySeqs->names[i]);
            string& name2 = *(databaseSeqs->names[j]);

            AlignerMTRAP aligner_mtrap(&name1, &seq1, &name2, &seq2, &globalOption);
            //const double metric = aligner_mtrap.calcPassMatrix();
            string* strAlignment = aligner_mtrap.genAlignment();

            vec_pString aligned;
            aligned.push_back(pString(genGappedSeq(seq1, *strAlignment, 'X')));
            aligned.push_back(pString(genGappedSeq(seq2, *strAlignment, 'Y')));

            delete strAlignment;

            // calc some scores
            //const double lnK = Log10(seq1.size() * seq2.size());
            //const double lnK2 = Log10(numDatabase);
            //const double lambda = -1.0;
            //const double ln2 = Log10(2);
            //const double bits = (lambda * metric - lnK) / ln2;
            //const double bits2 = (lambda * metric - lnK2) / ln2;
            //const double evalue = seq1.size() * seq2.size() * pow(2.0,-bits);
            //const double evalue2 = seq1.size() * seq2.size() * pow(2.0,-bits2);

            // calc EER
            Entropy::VT p1, p2;
            Entropy::VVT joint_p;
            Entropy::calcProb(p1, *aligned.at(0));
            Entropy::calcProb(p2, *aligned.at(1));
            Entropy::calcComplexProb(joint_p, *aligned.at(0), *aligned.at(1));
            double eer = Entropy::getEER2(p1, p2, joint_p);
            // output the results
            // ASTRAL SCOP HEADER: [sequence name 0..6]space[SCOP ID 8..14]
            assert(name1.size() > 16 && "The database may not be ASTRAL SCOP file.");
            assert(name2.size() > 16 && "The database may not be ASTRAL SCOP file.");
            string query_name = name1.substr(0,6-0+1);
            string subject_name = name2.substr(0,6-0+1);
            string query_family = name1.substr(8, 14-8+1);
            string subject_family = name2.substr(8, 14-8+1);
            ofs << ((query_name == subject_name) ? 1 : 0);
            ofs << ' ' << subject_name << ' ' << query_name;
            ofs << ' ' << ((query_family == subject_family)? 1 : 0);
            ofs << ' ' << subject_family << ' ' << query_family;
            ofs << ' ' << eer;
            ofs << ' ' << *aligned.at(0); // query
            ofs << ' ' << *aligned.at(1); // subject
            ofs << endl;

/*			ofs << "query:" << *name.at(0);
            ofs << "; target:" << *name.at(1);
            ofs << "; metric:" << metric;
            ofs << "; bits:" << bits;
            ofs << "; E:" << evalue;
            ofs << "; bits2:" << bits2;
            ofs << "; E2:" << evalue2;
            ofs << endl;
*/
        }
    }

    // === finalize ===
    ofs.close();
}
