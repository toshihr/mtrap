#ifndef SEQUENCES_H
#define SEQUENCES_H
#include <string>
#include "utility.h"

class Sequences {
public:
    vec_pString names;
    vec_pString seqs;
    VI index;	// the indeces of the input file

    Sequences ();
    Sequences (const std::string& filename, bool toUpperCase);
    void addSequence(const pString& name, const pString& seq, int i);
    void sort();
    void output(ostream& os, Sequences* mappingCDS = NULL) const;
    Sequences* genTranslatedSeqs(const VI& trans_type) const;
    vec_pString::iterator getIterOfSeq(int idx);

    Sequences* extract(const VI& indexSet);
};


VI* genNonGapIndexList (const std::string& s);
std::string* genGappedSeq(const std::string& seq, const std::string& alignment, char id);

#endif

