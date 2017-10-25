#include <string>
#include <iostream>
#include "utility.h"
#include "sequences.h"
#include "scorematrix.h"

using namespace std;

Sequences::Sequences () {
    names.clear();
    seqs.clear();
    index.clear();
}

Sequences::Sequences (const std::string& filename, bool toUpperCase) {
    VS s, n;

    names.clear();
    seqs.clear();
    index.clear();

    getMultiFASTA(n, s, filename, toUpperCase);
    for(size_t i=0; i<n.size(); i++)
    {
        names.push_back(pString(new string(n[i])));
        seqs.push_back(pString(new string('@' + s[i])));
        index.push_back(i);
    }

    // --- encode ---
    for(vec_pString::iterator s=seqs.begin(); s!=seqs.end(); ++s) ScoreMatrix::encode(**s);
}

void Sequences::addSequence(const pString& name, const pString& seq, int i) {
    seqs.push_back(seq);
    names.push_back(name);
    index.push_back(i);
}

// indexの番号順にソートする
void Sequences::sort() {
    vec_pString tn(names);
    vec_pString ts(seqs);
    for(size_t i=0; i<names.size(); i++)
    {
        names[index[i]] = tn[i];
        seqs[index[i]] = ts[i];
        index[i] = i;
    }
}


vec_pString::iterator Sequences::getIterOfSeq(int idx) {
    int i = distance(this->index.begin(), find(this->index.begin(), this->index.end(), idx));
    return this->seqs.begin() + i;
}


void Sequences::output(ostream& os, Sequences* mappingCDS) const {
    const ScoreMatrix::EncodedResidue EGAP = ScoreMatrix::encode(ScoreMatrix::GAP);

    for(size_t i=0; i< this->names.size(); i++)
    {
        // --- output the header ---
        os << '>' << *names[i] << endl;

        // --- output the sequence without prefix '@' ---
        string encodedAAs = (*seqs[i]).substr(1,(*seqs[i]).size()-1);
        string& outputSeq = encodedAAs;

        // mapping protein coding DNA sequence
        string gappedCDS = "";
        if(mappingCDS)
        {
            const int seqIdx = this->index[i];
            const string& CDS = **mappingCDS->getIterOfSeq(seqIdx);
            int CDSIndex = 1; // @

            for(string::const_iterator encodedAA = encodedAAs.begin(); encodedAA != encodedAAs.end(); ++encodedAA)
            {
                if(*encodedAA != EGAP)
                {
                    if(CDSIndex+3 <= CDS.size())
                    {
                        gappedCDS += CDS.substr(CDSIndex,3);
                    }else{
                        // if the length of CDS is not a multiple of 3, then the last incomplete codon will be appended to the last
                        gappedCDS += CDS.substr(CDSIndex);
                    }
                    CDSIndex += 3;
                }else{
                    gappedCDS += string(3,EGAP);
                }
            }

            outputSeq = gappedCDS;
        }

        ScoreMatrix::decode(outputSeq);
        for(size_t i=0; i<outputSeq.size(); i+=60)
        {
            os << outputSeq.substr(i,60) << endl;
        }
    }
}

Sequences* Sequences::genTranslatedSeqs(const VI& trans_type) const {
    Sequences* result = new Sequences();
    assert(result);

    // --- analyze genetic code ---
    // TODO
    VI genetic_codes;
    if(trans_type.size() == 1)
    {
        genetic_codes = VI(this->seqs.size(), trans_type[0]);
    }else if(trans_type.size() == this->seqs.size()){
        genetic_codes = VI(trans_type);
    }else{
        putError("wrong format: the number of genetic codes may be wrong.");
    }

    // --- translate ---
    for(std::size_t i=0; i<this->seqs.size(); i++)
    {
        const string& cds = this->seqs[i]->substr(1,this->seqs[i]->size()-1);
        const pString translatedAAs(new string(ScoreMatrix::encode('@') + ScoreMatrix::translate_CDS_to_AA(*this->names[i], cds, ScoreMatrix::geneticCode_AAs[genetic_codes[i]])));
        result->addSequence(this->names[i], translatedAAs, this->index[i]);
    }

    return result;
}


// probcons::MultiSequence::Project()相当
// indexSetで指定した配列群で新たな配列群を作成する
// 入力配列はアライメント済みであることを想定
// ギャップのみのサイトをカットする
Sequences* Sequences::extract(const VI& indexSet) {
    const ScoreMatrix::EncodedResidue GAP = ScoreMatrix::encode(ScoreMatrix::GAP);
    const int oldAlignmentSize = (*this->seqs[0]).size();
    Sequences* newAlignment = new Sequences();

    // generate the flags which sites are gapped
    vector<bool> gapSite(oldAlignmentSize, true);
    int numGapSite=0;
    // include the prefix '@'
    for(int c=0; c<oldAlignmentSize; c++)
    {
        for(VI::const_iterator i=indexSet.begin(); i!=indexSet.end(); ++i)
        {
            if((*this->seqs[*i])[c] != GAP)
            {
                gapSite[c] = false;
                numGapSite++;
                break;
            }
        }
    }

    if(numGapSite > 0)
    {
        // push_back the gap-site omitted sequence
        for(VI::const_iterator i=indexSet.begin(); i!=indexSet.end(); ++i)
        {
            string* newStr = genOmittedString(*this->seqs[*i], gapSite);
            newAlignment->addSequence(this->names[*i], pString(newStr), this->index[*i]);
        }
    }else{
        for(size_t i=0; i<indexSet.size(); i++)
        {
            newAlignment->addSequence(this->names[indexSet[i]], this->seqs[indexSet[i]], this->index[indexSet[i]]);
        }
    }

    return newAlignment;
}


// Returns a SafeVector<int> containing the indeces of every character in the sequence.
// For instance, if the data is "@ATGCC---GT--CA", the method returns
// {0,1,2,3,4,5,9,10,13,14}.
// probcons::Sequence::GetMapping()相当
VI* genNonGapIndexList (const std::string& s) {
    const ScoreMatrix::EncodedResidue GAP = ScoreMatrix::encode(ScoreMatrix::GAP);
    VI* ret = new VI(1, 0);
    // skip the prefix '@'
    for(size_t i=1; i < s.size(); i++)
    {
      if (s[i] != GAP) ret->push_back(i);
    }
    return ret;
}

// probcons::Sequence::AddGaps()相当
// seq = "@ATGCAGTCA"
// alignment = "XXXBBYYYBBYYXX"
// id = 'X'
// ->
// "@ATGCC---GT--CA"
// " XXXBBYYYBBYYXX"
std::string* genGappedSeq(const std::string& seq, const std::string& alignment, char id) {
    const ScoreMatrix::EncodedResidue GAP = ScoreMatrix::encode(ScoreMatrix::GAP);
    std::string* ret = new string();
    assert(ret);

    std::string::const_iterator seqi = seq.begin();
    ++seqi; ret->push_back(ScoreMatrix::encode('@'));
    for(std::string::const_iterator i=alignment.begin(); i!=alignment.end(); ++i)
    {
        if(*i == 'B' || *i == id)
        {
            ret->push_back(*seqi);
            ++seqi;
        }else{
            ret->push_back(GAP);
        }
    }

    return ret;
}

