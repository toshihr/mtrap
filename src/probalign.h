#ifndef PROBALIGN_H
#define PROBALIGN_H

#include <string>
#include "SafeVector.h"
#include "SparseMatrix.h"
#include "scorematrix.h"

#ifdef _WIN32
#define PROBALIGN_HUGE_VALL	HUGE_VAL
#else
#define PROBALIGN_HUGE_VALL	HUGE_VALL
#endif


class Probalign
{
private:
    typedef struct {
        std::string input;
        int matrix;
        int N;
        float T;
        float beta;
        char opt;			//can be 'P' or 'M'
        float gapopen;
        float gapext;
    } argument_decl;

    argument_decl argument;

    // --- options ---
    static const bool endgaps = false; // デフォルトはtrue
    static const bool PART_FULL_MEMORY = true;
    static const bool REVPART_FULL_MEMORY = true;

    float g_gap_open1, g_gap_open2, g_gap_ext1, g_gap_ext2;
    char aminos[26], matrixtype[20], bases[26];

    ScoreMatrix::AminoMatrix sub_matrix; // unnormalized matrix

    void init_arguments();
    void read_matrix(const std::string& matrixName, float beta);
    VPROFILE* doReverse(const SafeVector< std::string* > & sequences, const double termgapopen, const double termgapextend, long double **Zfm, const double d, const double e);
    long double** doForward(const SafeVector< std::string* > & sequences, const double termgapopen, const double termgapextend, const double d, const double e);

    // should not be implemented for omitting the unexpected behavior
    Probalign( const Probalign& rhs);
    Probalign& operator=(const Probalign& rhs);
public:
    Probalign();
    ~Probalign();

    Probalign* operator& (){ return this; };
    const Probalign* operator& () const{ return this; };

    VPROFILE* genProfile(std::string* seq1, std::string* seq2);

    // --- DEBUG ---
    void showPassMatrix(long double **Zm, long double **Zf, long double **Ze, int numRow, int numCol, std::string& A, std::string& B);
};

#endif // PROBALIGN_H
