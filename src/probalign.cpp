/* Probalign style alignment
* MTRAP版partiton functionが完成したため、このコードは本来いらない
* References:
* Probalign::PostProbs.cc, ReadMatrix.cc License=PUBLIC DOMAIN (based on the Hasegawa's paper)
*/

#include <iomanip>
#include <string>
#include <cassert>
#include "SafeVector.h"
#include "SparseMatrix.h"
#include "scorematrix.h"
#include "probalign.h"

using namespace std;


Probalign::Probalign() {
    init_arguments();
}

Probalign::~Probalign() {
}

void Probalign::read_matrix(const string& matrixName, float beta) {
    ScoreMatrix::loadMatrix(matrixName, false, this->sub_matrix);
    for(int r=0; r<ScoreMatrix::RANK; r++)
    {
        for(int c=0; c<ScoreMatrix::RANK; c++)
        {
            sub_matrix(r,c) = exp(beta * sub_matrix(r,c));
        }
    }

#ifdef DEBUG
    putLog("Probalign::read_matrix():");
    putLog("probalign matrix:");

    stringstream header;
    header << "  ";
    for(int c=0; c<ScoreMatrix::RANK; c++)
    {
        header << left << setw(5) << ScoreMatrix::RESIDUE[c];
    }
    putLog(header.str());

    for(int r=0; r<ScoreMatrix::RANK; r++)
    {
        stringstream line;
        line << ScoreMatrix::RESIDUE[r] << ' ';
        for(int c=0; c<ScoreMatrix::RANK; c++)
        {
            line << left << setw(5) << setprecision(2) << this->sub_matrix(r,c);
        }
        putLog(line.str());
    }
#endif

}

// corresponding to probalign::ReadMatrix.cc::init_arguments
void Probalign::init_arguments() {
    const float GAPOPEN = 0.0;
    const float GAPEXT = 0.0;
    const float TEMPERATURE = 5.0;
    //const int le = 4; // nucleotide
    const int le = 160; // protein
    assert(le == 4 || le == 160);

    // --- set the default values ---
    argument.N = 1;
    argument.input = "tempin";
    argument.matrix = le;
    argument.gapopen = GAPOPEN;
    argument.gapext = GAPEXT;
    argument.T = TEMPERATURE;
    argument.beta = 1.0f / TEMPERATURE;
    argument.opt = 'P';

    // --- set the gap costs ---
    float gap_open = 0, gap_ext = 0;
    if(le == 4)
    {
        gap_open = -4;
        gap_ext = -0.25;
        //TODO: read_matrix(nuc_id);
    }else if(le == 160){
        gap_open = -22;
        gap_ext = -1;
        read_matrix("CGONNET160", argument.beta);
    }
    if (argument.gapopen != 0.0 || argument.gapext != 0.0)
    {
        gap_open = -argument.gapopen;
        gap_ext = -argument.gapext;
    }
    argument.gapopen = gap_open;
    argument.gapext = gap_ext;
}


// Ai i=1,...,m
// Bj j=1,...,n
// forward:
// i: 0  1  2      m
// A: @ A1 A2 ... Am
// j: 0  1  2      n
// B: @ B1 B2 ... Bn
// backward:
// i:  0  1     m-1 m
// A: A1 A2 ... Am  @
// j:  0  1     n-1 n
// B: B1 B2 ... Bn  @
VPROFILE* Probalign::doReverse(const SafeVector< std::string* > & sequences, const double termgapopen, const double termgapextend, long double **Zfm, const double d, const double e) {
    long double **Zm = NULL;
    long double **Ze = NULL;
    long double **Zf = NULL;
    const int len0 = (*sequences[0]).size()-1; // without the prefix
    const int len1 = (*sequences[1]).size()-1;
    double probability;
    long double tempvar;
    double endgapopen, endgapextend;

    // --- initialize ---
    VPROFILE *posteriorPtr = new VPROFILE((len0 + 1) * (len1 + 1));
    VPROFILE & posterior = *posteriorPtr;
    VPROFILE::iterator ptr = posterior.begin();

    endgapopen = termgapopen;
    endgapextend = termgapextend;

    if (REVPART_FULL_MEMORY)
    {
        Ze = new long double *[len1 + 1];
        Zf = new long double *[len1 + 1];
        Zm = new long double *[len1 + 1];
        for (int i=0; i<=len1; i++)
        {
            Ze[i] = new long double[len0 + 1];
            Zf[i] = new long double[len0 + 1];
            Zm[i] = new long double[len0 + 1];
        }
    } else {
        Zm = new long double *[2];
        Ze = new long double *[2];
        Zf = new long double *[2];
        for (int i=0; i<=1; i++)
        {
            Zm[i] = new long double[len0 + 1];
            Ze[i] = new long double[len0 + 1];
            Zf[i] = new long double[len0 + 1];
        }
    }

    if (REVPART_FULL_MEMORY)
    {
        for (int i=0; i<=len1; i++)
        {
            for (int j=0; j<=len0; j++)
            {
                Zm[i][j] = 0.0;
                Zf[i][j] = 0.0;
                Ze[i][j] = 0.0;
            }
        }
    } else {
        for (int j=0; j<=len0; j++)
        {
            Zm[0][j] = 0;
            Zf[0][j] = 0;
            Ze[0][j] = 0;
            Zf[1][j] = 0;
            Ze[1][j] = 0;
            Zm[1][j] = 0;
        }
    }

    // --- main ---
    // fill the probability matrix with 0s
    for(int i=0; i<=len1; i++) for(int j=0; j<=len0; j++) ptr[j * (len1 + 1) + i] = 0;

    if(endgaps == false)
    {
        Zm[len1][len0] = 1;
        Ze[len1][len0] = Zf[len1][len0] = 0;
        Zf[len1 - 1][len0] = Zm[len1][len0] * d;
        Ze[len1][len0 - 1] = Zm[len1][len0] * d;

        // >=2ND ROW INIT
        if(REVPART_FULL_MEMORY)
        {
            for(int i=len1 - 2; i>=0; i--)
            {
                Zf[i][len0] = Zf[i + 1][len0] * e;
            }
        }

        // >=2ND COL INIT
        if(REVPART_FULL_MEMORY)
        {
            for(int j=len0 - 2; j>=0; j--)
            {
                Ze[len1][j] = Ze[len1][j + 1] * e;
            }
        } else {
            for(int j=len0 - 2; j>=0; j--)
            {
                Ze[0][j] = Ze[0][j + 1] * e;
            }
        }
    } else { // endgaps == true
        if (REVPART_FULL_MEMORY)
        {
            Zm[len1][len0] = 1;
            Ze[len1][len0] = Zf[len1][len0] = 0;
            Zf[len1 - 1][len0] = Zm[len1][len0] * endgapopen;
            Ze[len1][len0 - 1] = Zm[len1][len0] * endgapopen;

            // >=2ND ROW INIT
            for(int i=len1 - 2; i>=0; i--)
            {
                Zf[i][len0] = Zf[i + 1][len0] * endgapextend;
            }

            //M Iy= d+j*e

            // >=2ND COL INIT
            for(int j=len0 - 2; j>=0; j--)
            {
                Ze[len1][j] = Ze[len1][j + 1] * endgapextend;
            }
        } else {
            // in Zm
            // let:
            //  Zm(0) be the current row being filled/computed
            //  Zm(1) be the previous row

            Zm[1][len0] = 1;
            Ze[0][len0] = Zf[0][len0] = 0;
            Zf[1][len0] = Zm[1][len0] * endgapopen;
            Ze[0][len0 - 1] = Zm[1][len0] * endgapopen;

            // >=2ND COL INIT
            for(int j=len0 - 2; j>=0; j--)
            {
                Ze[0][j] = Ze[0][j + 1] * endgapextend;
            }
        } // REVPART_FULL_MEMORY
    } // endgaps

    for(int i=len1 - 1; i>=0; i--)
    {
        for (int j=len0 - 1; j>=0; j--)
        {
            const ScoreMatrix::EncodedResidue Si = (*sequences[1])[i+1];
            const ScoreMatrix::EncodedResidue Tj = (*sequences[0])[j+1];
            const double scorez = sub_matrix(Si,Tj);

            // endgaps modification aug 10
            double open0, extend0, open1, extend1;

            open0 = open1 = d;
            extend0 = extend1 = e;

            if(endgaps)
            {
                // check to see if one of the 2 sequences or both reach the end
                if(i==0)
                {
                    open0 = endgapopen;
                    extend0 = endgapextend;
                }
                if(j==0)
                {
                    open1 = endgapopen;
                    extend1 = endgapextend;
                }
            }

            if(REVPART_FULL_MEMORY)
            {
                // z computation
                Ze[i][j] = Zm[i][j + 1] * open0 + Ze[i][j + 1] * extend0;
                Zf[i][j] = Zm[i + 1][j] * open1 + Zf[i + 1][j] * extend1;
                Zm[i][j] = (Zm[i + 1][j + 1] + Zf[i + 1][j + 1] + Ze[i + 1][j + 1]) * scorez;
                // original probalign has bugs in these area
                // using paper index, Zm[i][j] -> Zm[i+1,j+1]
                // thus,
                // Zfm[i+1][j+1] * Zm[i+1][j+1]
                // = Z[i][j]*S(Ai+1,Bj+1) * Z[i+2][j+2] * S(Ai+1,Bj+1)
                // therefore divide scorez ( = S(Ai+1,Bj+1) ) is justified
                tempvar = Zfm[i + 1][j + 1] * Zm[i][j];
            } else {
                Zf[1][j] = Zm[1][j] * open1 + Zf[0][j] * extend1;
                Ze[1][j] = Zm[0][j + 1] * open0 + Ze[1][j + 1] * extend0;
                Zm[0][j] = (Zm[1][j + 1] + Zf[0][j + 1] + Ze[0][j + 1]) * scorez;

                tempvar = Zfm[i + 1][j + 1] * Zm[0][j];
            }
            // divide P(i,j) i.e. pairwise probability by denominator
            tempvar /= (scorez * Zfm[0][0]);
            probability = (double) tempvar;

            // store only noticable probabilities
            //if(probability <= 1.0 && probability >= 0.001)
            if(probability >= 0.001)
            {
                //algorithm goes...
                // validprob[i + 1][j + 1] = probability;
                ptr[(j + 1) * (len1 + 1) + (i + 1)] = min(1.0, probability);
            }
        } // end of for j

        if (REVPART_FULL_MEMORY == false)
        {
            for(int t=0; t<=len0; t++)
            {
                Ze[0][t] = Ze[1][t];
                Ze[1][t] = 0;

                Zf[0][t] = Zf[1][t];
                Zf[1][t] = 0;

                Zm[1][t] = Zm[0][t];
                Zm[0][t] = 0;
            }
            Zf[0][len0] = 1;
        }
    } // end of for i

#ifdef DEBUG
    putLog("Probalign::reverse matrix:");
    showPassMatrix(Zm, Zf, Ze, len1+1, len0+1, *sequences[0], *sequences[1]);
#endif

    // --- finalize ---
    if (REVPART_FULL_MEMORY)
    {
        for(int i=0; i<=len1; i++)
        {
            delete(Zm[i]);
            delete(Zf[i]);
            delete(Ze[i]);
        }
    } else {
        delete(Zf[0]);
        delete(Ze[0]);
        delete(Zm[0]);

        delete(Zm[1]);
        delete(Zf[1]);
        delete(Ze[1]);
    }

    for(int i=0; i<=len1; i++) delete(Zfm[i]);
    delete(Zf);
    delete(Ze);
    delete(Zm);
    delete(Zfm);

    posterior[0] = 0;

    return (posteriorPtr);
}

/*
 * row-wise calculation of Z matrix
 */
long double** Probalign::doForward(const SafeVector< std::string* > & sequences, const double termgapopen, const double termgapextend, const double d, const double e) {
    const int len0 = (*sequences[0]).size()-1; // without the prefix
    const int len1 = (*sequences[1]).size()-1;
    long double **Zm = NULL, **Zf = NULL, **Ze = NULL, zz = 0;
    double endgapopen, endgapextend;

    // --- initialize ---
    endgapopen = termgapopen;
    endgapextend = termgapextend;

    if(PART_FULL_MEMORY)
    {
        Zf = new long double *[len1 + 1];
        Ze = new long double *[len1 + 1];
        Zm = new long double *[len1 + 1];

        for(int i=0; i<=len1; i++)
        {
            Zf[i] = new long double[len0 + 1];
            Ze[i] = new long double[len0 + 1];
            Zm[i] = new long double[len0 + 1];
        }
    } else {
        Zm = new long double *[len1 + 1];
        Ze = new long double *[2];
        Zf = new long double *[2];
        for (int i=0; i<=len1; i++)
        {
            Zm[i] = new long double[len0 + 1];
        }
        Ze[0] = new long double[len0 + 1];
        Zf[0] = new long double[len0 + 1];
        Ze[1] = new long double[len0 + 1];
        Zf[1] = new long double[len0 + 1];
    }

    if (PART_FULL_MEMORY)
    {
        for(int i=0; i<=len1; i++)
        {
            for(int j=0; j<=len0; j++)
            {
                Zm[i][j] = 0.00;
                Zf[i][j] = 0.00;
                Ze[i][j] = 0.00;
            }
        }
    } else {
        for(int i=0; i<=len1; i++)
        {
            for(int j=0; j<=len0; j++)
            {
                Zm[i][j] = 0;
            }
        }
        for (int j=0; j<=len0; j++)
        {
            Zf[0][j] = 0;
            Ze[0][j] = 0;
            Zf[1][j] = 0;
            Ze[1][j] = 0;
        }
    }

    if (endgaps == false)
    {
        Zm[0][0] = 1.00;

        Zf[0][0] = Ze[0][0] = 0;
        Zf[1][0] = Zm[0][0] * d; // gap open
        Ze[0][1] = Zm[0][0] * d; // gap open
        // >=2ND ROW INIT
        if(PART_FULL_MEMORY)
        {
            for(int i=2; i<=len1; i++)
            {
                Zf[i][0] = Zf[i - 1][0] * e; // gap ext
            }
        }

        // >=2ND COL INIT
        for (int j=2; j<=len0; j++)
        {
            Ze[0][j] = Ze[0][j - 1] * e; // gap ext
        }
    } else {
        Zm[0][0] = 1.00;
        Zf[0][0] = Ze[0][0] = 0;
        Zf[1][0] = Zm[0][0] * endgapopen;
        Ze[0][1] = Zm[0][0] * endgapopen;

        // >=2ND ROW INIT
        if(PART_FULL_MEMORY)
        {
            for(int i=2; i<=len1; i++)
            {
                Zf[i][0] = Zf[i - 1][0] * endgapextend;
            }
        }

        // >=2ND COL INIT
        for(int j=2; j<=len0; j++)
        {
            Ze[0][j] = Ze[0][j - 1] * endgapextend;
        }
    }

    // --- main ---
    for(int i=1; i<=len1; i++)
    {

        for (int j=1; j<=len0; j++)
        {
            const ScoreMatrix::EncodedResidue Si = (*sequences[1])[i];
            const ScoreMatrix::EncodedResidue Tj = (*sequences[0])[j];
            const double score = sub_matrix(Si,Tj);

            double open0, extend0, open1, extend1;

            open0 = open1 = d;
            extend0 = extend1 = e;

            if(endgaps)
            {
                if(i == len1)
                {
                    open0 = endgapopen;
                    extend0 = endgapextend;

                }
                if (j == len0)
                {
                    open1 = endgapopen;
                    extend1 = endgapextend;
                }
            }

            // オリジナルのここにはコメントあり 650行目

            if(PART_FULL_MEMORY)
            {
                Ze[i][j] = Zm[i][j - 1] * open0 + Ze[i][j - 1] * extend0;
                assert(Ze[i][j] < PROBALIGN_HUGE_VALL && "ERROR: huge val error for Ze");

                Zf[i][j] = Zm[i - 1][j] * open1 + Zf[i - 1][j] * extend1;
                assert(Zf[i][j] < PROBALIGN_HUGE_VALL && "ERROR: huge val error for Zf");

                Zm[i][j] = (Zm[i - 1][j - 1] + Ze[i - 1][j - 1] + Zf[i - 1][j - 1]) * score;
                assert(Zm[i][j] < PROBALIGN_HUGE_VALL && "ERROR: huge val error for Zm");

                zz = Zm[i][j] + Ze[i][j] + Zf[i][j];
            } else {
                Ze[1][j] = Zm[i][j - 1] * open0 + Ze[1][j - 1] * extend0;
                assert(Ze[1][j] < PROBALIGN_HUGE_VALL && "ERROR: huge val error for zE");

                Zf[1][j] = Zm[i - 1][j] * open1 + Zf[0][j] * extend1;
                assert(Zf[1][j] < PROBALIGN_HUGE_VALL && "ERROR: huge val error for zF");

                // TODO: ここはijでいいの？1jでは？
                Zm[i][j] = (Zm[i - 1][j - 1] + Ze[0][j - 1] + Zf[0][j - 1]) * score;
                assert(Zm[i][j] < PROBALIGN_HUGE_VALL && "ERROR: huge val error for zM");

                zz = Zm[i][j] + Ze[1][j] + Zf[1][j];
            }
        } // end of for j

        if(PART_FULL_MEMORY == false)
        {
            for (int t=0; t<=len0; t++)
            {
                Ze[0][t] = Ze[1][t];
                Ze[1][t] = 0;

                Zf[0][t] = Zf[1][t];
                Zf[1][t] = 0;
            }

            Zf[1][0] = 1;
        }
    } // end of for i

    // store the sum of zm zf ze (m,n)s in zm's 0,0 th position
    Zm[0][0] = zz;

#ifdef DEBUG
    putLog("Probalign::forward matrix:");
    showPassMatrix(Zm, Zf, Ze, len1+1, len0+1, *sequences[0], *sequences[1]);
#endif

    // --- finalize ---
    if(PART_FULL_MEMORY)
    {
        for(int i=0; i<=len1; i++)
        {
            delete(Zf[i]);
            delete(Ze[i]);
        }
    } else {
        delete(Zf[0]);
        delete(Ze[0]);
        delete(Zf[1]);
        delete(Ze[1]);
    }

    delete(Zf);
    delete(Ze);

    return Zm;
}

// corresponding to Probalign::PostProbs.cc::ComputePostProbs
VPROFILE* Probalign::genProfile(std::string* seq1, std::string* seq2) {
    // --- initialize ---
    SafeVector< std::string* > sequences;
    sequences.push_back(seq1);
    sequences.push_back(seq2);

    double gap_open = argument.gapopen;
    double gap_ext = argument.gapext;
    double beta = argument.beta;

    double termgapopen = exp(beta * 0.0);
    double termgapextend = exp(beta * 0.0);
#ifdef DEBUG
    putLog("probalign beta = " + DS(beta));
    putLog("probalign pre gap open = " + DS(gap_open));
    putLog("probalign pre gap ext  = " + DS(gap_ext));
#endif
    gap_open = exp(beta * gap_open);
    gap_ext = exp(beta * gap_ext);
#ifdef DEBUG
    putLog("probalign gap open = " + DS(gap_open));
    putLog("probalign gap ext  = " + DS(gap_ext));
#endif
    // --- main ---
    long double **MAT1;
    MAT1 = doForward(sequences, termgapopen, termgapextend, gap_open, gap_ext);
    return doReverse(sequences, termgapopen, termgapextend, MAT1, gap_open, gap_ext);
}

void Probalign::showPassMatrix(long double **Zm, long double **Zf, long double **Ze, int numRow, int numCol, std::string& A, std::string& B)
{
    const double ZERO = 0.0;

    for(int r=0; r<numRow; r++)
    {
        stringstream line;

        for(int c=0; c<numCol; c++)
        {
            line << "(";
            for(int i=0; i<3; i++)
            {
                long double m;
                switch(i){
                case 0: m = Zm[r][c]; break;
                case 1: m = Zf[r][c]; break;
                case 2: m = Ze[r][c]; break;
                }
                if(abs(ZERO - m) < 1.0e-15)
                {
                    line << left << setw(10) << "---";
                }else{
                    line << scientific << left << setw(10) << setprecision(2) << m;
                }
            }
            line << ") ";
        }
        putLog(line.str());
    }
    putLog("");

}

