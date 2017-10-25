#ifndef _DISTANCE_H_
#define _DISTANCE_H_

#include <string>
#include "entropy.h"
#include "scorematrix.h"
#include "utility.h"

namespace Distance {

extern bool USE_ENCODED_GAP;
extern char GAP;
extern ScoreMatrix::EncodedResidue EGAP;

/*
 * GID^    : Global sequence identity
 * LID^    : Local sequence identity
 * PDIST   : Proportion of residue sites ( = 1 - LID )
 * EER     : Entropy evolution rate
 * EER2    : Entropy evolution rate2 (satisfy triangle inequality)
 * ECD*    : Entropic chaos degree
 * ECDR    : = (ECDR(A,B)+ECDR(B,A))/2
 * ECDR2   : = (ECD(A,B)/S(B)+ECD(B,A)/S(A))/2
 * MUTUAL* : Mutual information
 *
 *  ^: this is not a distance i.e. the value 1 means most similar.
 *  *: is not normalized
 */
typedef enum { GID, LID, PDIST, EER, EER2, ECD, ECDR, ECDR2, MUTUAL } MODE;

template<typename T>
T calcDistance(const string& s1, const string& s2, MODE mode)
{
    assert(s1.size() == s2.size() && "size is wrong!!");

    if( mode < EER )
    {
        switch(mode)
        {
        case GID:
        {
            const int len = static_cast<int>(s1.size());
            int same = 0;
            string::const_iterator i1 = s1.begin();
            string::const_iterator i2 = s2.begin();
            for(;i1 != s1.end(); ++i1, ++i2)
            {
                if(*i1 == *i2) same++;
            }
            return static_cast<T>(same) / len;
            break;
        }
        case LID:
        {
            const char G = ((USE_ENCODED_GAP) ? EGAP : GAP);
            int len = 0;
            int same = 0;
            string::const_iterator i1 = s1.begin();
            string::const_iterator i2 = s2.begin();
            for(;i1 != s1.end(); ++i1, ++i2)
            {
                if(*i1 == G || *i2 == G) continue;
                if(*i1 == *i2) same++;
                len++;
            }
            return static_cast<T>(same) / len;
            break;
        }
        case PDIST:
        {
            const char G = ((USE_ENCODED_GAP) ? EGAP : GAP);
            int len = 0;
            int same = 0;
            string::const_iterator i1 = s1.begin();
            string::const_iterator i2 = s2.begin();
            for(;i1 != s1.end(); ++i1, ++i2)
            {
                if(*i1 == G || *i2 == G) continue;
                if(*i1 == *i2) same++;
                len++;
            }
            return 1 - static_cast<T>(same) / len;
            break;
        }
        default:
        {
            putError("Not implemented.");
            break;
        }
        }
    }else{
        // maximum entropy
        const T P = static_cast<T>(1) / 21;
        const T maxS = -21 * P * Log2(P);

        Entropy::VT prob1, prob2;
        Entropy::VVT complexProb;

        Entropy::calcProb(prob1, s1);
        Entropy::calcProb(prob2, s2);
        Entropy::calcComplexProb(complexProb, s1, s2);

        switch(mode)
        {
        case EER:
        {
            return static_cast<T>(Entropy::getEER(prob1, prob2, complexProb));
            break;
        }
        case EER2:
        {
            return static_cast<T>(Entropy::getEER2(prob1, prob2, complexProb));
            break;
        }
        case ECD:
        {
            return static_cast<T>(Entropy::getECD(prob1, prob2, complexProb));
            break;
        }
        case ECDR:
        {
            const Entropy::T v1 = Entropy::getECD(prob1, prob2, complexProb);
            const Entropy::T v2 = Entropy::getECD(prob2, prob1, complexProb);
            return static_cast<T>((v1+v2) / maxS / 2);
            break;
        }
        default:
        {
            putError("Not implemented.");
            break;
        }
        }
    }

    putError("Not implemented.");
    return 0;
}


}

#endif

