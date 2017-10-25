#pragma once
#include <string>
//#include <vector>
#include <sstream>
#include <cmath>
#include <map>

#ifdef HAVE_SYS_TR1_MEMORY_H
#include <tr1/memory>
#else
#include <memory>
#endif

#include "SafeVector.h"
#include "SparseMatrix.h"

#define EPS 1.0e-15

typedef std::vector<double> vec_double;
typedef std::vector< vec_double > vmat_double;
typedef std::vector<std::string> vec_string;
typedef std::vector< vec_string > vmat_string;
typedef std::vector<int> vec_int;
typedef std::vector< vec_int > vmat_int;

#ifdef HAVE_SYS_TR1_MEMORY_H
typedef std::tr1::shared_ptr<std::string> pString;
#else
typedef std::shared_ptr<std::string> pString;
#endif
typedef std::vector<pString> vec_pString;

#ifdef HAVE_DECL_ISINF
#define Isinf(x) (isinf(x))
#else
#define Isinf(x) (!_finite(x))
#endif

typedef SafeVector<std::string> VS;
typedef SafeVector<VS> VVS;
typedef SparseMatrix* pSparseMatrix;
typedef std::vector<pSparseMatrix> VpSM;
typedef std::vector<VpSM> VVpSM;
typedef SafeVector<char> VC;
typedef SafeVector<VC> VVC;

extern VS fileSearchDirList;

int getIndex(const vec_pString& v, const std::string& s);
void strReplace(std::string& str, const std::string& from, const std::string& to);
void adjustDirDelim(std::string& str);
void initFileSearchDir();
std::string fixFilePath(const std::string& str);
void putError(std::string& str, bool finish = true);
void putError(const std::string& str, bool finish = true);
void putLog(std::string& str);
void putLog(const std::string& str);
void putWarning(std::string& str);
void putWarning(const std::string& str);
void mySplit(const std::string& str, const std::string& delim, std::vector<std::string>& parts);
bool getMultiFASTA(vmat_string& fasta, const std::string& filename, bool toUpperCase, bool omitLowerCase = false);
bool getMultiFASTA(vec_string& names, vec_string& seqs, const std::string& filename, bool toUpperCase, bool omitLowerCase = false);

template<typename T>
T string2binary(const std::string& text) {
    std::stringstream is(text);
    T value;
    is >> value;
    return value;
}

template<typename X>
std::string binary2string(X value) {
    std::stringstream os;
    os << value;
    return os.str();
}

template<typename T>
T SB(const std::string& text) {
    std::stringstream is(text);
    T value;
    is >> value;
    return value;
}

template<typename X>
std::string BS(X value) {
    std::stringstream os;
    os << value;
    return os.str();
}

inline std::string IS(int value) { return BS<int>(value); }
inline std::string DS(double value) { return BS<double>(value); }
inline int SI(const std::string& value) { return SB<int>(value); }
inline double SD(const std::string& value) { return SB<double>(value); }

inline double Log2(double n) { return log(n)/log(2.0); }
inline double Log10(double n) { return log(n)/log(10.0); }
inline double Round(double x) { if(x>0.0){return floor(x+0.5);}else{return -1.0*floor(fabs(x)+0.5);} }

template<typename T>
std::pair<T,T> myMinmax(T i1, T i2) {
    if(i1 <= i2)
    {
        return make_pair(i1,i2);
    }else{
        return make_pair(i2,i1);
    }
}

template<typename T>
void init4Vec(SafeVector< SafeVector< SafeVector< SafeVector<T> >  > >& v, int n) {
    v.resize(n);
    for(int i1=0;i1<n; ++i1)
    {
        v.at(i1).resize(n);
        for(int i2=0;i2<n; ++i2)
        {
            v.at(i1).at(i2).resize(n);
            for(int i3=0;i3<n; ++i3)
            {
                v.at(i1).at(i2).at(i3).resize(n, 0);
            }
        }
    }
}

template<typename T>
void init2Vec(SafeVector< SafeVector<T> >& v, int n) {
    v.resize(n);
    for(int i1=0;i1<n; ++i1)
    {
        v.at(i1).resize(n, 0);
    }
}

template<typename T>
void init1Vec(SafeVector<T>& v, int n) {
    v.resize(n, 0);
}

template<typename T, typename C>
inline pair<T,C> choiceMax(T d, T l , T u, C D, C L, C U) {
    if(d>=l)
    {
        if(d>=u)
        {	// d > u,l
            return make_pair<T,C>(d, D);
        }else{
            // u > d > l
            return make_pair<T,C>(u, U);
        }
    }else{
        if(l>=u)
        {	// l > d,u
            return make_pair<T,C>(l, L);
        }else{
            // u > l > d
            return make_pair<T,C>(u, U);
        }
    }

}

template<typename T, typename C>
inline pair<T,C> choiceMin(T d, T l , T u, C D, C L, C U) {
    if(d<l)
    {
        if(d<u)
        {	// d < u,l
            return make_pair<T,C>(d, D);
        }else{
            // u < d < l
            return make_pair<T,C>(u, U);
        }
    }else{
        if(l<u)
        {	// l < d,u
            return make_pair<T,C>(l, L);
        }else{
            // u < l < d
            return make_pair<T,C>(u, U);
        }
    }

}

template<typename T, typename C>
inline pair<T,C> choiceMin(T d, T l, C D, C L) {
    if(d<l)
    {
        return make_pair<T,C>(d,D);
    }else{
        return make_pair<T,C>(l,L);
    }
}

template<typename T>
inline const T& choiceMin(const T& a, const T& b, const T& c)
{
    //           (a<b)      a<b,c  c<a<b
    return (a < b) ? ((a < c) ? a    : c)
                   : ((b < c) ? b    : c);
    //           (b<a)      b<a,c  c<b<a
}

template<typename T>
inline const T& choiceMin(const T& a, const T& b)
{
    return (a < b) ? a : b;
}

template<typename T>
inline const T& choiceMax(const T& a, const T& b, const T& c)
{
    //           (a>b)      a>b,c  c>a>b
    return (a > b) ? ((a > c) ? a    : c)
                   : ((b > c) ? b    : c);
    //           (b>a)      b>a,c  c>b>a
}

template<typename T>
inline const T& choiceMax(const T& a, const T& b)
{
    return (a > b) ? a : b;
}

std::string* genOmittedString(const std::string& s, const std::vector<bool>& flag);


