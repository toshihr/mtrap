#include <iostream>
#include "entropy.h"
#include "utility.h"

using namespace std;

bool Entropy::calcComplexProb(Entropy::VVT& complexProb, const std::string& s1, const std::string& s2){
    //サイズの初期化
    complexProb.clear();
    complexProb.resize(PROB_TABLE_SIZE);
    for(int i=0; i<PROB_TABLE_SIZE; i++)complexProb.at(i).resize(PROB_TABLE_SIZE, ZERO);

    int l = (s1.size() < s2.size())?(int)s1.size():(int)s2.size();

    if(s1.size() != s2.size())cerr << "calcComplexProb: wrong length." << endl;

    //カウント
    for(int i=0; i<l; i++)complexProb.at(s1.at(i)).at(s2.at(i)) += ONE;
    //全部の要素を要素数で割って確率に変換する
    for(int i=0; i<PROB_TABLE_SIZE; i++)for(int j=0; j<PROB_TABLE_SIZE; j++)complexProb.at(i).at(j) /= l;

    return true;
}

bool Entropy::calcProb(Entropy::VT& prob, const std::string s){
    //頻度表の初期化
    prob.clear();
    for(int i=0; i<PROB_TABLE_SIZE; i++)prob.push_back(ZERO);

    //頻度のカウント
    for(int i=0; i<(int)s.size(); i++)prob.at(s.at(i)) += ONE;

    //確率計算
    for (int i=0; i<PROB_TABLE_SIZE; i++)prob.at(i) /= static_cast<T>(s.size());

    return true;
}

Entropy::T Entropy::getEntropy(const Entropy::VT& prob){
    T res = ZERO;

    //エントロピーの計算
    for(VT::const_iterator i = prob.begin(); i!=prob.end(); ++i)
    {
        if(*i == ZERO)continue;
        res -= *i * Log2(*i);
    }

    //誤差による負の値の対処
    return max(res, ZERO);
}

Entropy::T Entropy::getMutualEntropy(const Entropy::VT& prob1, const Entropy::VT& prob2, const Entropy::VVT& complexProb){
    T res = 0.0;

    //エントロピーの計算
    for (size_t i=0; i<PROB_TABLE_SIZE; i++)for (size_t j=0; j<PROB_TABLE_SIZE; j++)
    {
        if(prob1.at(i) == ZERO || prob2.at(j) == ZERO || complexProb.at(i).at(j) == ZERO)continue;
        res += complexProb.at(i).at(j) * Log2(complexProb.at(i).at(j)/prob1.at(i)/prob2.at(j));
    }

    //誤差による負の値の対処
    return max(res, ZERO);
}

Entropy::T Entropy::getEER(const Entropy::VT& prob1, const Entropy::VT& prob2, const Entropy::VVT& complexProb){
    T I = getMutualEntropy(prob1, prob2, complexProb);
    const T e1 = getEntropy(prob1);
    const T e2 = getEntropy(prob2);
    const T ret = ONE - (I/e1 + I/e2)/2;

    if(max(max(max(e1,e2),I), ZERO) > ZERO)
    {
        return max(ret, ZERO);
    }else{
        //0除算によるnanが発生するとき
        return ZERO;
    }
}

Entropy::T Entropy::getEER2(const Entropy::VT& prob1, const Entropy::VT& prob2, const Entropy::VVT& complexProb){
    const T I = getMutualEntropy(prob1, prob2, complexProb);
    const T e1 = getEntropy(prob1);
    const T e2 = getEntropy(prob2);
    const T ret = ONE - I/(e1 + e2 - I);

    if(max(max(max(e1,e2),I), ZERO) > ZERO)
    {
        return max(ret, ZERO);
    }else{
        //0除算によるnanが発生するとき
        return ZERO;
    }
}

Entropy::T Entropy::getECD(const Entropy::VT& prob1, const Entropy::VT& prob2, const Entropy::VVT& complexProb){
    T I = getMutualEntropy(prob1, prob2, complexProb);

    return getEntropy(prob1) - I;
}
