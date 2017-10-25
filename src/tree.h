#pragma once

#include "SafeVector.h"
#include "utility.h"
#include "SparseMatrix.h" // included for VPROFILE,VVPROFILE


// 木構造における葉，節
// pLeft = pRight = NULL: 葉, pLeft != NULL && pRight != NULL: 節, 片方だけNULLは許さない
// id: 葉のID．葉の名前に相当．木構造内でユニークである必要がある．節のIDは-1
class tr_node{
private:
public:
    tr_node* pLeft;
    tr_node* pRight;
    int id;	//この節の名前ID
    int numLeafs; // この節に葉が何個属しているか

    // variables set by the create tree method
    double distToParent; // 親までの距離

    // 葉の重み ClustalW風
    double weight;

    tr_node(); // 葉か節かは未定
    tr_node(tr_node* pL, tr_node* pR); // 節を作成
    tr_node(int i);	// 葉を作成 i>=0
    ~tr_node();
    void free();
    void getIDList(VI& vec_ID); // 自分が節の時，左サブツリー -> 右サブツリー の順で葉のIDを列挙する．自分が葉の時は自分のIDを返す．
    std::string toString(); // 節の時，(左サブツリー,右サブツリー)を出力　葉の時，idを出力
    std::string toString(const VS& name); //name[id]を用いて木構造を出力
    void getRPN(VI& vec_ID, const int DEL);	//逆ポーランド記法でIDリストを返す．節の時，左->右->DELと出力，葉の時idを出力
    void setClustalWeight(VPROFILE& weights, PROFILE_TYPE parentWeight = PROFILE_ZERO);
};

// legacy code
PROFILE_TYPE getUPGMADistance(const tr_node& tree, const VVPROFILE& distMatrix);
tr_node* do_upgma(SafeVector<tr_node*>& node_list, const VVPROFILE& distMatrix);

tr_node* genTreeUPGMA(const VVPROFILE& distMatrix); // 新しいUPGMAのコード．こちらのほうが速くてエレガント

