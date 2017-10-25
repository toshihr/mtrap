#include <iostream>
#include <map>
#include <algorithm>
#include <cassert>
#include <cfloat>
#include <string>
#include <sstream>
#include "SafeVector.h"
#include "utility.h"
#include "tree.h"
#include "SparseMatrix.h" // included for VPROFILE,VVPROFILE

using namespace std;

tr_node::tr_node()
    : pLeft(NULL), pRight(NULL), id(-1), numLeafs(0), distToParent(0.0), weight(0.0)
{
}

tr_node::tr_node(tr_node* pL, tr_node* pR)
    : pLeft(pL), pRight(pR), id(-1), numLeafs(pL->numLeafs + pR->numLeafs), distToParent(0.0), weight(0.0)
{
}

tr_node::tr_node(int i)
    : pLeft(NULL), pRight(NULL), id(i), numLeafs(1)
{
    assert(id>=0 && "leafs must have id >= 0");
}

tr_node::~tr_node(){
}

void tr_node::free(){
    if(pLeft != NULL)
    {
        delete pLeft;
        pLeft = NULL;
    }

    if(pRight != NULL)
    {
        delete pRight;
        pRight = NULL;
    }
}

void tr_node::getIDList(VI& vec_ID){
    if(pLeft != NULL)pLeft->getIDList(vec_ID);
    if(pRight != NULL)pRight->getIDList(vec_ID);
    if(pLeft == NULL && pRight == NULL)vec_ID.push_back(id);
}

string tr_node::toString(){
//	if(pLeft != NULL && pRight != NULL)
    if(id < 0)
    {
        return "(" + pLeft->toString() + "," + pRight->toString() + ")";
    } else {
        return binary2string(id);
    }
}

string tr_node::toString(const VS& name){
    if(pLeft != NULL && pRight != NULL)
    {
        return "(" + pLeft->toString() + "," + pRight->toString() + ")";
    } else {
        return name.at(id);
    }
}


// pLeftにあるIDListとpRightにあるIDListを用いて
// グループ間距離を計算する
PROFILE_TYPE getUPGMADistance(const tr_node& tree, const VVPROFILE& distMatrix){
    assert(tree.pLeft);
    assert(tree.pRight);

    VI id_left;
    VI id_right;
    tree.pLeft->getIDList(id_left);
    tree.pRight->getIDList(id_right);

    /*
    cout << "Left IDs = ";
    for(int i=0; i<(int)id_left.size(); i++)cout << id_left.at(i) << ",";
    cout << endl;

    cout << "Right IDs = ";
    for(int i=0; i<(int)id_right.size(); i++)cout << id_right.at(i) << ",";
    cout << endl;
    */

    //UPGMA法で距離を計算
    PROFILE_TYPE dist = PROFILE_ZERO;
    for(int i=0; i<(int)id_left.size(); i++)
    {
        for(int j=0; j<(int)id_right.size(); j++)
        {
            int d1 = id_left.at(i);
            int d2 = id_right.at(j);
            dist += distMatrix.at(d1).at(d2);
            //cout << "(" << id_left.at(i) << "," << id_right.at(j) << ")";
        }
    }
    dist /= (int)(id_left.size() * id_right.size());
    //cout << "/" << (int)(id_left.size() * id_right.size()) << endl;

    return dist;
}

//ノードのリストを１段階だけまとめる
//これを繰り返して最終的に１つにまとまると完成
//距離値が同じになるものが３つ以上ある場合、後ろの組が優先される
// 1. Ohya Lab.流: すべての葉間の距離の平均をだす．毎回すべての葉までの平均値をだすため遅い
// d([[a,b],c],d) = (d(a,d)+d(b,d)+d(c,d)) / (3*1), すべての葉までの平均を計算する
// 毎度，葉のリストを構築し，平均を求めている．ここでn!程度の計算量がかかる
// 2. 世間流
// d([[a,b],c],d) = (d([a,b],d)*2+d(c,d)*1) / (2+1), ここでd([a,b],d)はすでに計算してあるため再計算は必要ない
// 毎度距離行列を全部作りなおすとn(n-1)/2の階乗つまりn!(n-1)!/2程度計算コストがかかることになりとんでもないことになる
//
// node_list: まとめたい葉及び節の集合．この関数をよびだすごとにこの集合は書き換えられる
// distMatrix: 距離行列(対称行列)
tr_node* do_upgma(SafeVector<tr_node*>& node_list, const VVPROFILE& distMatrix){
    //全ての要素間の距離を計算し、最小値を取る組も算出
    const int nsize = (int)node_list.size();
    PROFILE_TYPE min;
    int min_i, min_j;
    bool minFlag = false;
    VVPROFILE distance(nsize, VPROFILE(nsize, PROFILE_ZERO));

    for(int i=0; i<nsize-1; i++)
    {
        for(int j=i+1; j<nsize; j++)
        {
            const tr_node node(node_list.at(i), node_list.at(j));
            const PROFILE_TYPE d = getUPGMADistance(node, distMatrix);
            distance.at(i).at(j) = d;
            if(minFlag == false || (minFlag == true && d < min))
            {
                min = d;
                min_i = i;
                min_j = j;
                minFlag = true;
            }
        }
    }

    //最小値を取る組をグループにまとめる
    tr_node* pNode1 = node_list.at(min_i);
    tr_node* pNode2 = node_list.at(min_j);
    tr_node* pNewNode = new tr_node(pNode1, pNode2);
    node_list.erase(remove(node_list.begin(), node_list.end(), pNode1));
    node_list.erase(remove(node_list.begin(), node_list.end(), pNode2));
    node_list.push_back(pNewNode);

    //グループ化を再帰的に実行
    if((int)node_list.size() >= 2)return do_upgma(node_list, distMatrix);

    return *node_list.begin();
}

/*!
* 逆ポーランド記法でIDリストを返す
* 帰りがけ順
* @param[out] vec_ID 逆ポーランド記法でのIDリスト
* @param[in] DEL デリミターで用いる値
*
*/
void tr_node::getRPN(VI& vec_ID, const int DEL){
    if(pLeft != NULL && pRight != NULL)	// node
    {
        pLeft->getRPN(vec_ID, DEL);
        pRight->getRPN(vec_ID, DEL);
        vec_ID.push_back(DEL);
    }else{	// leaf
        vec_ID.push_back(id);
    }
}

/* ClustalW風Sequence Weightを設定
* root node->setWeight()を呼び出すことですべての節，葉のWeightが設定される
* 最終的に葉のweightは　枝長 / 属している葉の数　をルートから足しあわせたもの　となる
* weights: 葉のweightを返すコンテナ(root node->numLeafs数で初期化されている必要あり)
*/
void tr_node::setClustalWeight(VPROFILE& weights, PROFILE_TYPE parentWeight){
/*
#ifdef DEBUG
    stringstream ss;
    ss << this->id << ":parentWeight=" << parentWeight << ",distToParent/numLeafs=" << this->distToParent << "/" << this->numLeafs;
    putLog(ss.str());
#endif
*/
    this->weight = parentWeight + this->distToParent / this->numLeafs;
    if(pLeft != NULL) pLeft->setClustalWeight(weights, this->weight);
    if(pRight != NULL) pRight->setClustalWeight(weights, this->weight);
    if(id >= 0)weights[id] = this->weight;
    //	cout << this->toString() << " = " << this->weight << endl;
}

/* === UPGMA method ===
* 入力：
* distMatrix: 葉の間の距離行列(対象行列かつdistMatrix[i][j]>=0，distMatrix[i][j] i<jの値を用いる)
* 出力：
* ルート節へのポインタ．葉は距離行列におけるインデックスで構成される．
* アルゴリズム：
* 世間流で行い，距離行列の再計算時，新たなグループは距離行列のインデックスの小さい方を利用する
* 参考：probcons::EvolutionaryTree.h
*/
tr_node* genTreeUPGMA(const VVPROFILE& distMatrix){
    const int nsize = (int)distMatrix.size();
    VVPROFILE distances(nsize, VPROFILE(nsize, PROFILE_ZERO));
    SafeVector<tr_node*> nodes(nsize, NULL);
    SafeVector<bool> valid (nsize, true);

    // make a copy of the distance matrix
    for(int i=0; i<nsize-1; i++)
    {
        for(int j=i+1; j<nsize; j++)
        {
            distances[i][j] = distMatrix[i][j];
            assert(distances[i][j] >= 0 && "distance must be positive value.");
        }
    }

    // create all the leaf nodes
    for(int i=0; i<nsize; i++)
    {
        nodes[i] = new tr_node(i);
        assert(nodes[i]);
        nodes[i]->distToParent = 0.0;
    }

    // repeat until only a single node left
    for(int numLeafsLeft = nsize; numLeafsLeft > 1; numLeafsLeft--)
    {
        PROFILE_TYPE minValue = PROFILE_MAX;
        pair<int,int> minPairAB;

        // find the closest pair
        for(int i=0; i<nsize-1; i++) if(valid[i])
        {
            for(int j=i+1; j<nsize; j++)  if(valid[j])
            {
                if(distances[i][j] < minValue)
                {
                    minValue = distances[i][j];
                    minPairAB = make_pair(i,j);
                }
            }
        }

        const int numLeafsA = nodes[minPairAB.first]->numLeafs;
        const int numLeafsB = nodes[minPairAB.second]->numLeafs;

        // set the branch length from A,B to the parent node (A,B)
        nodes[minPairAB.first]->distToParent = minValue / 2.0;
        nodes[minPairAB.second]->distToParent = minValue / 2.0;

//		cout << "min pair = " << nodes[minPairAB.first]->toString() << "(num.leafs=" << numLeafsA << ")," << nodes[minPairAB.second]->toString() << "(num.leafs=" << numLeafsB << ")" << endl;

        // merge the closest pair
        tr_node *newParent = new tr_node(nodes[minPairAB.first], nodes[minPairAB.second]);
        nodes[minPairAB.first] = newParent;
        nodes[minPairAB.second] = NULL;

        // update the distance d([A,B],i) = ( d(A,i)*n(A) + d(B,i)*n(B) ) / ( n(A) + n(B) )
        for(int i=0; i<nsize; i++) if(i != minPairAB.first && i != minPairAB.second && valid[i])
        {
            const pair<int,int> iA = myMinmax(i,minPairAB.first);
            const pair<int,int> iB = myMinmax(i,minPairAB.second);

            const PROFILE_TYPE d = (distances[iA.first][iA.second] * numLeafsA
                                  + distances[iB.first][iB.second] * numLeafsB)
                                  / (numLeafsA + numLeafsB);
            distances[iA.first][iA.second] = d;

//			cout << "distance between " << newParent->toString() << " and " << nodes[i]->toString() << " = " << d << endl;
        }

        // mark the second node entry as no longer valid
        valid[minPairAB.second] = false;
    }

    assert(nodes[0]);
    return nodes[0];
}

