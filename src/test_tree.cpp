#include "SafeVector.h"
#include "SparseMatrix.h"
#include "tree.h"

int main(){
    VS inStr;
    //test
    //(((a,b),(c,d)),e)
    inStr.push_back("adasfasadsfa");
    inStr.push_back("badsfasdfasd");
    inStr.push_back("casdfadsfafa");
    inStr.push_back("dadfasdfasdf");
    inStr.push_back("easdfasfasas");

    tr_node* pNode_a = new tr_node(0);
    tr_node* pNode_b = new tr_node(1);
    tr_node* pNode_c = new tr_node(2);
    tr_node* pNode_d = new tr_node(3);
    tr_node* pNode_e = new tr_node(4);


//	tr_node* pNode_ab = new tr_node(pNode_a, pNode_b);
//	tr_node* pNode_ab_c = new tr_node(pNode_ab, pNode_c);
//	tr_node* pNode_de = new tr_node(pNode_d, pNode_e);
//	tr_node* pRoot = new tr_node(pNode_ab_c, pNode_de);


    VVPROFILE distMatrix(5,VPROFILE(5, PROFILE_ZERO));
    VVPROFILE distMatrixSame(5, VPROFILE(5, PROFILE_ZERO));
    distMatrix.at(0).at(0) = 0;
    distMatrix.at(0).at(1) = 3;
    distMatrix.at(0).at(2) = 5;
    distMatrix.at(0).at(3) = 7;
    distMatrix.at(0).at(4) = 8;
    distMatrix.at(1).at(0) = 3;
    distMatrix.at(1).at(1) = 0;
    distMatrix.at(1).at(2) = 6;
    distMatrix.at(1).at(3) = 8;
    distMatrix.at(1).at(4) = 8;
    distMatrix.at(2).at(0) = 5;
    distMatrix.at(2).at(1) = 6;
    distMatrix.at(2).at(2) = 0;
    distMatrix.at(2).at(3) = 9;
    distMatrix.at(2).at(4) = 10;
    distMatrix.at(3).at(0) = 7;
    distMatrix.at(3).at(1) = 8;
    distMatrix.at(3).at(2) = 9;
    distMatrix.at(3).at(3) = 0;
    distMatrix.at(3).at(4) = 7;
    distMatrix.at(4).at(0) = 8;
    distMatrix.at(4).at(1) = 8;
    distMatrix.at(4).at(2) = 10;
    distMatrix.at(4).at(3) = 7;
    distMatrix.at(4).at(4) = 0;
    for(size_t i=0;i<distMatrix.size();i++)
    {
        for(size_t j=0;j<distMatrix.size();j++)
        {
            distMatrixSame[i][j] = distMatrix[i][j];
        }
    }

    SafeVector<tr_node*> node_list;
    node_list.push_back(pNode_a);
    node_list.push_back(pNode_b);
    node_list.push_back(pNode_c);
    node_list.push_back(pNode_d);
    node_list.push_back(pNode_e);

    tr_node* pRoot = do_upgma(node_list, distMatrix);
    cout << "old UPGMA tree:" << pRoot->toString() << endl;
    tr_node* pRoot2 = genTreeUPGMA(distMatrixSame);
    cout << "new UPGMA tree:" << pRoot2->toString() << endl;
    VPROFILE seqWeights(distMatrix.size());
    pRoot2->setClustalWeight(seqWeights);
    for(int i=0; i<(int)seqWeights.size(); i++)
    {
        cout << seqWeights[i] << ",";
    }
    cout << endl;

    VI rpn;
    pRoot->getRPN(rpn, -1);
    for(int i=0; i<(int)rpn.size(); i++)
    {
        cout << rpn.at(i) << ",";
    }
    cout << endl;
    cout << getUPGMADistance(*pRoot, distMatrix) << endl;

    pRoot->free();

    delete pRoot;
    delete pRoot2;

    cout << "finished." << endl;

    return 0;
}

