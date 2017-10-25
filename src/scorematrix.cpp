/* スコア行列に関するクラス
 * スコア: PAM, BLOSUM等
 * PAMって何？ 250?350?得体がしれないのでいつかはずすべき
 * 距離: スコア行列を[0,1]にスケーリングしたもの
 * ギャップスコアに関しては管理クラス等含め調整の必要あり
 * なぜサイズは25?
 * distance版PAMはDISTPAMモードとかとして区別するのもよいかも
 */
#include <iostream>
#include <iomanip>
#include <cassert>
#include <map>
#include <string>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include "scorematrix.h"
#include "tqmatrix.h"
#include "weightedmatrix.h"
#include "utility.h"
#include "config.h"

using namespace std;

void ScoreMatrix::init()
{
    initTQ();
    initTable();
}

//! load & set the ramachandran matrix 将来的にはweightedモードに吸収
/*!

    \param rmode
*/
void ScoreMatrix::setRamachandran(const string& fileName)
{
    loadTriangleFormat(fileName, ScoreMatrix::rmatrix);
    // Average B,Z,X at Ramachandran quantity
    for(int i=0; i<RANK; i++)
    {
        // B,Z,Xの平均化の順番で結果が若干変わるが気にしない
        // M(B,aa) = (M(D,aa) + M(N,aa)) / 2
        if( i != RESIDUE_INDEX['B'] )
        {
            rmatrix(RESIDUE_INDEX['B'], i) = (rmatrix(RESIDUE_INDEX['D'], i) + rmatrix(RESIDUE_INDEX['N'], i)) / 2.0;
            rmatrix(i, RESIDUE_INDEX['B']) = rmatrix(RESIDUE_INDEX['B'], i);
        }
        // M(Z,aa) = (M(E,aa) + M(Q,aa)) / 2
        if( i != RESIDUE_INDEX['Z'] )
        {
            rmatrix(RESIDUE_INDEX['Z'], i) = (rmatrix(RESIDUE_INDEX['E'], i) + rmatrix(RESIDUE_INDEX['Q'], i)) / 2.0;
            rmatrix(i, RESIDUE_INDEX['Z']) = rmatrix(RESIDUE_INDEX['Z'], i);
        }
        // M(X,aa) = ΣM(j,aa) / (RANK-1)   jはX以外
        if( i != RESIDUE_INDEX['X'] )
        {
            double sum=0.0;
            for(int j=0; j<RANK; j++) if(j != RESIDUE_INDEX['X']) sum += rmatrix(j,i);
            rmatrix(RESIDUE_INDEX['X'], i) = sum / (RANK-1);
            rmatrix(i, RESIDUE_INDEX['X']) = rmatrix(RESIDUE_INDEX['X'], i);
        }
    }
}

void ScoreMatrix::initTable()
{
    for(int i=0; i<(256^sizeof(char)); i++)RESIDUE_INDEX[i] = 127;
    for(int i=0; i<RANK+2; i++)RESIDUE_INDEX[RESIDUE[i]] = i;
}


void ScoreMatrix::saveEMBOSSformat(const string& fileName, const AminoMatrix& mat, const double& iscale, const string& database, const double& mutualEntropy)
{
    const int rmax = ScoreMatrix::RANK+1;
    const int rmax_extend = ScoreMatrix::RANK;

    ofstream ofs(fileName.c_str());
    // --- header ---
    ofs << "#  Matrix made by mtrap" << endl;
    ofs << "#" << endl;
    ofs << "#  BLOSUM Clustered Scoring Matrix in 1/" << iscale << " Bit Units" << endl;
    ofs << "#  Database = " << database << endl;
    ofs << "#  Cluster Percentage: >= ??" << endl;
    ofs << "#  Entropy = " << mutualEntropy << endl;
    // --- first row ---
    ofs << " ";
    for(int i1=0; i1<rmax; ++i1)
    {
        if(ScoreMatrix::countModeNoRound == false)
        {
            ofs << "  " << ScoreMatrix::decode(i1);
        }else{
            ofs << "      " << ScoreMatrix::decode(i1);
        }
    }
    ofs << endl;
    // --- main ---
    for(int i1=0; i1<rmax_extend; ++i1)
    {
        ofs << ScoreMatrix::decode(i1);
        for(int i2=0; i2<rmax_extend; ++i2)
        {
            if(ScoreMatrix::countModeNoRound == false)
            {
                ofs << " " << std::setw(2) << std::right << static_cast<int>(Round(mat(i1,i2)));
            }else{
                //ofs.form(" %2.3f",mat(i1,i2));
                ofs << " " << std::setw(6) << setiosflags(ios::fixed) << std::setprecision(3) << mat(i1,i2);
            }
        }
        if(ScoreMatrix::countModeNoRound == false)
        {
            ofs << "  0" << endl;
        }else{
            ofs << "  0.000" << endl;
        }
    }
    // --- last row ---
    ofs << ScoreMatrix::GAP;
    for(int i1=0; i1<rmax; ++i1)
    {
        if(ScoreMatrix::countModeNoRound == false)
        {
            ofs << "  0";
        }else{
            ofs << "  0.000";
        }
    }
    ofs << endl;

    ofs.close();
}

void ScoreMatrix::loadEMBOSSformat(const string& fileName, AminoMatrix& mat)
{
    string AAlist = "@";
    string buff;
    bool existBZX = true;

    ifstream ifs(fileName.c_str());
    if(!ifs) putError(fileName + " is not found.", true);

    mat.setConstant(0);

    // omit the head comment line.
    while(ifs && getline(ifs, buff) && buff.at(0) == '#') ;

    // construct AA list from the first row
    {
        VS list;
        mySplit(buff, ", ", list);
        for(VS::iterator i=list.begin(); i!=list.end(); ++i)
        {
            const unsigned char AA = fixAA(i->at(0));
            if(AA != ScoreMatrix::GAP && AA != 0)
            {
                AAlist.push_back(AA);
            }
        }
    }

    // treat B,Z,X
    if(find(AAlist.begin(),AAlist.end(),'B') == AAlist.end())
    {
        existBZX = false;
    }

    while(ifs && getline(ifs, buff))
    {
        VS list;
        mySplit(buff, ", ", list);
        if(list.size() == 0)break;

        // get a residue from the first column
        const unsigned char AA1 = fixAA(list[0][0]);

        // if the residue is gap then skip this turn
        if(AA1 == ScoreMatrix::GAP || AA1 == 0) continue;

        for(int c=1; c<(int)list.size(); c++)
        {
            const unsigned char AA2 = AAlist[c];
            if(AA2 == ScoreMatrix::GAP || AA2 == 0) continue;

            const double value = SD(list[c]);

            mat(RESIDUE_INDEX[AA1], RESIDUE_INDEX[AA2]) = value;
            mat(RESIDUE_INDEX[AA2], RESIDUE_INDEX[AA1]) = value;
        }
    }

    ifs.close();
}


void ScoreMatrix::saveTriangleFormat(const string& fileName, const AminoMatrix& mat)
{

}


void ScoreMatrix::loadTriangleFormat(const string& fileName, AminoMatrix& mat)
{
#ifdef DEBUG
    cout << "DEBUG: [loadTriangleFormat] fileName=" << fileName << endl;
#endif

    ifstream ifs(fileName.c_str());
    string buff;
    vector<string> AAlist;

    mat.setConstant(ScoreMatrix::NQ);

    // omit the head comment line
    while(ifs && getline(ifs, buff) && buff.at(0) == '#') ;

    // construct value list
    do
    {
        vector<string> list;
        vector<string> vs;
        mySplit(buff, ", ", vs);
        {
            vector<string>::iterator i=vs.begin();
            AAlist.push_back(*i);
            i++;
            for(; i!=vs.end(); ++i)list.push_back(*i);
        }
        const unsigned char AA1 = fixAA(AAlist.rbegin()->at(0));
        if(AA1 == ScoreMatrix::GAP) continue;

        for(int i=0; i<(int)list.size(); i++)
        {
            const unsigned char AA2 = fixAA(AAlist.at(i).at(0));
            if(AA2 == ScoreMatrix::GAP) continue;

            const double value = string2binary<double>(list.at(i));
            mat(RESIDUE_INDEX[AA1], RESIDUE_INDEX[AA2]) = value;
            mat(RESIDUE_INDEX[AA2], RESIDUE_INDEX[AA1]) = value;
        }
    } while(ifs && getline(ifs, buff));

    ifs.close();

    // debug
    assert( !((abs(mat.array() - ScoreMatrix::NQ) < EPS).any()));
}


/*!
    指定された行列をセット
    translateToDistance: 行列およびギャップコストを距離[0,1]に変換するかどうかのフラグ
                         false: Weighted Matrixの計算時あるいは最初から距離距離行列な場合のみ false で利用
    @return 変換に成功したらtrueを返す．
*/
bool ScoreMatrix::loadMatrix(const string& fileName, bool translateToDistance, AminoMatrix& mat){
    ScoreMatrix::MODE m;
    if(fileName == "CID") m = ScoreMatrix::CID;
    else if(fileName == "CPAM") m = ScoreMatrix::CPAM;
    else if(fileName == "CGONNET40") m = ScoreMatrix::CGONNET40;
    else if(fileName == "CGONNET80") m = ScoreMatrix::CGONNET80;
    else if(fileName == "CGONNET120") m = ScoreMatrix::CGONNET120;
    else if(fileName == "CGONNET160") m = ScoreMatrix::CGONNET160;
    else if(fileName == "CGONNET250") m = ScoreMatrix::CGONNET250;
    else if(fileName == "CGONNET300") m = ScoreMatrix::CGONNET300;
    else if(fileName == "CGONNET350") m = ScoreMatrix::CGONNET350;
    else m = ScoreMatrix::ORIGINAL;

    // TODO: Fix!! なぜかコンストラクタで作成したテーブルが壊れるため，ここで再度構築しなおす
    initTable();
    mode = m;
    bool d = false;
    // スコア行列を読み込み
    switch(m){
    case CID:
        mat.setIdentity();
        break;
    case CPAM:
        for(int y=0; y<RANK; y++)
        {
            for(int x=0; x<RANK; x++)
            {
                mat(x,y) = PAM[x][y];
            }
        }
        break;
    case CGONNET40:
        for(int y=0; y<RANK; y++)
        {
            for(int x=0; x<RANK; x++)
            {
                unsigned int index = (1+max(x,y))*max(x,y)/2 + min(x,y);
                unsigned char c1 = clustalOrder[x];
                unsigned char c2 = clustalOrder[y];
                mat(RESIDUE_INDEX[c1], RESIDUE_INDEX[c2]) = ((double)gon40[index]) / 10;
            }
        }
        break;
    case CGONNET80:
        for(int y=0; y<RANK; y++)
        {
            for(int x=0; x<RANK; x++)
            {
                unsigned int index = (1+max(x,y))*max(x,y)/2 + min(x,y);
                unsigned char c1 = clustalOrder[x];
                unsigned char c2 = clustalOrder[y];
                mat(RESIDUE_INDEX[c1], RESIDUE_INDEX[c2]) = ((double)gon80[index]) / 10;
            }
        }
        break;
    case CGONNET120:
        for(int y=0; y<RANK; y++)
        {
            for(int x=0; x<RANK; x++)
            {
                unsigned int index = (1+max(x,y))*max(x,y)/2 + min(x,y);
                unsigned char c1 = clustalOrder[x];
                unsigned char c2 = clustalOrder[y];
                mat(RESIDUE_INDEX[c1], RESIDUE_INDEX[c2]) = ((double)gon120[index]) / 10;
            }
        }
        break;
    case CGONNET160:
        for(int y=0; y<RANK; y++)
        {
            for(int x=0; x<RANK; x++)
            {
                unsigned int index = (1+max(x,y))*max(x,y)/2 + min(x,y);
                unsigned char c1 = clustalOrder[x];
                unsigned char c2 = clustalOrder[y];
                mat(RESIDUE_INDEX[c1], RESIDUE_INDEX[c2]) = ((double)gon160[index]) / 10;
            }
        }
        break;
    case CGONNET250:
        for(int y=0; y<RANK; y++)
        {
            for(int x=0; x<RANK; x++)
            {
                unsigned int index = (1+max(x,y))*max(x,y)/2 + min(x,y);
                unsigned char c1 = clustalOrder[x];
                unsigned char c2 = clustalOrder[y];
                mat(RESIDUE_INDEX[c1], RESIDUE_INDEX[c2]) = ((double)gon250[index]) / 10;
            }
        }
        break;
    case CGONNET300:
        for(int y=0; y<RANK; y++)
        {
            for(int x=0; x<RANK; x++)
            {
                unsigned int index = (1+max(x,y))*max(x,y)/2 + min(x,y);
                unsigned char c1 = clustalOrder[x];
                unsigned char c2 = clustalOrder[y];
                mat(RESIDUE_INDEX[c1], RESIDUE_INDEX[c2]) = ((double)gon120[index]) / 10;
            }
        }
        break;
    case CGONNET350:
        for(int y=0; y<RANK; y++)
        {
            for(int x=0; x<RANK; x++)
            {
                unsigned int index = (1+max(x,y))*max(x,y)/2 + min(x,y);
                unsigned char c1 = clustalOrder[x];
                unsigned char c2 = clustalOrder[y];
                mat(RESIDUE_INDEX[c1], RESIDUE_INDEX[c2]) = ((double)gon350[index]) / 10;
            }
        }
        break;
    case ORIGINAL:
        if(isWeightedMatrixFile(fileName))
        {
            loadWeightedMatrix(fileName);
            d = true;
        }else if(isTriangleDistanceMatrixFile(fileName))
        {
            loadTriangleFormat(fileName, mat);
            d = true;
        }else{
            loadEMBOSSformat(fileName, mat);
        }
        break;
    default: // ファイルから読み込み
        assert(true && "[ScoreMatrix::loadMatrix]");
        ;
    }

    // 拡張アミノ酸表記B=D or N, Z=E or Q, Xへの対応．ギャップとの距離以外を平均化
    if(averageExtendedAA)
    {
        for(int i=0; i<RANK; i++)
        {
            // B,Z,Xの平均化の順番で結果が若干変わるが気にしない
            // M(B,aa) = (M(D,aa) + M(N,aa)) / 2
            cout << RESIDUE_INDEX['B'] << endl;
            if( i != RESIDUE_INDEX['B'] )
            {
                mat(RESIDUE_INDEX['B'], i) = (mat(RESIDUE_INDEX['D'], i) + mat(RESIDUE_INDEX['N'], i)) / 2;
                mat(i, RESIDUE_INDEX['B']) = mat(RESIDUE_INDEX['B'], i);
            }
            // M(Z,aa) = (M(E,aa) + M(Q,aa)) / 2
            if( i != RESIDUE_INDEX['Z'] )
            {
                mat(RESIDUE_INDEX['Z'], i) = (mat(RESIDUE_INDEX['E'], i) + mat(RESIDUE_INDEX['Q'], i)) / 2;
                mat(i, RESIDUE_INDEX['Z']) = mat(RESIDUE_INDEX['Z'], i);
            }
            // M(X,aa) = ΣM(j,aa) / (RANK-1)   jはX以外
            if( i != RESIDUE_INDEX['X'] )
            {
                double sum=0.0;
                for(int j=0; j<RANK; j++) if(j != RESIDUE_INDEX['X']) sum += mat(j, i);
                mat(RESIDUE_INDEX['X'], i) = sum / (RANK-1);
                mat(i, RESIDUE_INDEX['X']) = mat(RESIDUE_INDEX['X'], i);
            }
        }
    }

//	if(isWeightedMatrixFile(fileName))
//	{
//		// REMARK: should not normalize when this is weighted matrix, because of double normalizing.
//	}else{
        if(d == false && translateToDistance)
        {
            normalize(mat, score_min, score_max, score_ave, TRANS_A, TRANS_B);
            distance = true;
            normalized = true;
        }else if(d == true){
            distance = true;
            if(isWeightedMatrixFile(fileName))
            {
                normalized = true;
            }else{
                // REMARK: MRAMA1, MRAMA2
                normalized = false;
            }
        }else{ // d == false && translateToDistance == false
            distance = false;
            normalized = false;
        }
//	}
    return true;
}


void ScoreMatrix::normalize(AminoMatrix& mat, double& score_min, double& score_max, double& score_ave, double& TRANS_A, double& TRANS_B)
{
    getNormalizedParam(mat, score_min, score_max, score_ave, TRANS_A, TRANS_B);

    // score to disntance conversion.
    for(int y=0; y<RANK; y++)
    {
        for(int x=0; x<RANK; x++)
        {
            mat(x,y) = score_to_dist(mat(x,y), 1);

            // --- Mix Ramachandran distance --- (Only for distance mode)
            assert( mat != rmatrix );
            mat(x,y) = (1.0-gamma)*mat(x,y) + gamma*rmatrix(x,y);
        }
    }
}



bool ScoreMatrix::testScore(double target, double base)
{
    if((distance == true && target <= base) || (distance == false && target >= base))
    {
        return true;
    }else{
        return false;
    }
}


// スコアの和にも対応
// D1+D2=(S1A+B)+(S2A+B)=(S1+S2)A+2B
double ScoreMatrix::score_to_dist(double score, int pairs)
{
    double ret = score * TRANS_A + TRANS_B * pairs;
    //if(ret < 0){cout << "warning: fix " << ret << " to 0." << endl; ret = 0;}
    //if(ret > 1){cout << "warning: fix " << ret << " to 1." << endl; ret = 1;}
    //assert(distance_min <= ret && ret <= distance_max);

    return ret;
}

double ScoreMatrix::dist_to_score(double distance, int pairs)
{
    // TODO: FIX!! rama terms at dist_to_score, score_to_dist
    // distance = (1-gamma)*preDist + gamma*rama
    // preDist = (distance - gamma*rama)/(1-gamma)
    assert(false && "[ScoreMatrix::dist_to_score] Should be fixed.");

    double ret = (distance - TRANS_B * pairs) / TRANS_A;
    return ret;
}

void ScoreMatrix::getNormalizedParam(AminoMatrix& mat, double& score_min, double& score_max, double& score_ave, double& TRANS_A, double& TRANS_B)
{
    // スコアの最大値・最小値を検索
    score_min = mat.minCoeff();
    score_max = mat.maxCoeff();
    score_ave = 0;
    for(int y=0; y<RANK; y++)
    {
        for(int x=0; x<RANK; x++)
        {
            if(y<x)score_ave+=mat(x,y);
        }
    }
    score_ave /= (RANK*(RANK-1))/2;

    /* === スコア値を距離値に変換(線形変換) ===
     * アルゴリズム的バグ：TRANS_B = 0 以外では正しく最適解が求められなくなる
     *距離=スコア*a+bとしたとき
     *score min <-> dist max (1)
     *score max <-> dist min (0)
     *傾きa = (Dmin - Dmax) / (Smax - Smin)
     *切片b = -a*Smax + Dmin
     *distance = S*a+b
     *S = (distance - b)/a
     */
    if(TRANS_MODE == 0)
    {
        TRANS_A = (double)(distance_min - distance_max) / (double)(score_max - score_min);
        TRANS_B = -TRANS_A*score_max + distance_min;
    }else if(TRANS_MODE == 1){
        /* === スケーリングのみの変換 ===
         * score min <-> metric max = score min * -1/(score max - score min)
         * score max <-> metric min = score max * -1/(score max - score min)
         * Note that score max - score min = 1.
         */
        TRANS_A = -1.0 / (double)(score_max - score_min);
        TRANS_B = 0;
    }else if(TRANS_MODE == 2){
        TRANS_A = -1.0;
        TRANS_B = 0;
    }
}


bool ScoreMatrix::isTriangleDistanceMatrixFile(const string& fileName)
{
    ifstream ifs(fileName.c_str());
    if(!ifs) return false;
    string buff;
    assert(ifs);
    getline(ifs, buff);
    return (buff == "# This is a triangle format ramachandran matrix type 1 for MTRAP.");
}


// ギャップ文字の統一
char ScoreMatrix::fixAA(const char c)
{
    if(c == '-' || c == '*' || c == '~')
    {
        return GAP;
    }

    return c;
}


void ScoreMatrix::showScoreMatrix(AminoMatrix& mat)
{
    stringstream header;
    header << "  ";
    for(int c=0; c<RANK; c++)
    {
        header << left << setw(5) << ScoreMatrix::RESIDUE[c];
    }
    putLog(header.str());

    for(int r=0; r<RANK; r++)
    {
        stringstream line;
        line << ScoreMatrix::RESIDUE[r] << ' ';
        for(int c=0; c<RANK; c++)
        {
            line << left << setw(5) << setprecision(2) << mat(r,c);
        }
        putLog(line.str());
    }

}


//! string to MODE
/*!

    \param name

  \return
*/
ScoreMatrix::MODE ScoreMatrix::StringToMODE(const string& name)
{
    if(name == "CID") return ScoreMatrix::CID;
    else if(name == "CPAM") return ScoreMatrix::CPAM;
    else if(name == "CGONNET40") return ScoreMatrix::CGONNET40;
    else if(name == "CGONNET80") return ScoreMatrix::CGONNET80;
    else if(name == "CGONNET120") return ScoreMatrix::CGONNET120;
    else if(name == "CGONNET250") return ScoreMatrix::CGONNET250;
    else if(name == "CGONNET300") return ScoreMatrix::CGONNET300;
    else if(name == "CGONNET350") return ScoreMatrix::CGONNET350;
    else return ScoreMatrix::ORIGINAL;
}

std::string ScoreMatrix::translate_CDS_to_AA(const std::string& name, const std::string& encoded_cds, const char trans_table[])
{
    map<ScoreMatrix::EncodedResidue,int> NA_index;
    // NA
    NA_index[ScoreMatrix::encode('T')] = 0;
    NA_index[ScoreMatrix::encode('C')] = 1;
    NA_index[ScoreMatrix::encode('A')] = 2;
    NA_index[ScoreMatrix::encode('G')] = 3;
    NA_index[ScoreMatrix::encode('U')] = 0;
    // gap
    NA_index[ScoreMatrix::encode(ScoreMatrix::GAP)] = -1;
    // extended NA symbols
    NA_index[ScoreMatrix::encode('R')] = -1;
    NA_index[ScoreMatrix::encode('Y')] = -1;
    NA_index[ScoreMatrix::encode('M')] = -1;
    NA_index[ScoreMatrix::encode('K')] = -1;
    NA_index[ScoreMatrix::encode('S')] = -1;
    NA_index[ScoreMatrix::encode('W')] = -1;
    NA_index[ScoreMatrix::encode('B')] = -1;
    NA_index[ScoreMatrix::encode('H')] = -1;
    NA_index[ScoreMatrix::encode('V')] = -1;
    NA_index[ScoreMatrix::encode('D')] = -1;
    NA_index[ScoreMatrix::encode('N')] = -1;

    string encodedAAseq = "";

    if(encoded_cds.size() % 3 != 0)
    {
        putLog(string("Warning: the sequence ") + name + string(" may not be CDS (protein coding DNA sequence). sequence length should be a multiple of 3."));
    }

    const int multiple3 = (encoded_cds.size() / 3) * 3;
    for(int idx=0; idx<multiple3; idx+=3)
    {
        const int index1 = NA_index[encoded_cds.at(idx+0)];
        const int index2 = NA_index[encoded_cds.at(idx+1)];
        const int index3 = NA_index[encoded_cds.at(idx+2)];
        const int index = (index1 == -1 || index2 == -1 || index3 == -1) ? (-1) : (index1 * 16 + index2 * 4 + index3);
        const ScoreMatrix::EncodedResidue AA = (index != -1) ? ScoreMatrix::encode(trans_table[index]) : ScoreMatrix::encode('X');
        encodedAAseq += AA;
    }

    if(encoded_cds.size() % 3 != 0) encodedAAseq += ScoreMatrix::encode('X');

    return encodedAAseq;
}


//==============================================================================
// 変数
//==============================================================================
namespace ScoreMatrix {
    MODE mode;										// 現在のスコア行列
    MODE rmode;										// Ramachandran
    bool distance;								// 距離モードかを表すフラグ
    bool normalized;

    AminoMatrix amatrix;
    AminoMatrix rmatrix;
    AminoMatrix pmatrix;
    double score_min;
    double score_max;
    double score_ave;
    bool averageExtendedAA;
    double gamma;

    double TRANS_A;	//!< ギャップスコア-距離変換で利用するための傾き
    double TRANS_B;	//!< ギャップスコア-距離変換で利用するための切片
    int TRANS_MODE; // 0: A,B両方あり 1:幅1, B=0  2: 反転のみA=-1, B=0
}

//==============================================================================
// 定数
//==============================================================================
ScoreMatrix::EncodedResidue ScoreMatrix::RESIDUE_INDEX[256^sizeof(char)];

const char ScoreMatrix::GAP = '-';
const double ScoreMatrix::NQ = 1<<31;

const double ScoreMatrix::distance_min = 0.0;
const double ScoreMatrix::distance_max = 1.0;

// Reference of the genetic code: http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
// Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
// Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
// Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
const char ScoreMatrix::geneticCode_AAs[4][65] =
{
    // here X means STOP codon
    "FFLLSSSSYYXXCCXWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", // 0: The Standard Code (transl_table=1)
    "FFLLSSSSYYXXCCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSXXVVVVAAAADDEEGGGG", // 1: The Vertebrate Mitochondrial Code (transl_table=2)
    "FFLLSSSSYYXXCCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG", // 2: The Invertebrate Mitochondrial Code (transl_table=5)
    "FFLLSSSSYYXXCCXWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", // 3: The Bacterial, Archaeal and Plant Plastid Code (transl_table=11)
};
const char ScoreMatrix::geneticCode_Starts[4][65] =
{
    "---M---------------M---------------M----------------------------",
    "--------------------------------MMMM---------------M------------",
    "---M----------------------------MMMM---------------M------------",
    "---M---------------M------------MMMM---------------M------------",
};

const char ScoreMatrix::RESIDUE[RANK+3] =      "ARNDCQEGHILKMFPSTWYVBZX-@"; // 23 amino acids + gap + prefix + null
const char ScoreMatrix::clustalOrder[RANK+1] = "ABCDEFGHIKLMNPQRSTVWXYZ"; //23 amino acids + null
const int ScoreMatrix::PAM[RANK][RANK] =
{ // A   R   N   D   C   Q   E
    {  5, -4, -2, -1, -4, -2, -1,0,-4,-2,-4,-4,-3,-6,0,1,1,-9,-5,-1,-1,-1,-2},
    { -4,  8, -3, -6, -5,  0, -5,-6,0,-3,-6,2,-2,-7,-2,-1,-4,0,-7,-5,-4,-2,-3},
    { -2, -3,  6,  3, -7, -1,  0,-1,1,-3,-5,0,-5,-6,-3,1,0,-6,-3,-5,5,-1,-2},
    { -1, -6,  3,  6, -9,  0,  3,-1,-1,-5,-8,-2,-7,-10,-4,-1,-2,-10,-7,-5,5,2,-3},
    { -4, -5, -7, -9,  9, -9, -9,-6,-5,-4,-10,-9,-9,-8,-5,-1,-5,-11,-2,-4,-8,-9,-6},
    { -2,  0, -1,  0, -9,  7,  2,-4,2,-5,-3,-1,-2,-9,-1,-3,-3,-8,-8,-4,-1,5,-2},
    { -1, -5,  0,  3, -9,  2,  6,-2,-2,-4,-6,-2,-4,-9,-3,-2,-3,-11,-6,-4,2,5,-3},
    {  0, -6, -1, -1, -6, -4, -2,6,-6,-6,-7,-5,-6,-7,-3,0,-3,-10,-9,-3,-1,-3,-3},
    { -4,  0,  1, -1, -5,  2, -2,-6,8,-6,-4,-3,-6,-4,-2,-3,-4,-5,-1,-4,0,1,-3},
    { -2, -3, -3, -5, -4, -5, -4,-6,-6,7,1,-4,1,0,-5,-4,-1,-9,-4,3,-4,-4,-3},
    { -4, -6, -5, -8,-10, -3, -6,-7,-4,1,6,-5,2,-1,-5,-6,-4,-4,-4,0,-6,-4,-4},
    { -4,  2,  0, -2, -9, -1, -2,-5,-3,-4,-5,6,0,-9,-4,-2,-1,-7,-7,-6,-1,-2,-3},
    { -3, -2, -5, -7, -9, -2, -4,-6,-6,1,2,0,10,-2,-5,-3,-2,-8,-7,0,-6,-3,-3},
    { -6, -7, -6,-10, -8, -9, -9,-7,-4,0,-1,-9,-2,8,-7,-4,-6,-2,4,-5,-7,-9,-5},
    {  0, -2, -3, -4, -5, -1, -3,-3,-2,-5,-5,-4,-5,-7,7,0,-2,-9,-9,-3,-4,-2,-3},
    {  1, -1,  1, -1, -1, -3, -2,0,-3,-4,-6,-2,-3,-4,0,5,2,-3,-5,-3,0,-2,-1},
    {  1, -4,  0, -2, -5, -3, -3,-3,-4,-1,-4,-1,-2,-6,-2,2,6,-8,-4,-1,-1,-3,-2},
    { -9,  0, -6,-10,-11, -8,-11,-10,-5,-9,-4,-7,-8,-2,-9,-3,-8,13,-3,-10,-7,-10,-7},
    { -5, -7, -3, -7, -2, -8, -6,-9,-1,-4,-4,-7,-7,4,-9,-5,-4,-3,9,-5,-4,-7,-5},
    { -1, -5, -5, -5, -4, -4, -4,-3,-4,3,0,-6,0,-5,-3,-3,-1,-10,-5,6,-5,-4,-2},
    { -1, -4,  5,  5, -8, -1,  2,-1,0,-4,-6,-1,-6,-7,-4,0,-1,-7,-4,-5,5,1,-2},
    { -1, -2, -1,  2, -9,  5,  5,-3,1,-4,-4,-2,-3,-9,-2,-2,-3,-10,-7,-4,1,5,-3},
    { -2, -3, -2, -3, -6, -2, -3,-3,-3,-3,-4,-3,-3,-5,-3,-1,-2,-7,-5,-2,-2,-3,-3},
};

const short ScoreMatrix::gon40[]={
  92,
   0,   0,
 -31,   0, 163,
 -56,   0,-135, 111,
 -37,   0,-140,  16, 105,
 -92,   0, -64,-152,-143, 126,
 -32,   0, -91, -51, -76,-152, 105,
 -65,   0, -67, -41, -40, -50, -81, 145,
 -76,   0, -87,-150,-106, -39,-158, -94, 104,
 -54,   0,-132, -47, -13,-127, -79, -34, -86, 103,
 -68,   0, -85,-155,-108, -13,-141, -85,   5, -85,  89,
 -45,   0, -63,-130, -80, -16,-114, -60,  10, -57,  16, 140,
 -62,   0, -83,   6, -38,-104, -40,  -7, -99, -20,-112, -91, 115,
 -37,   0,-137, -69, -60,-128, -87, -71,-108, -62, -83,-119, -78, 124,
 -43,   0,-113, -32,  10,-100, -71,   0, -91,   2, -60, -35, -25, -46, 118,
 -61,   0, -86, -77, -50,-130, -69, -31,-103,  19, -84, -81, -47, -73,  -6, 112,
   0,   0, -35, -36, -41,-111, -37, -48, -95, -43, -95, -64, -11, -35, -35, -51,  99,
 -25,   0, -59, -47, -52, -90, -85, -46, -51, -34, -78, -44, -27, -42, -39, -52,  13, 100,
 -22,   0, -43,-133, -74, -58,-122, -98,  28, -82, -18, -22,-103, -86, -79, -88, -74, -25,  97,
-120,   0, -68,-171,-131,  -6,-108, -70, -93,-127, -71, -72,-119,-149, -87, -63, -98,-120,-115, 181,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
 -95,   0, -56, -98,-107,  31,-129,   5, -76, -88, -64, -66, -62,-106, -81, -75, -69, -87, -73,   1,   0, 135,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0};

const short ScoreMatrix::gon80[]={
  75,
   0,   0,
 -10,   0, 154,
 -31,   0, -93,  96,
 -17,   0, -94,  31,  88,
 -64,   0, -39,-111,-102, 114,
 -11,   0, -61, -26, -47,-115,  97,
 -39,   0, -43, -17, -17, -26, -53, 127,
 -43,   0, -54,-106, -73, -15,-114, -64,  86,
 -30,   0, -88, -21,   4, -89, -50, -12, -59,  85,
 -43,   0, -55,-109, -75,   7,-104, -57,  22, -58,  77,
 -26,   0, -39, -88, -53,   3, -83, -38,  25, -37,  31, 117,
 -34,   0, -55,  21, -13, -75, -18,   9, -71,  -2, -79, -62,  97,
 -16,   0, -93, -42, -35, -93, -58, -45, -75, -37, -58, -78, -48, 114,
 -22,   0, -76,  -9,  23, -70, -44,  14, -60,  17, -39, -19,  -6, -24,  95,
 -36,   0, -60, -44, -23, -90, -43, -10, -71,  33, -58, -53, -22, -45,  11,  97,
  14,   0, -15, -14, -19, -77, -16, -25, -62, -20, -64, -41,   5, -14, -15, -27,  78,
  -5,   0, -34, -24, -27, -62, -52, -24, -28, -15, -49, -25,  -7, -20, -18, -27,  25,  81,
  -6,   0, -21, -89, -51, -31, -86, -65,  41, -54,   3,   1, -69, -57, -51, -60, -43,  -9,  80,
 -87,   0, -43,-124, -98,  16, -81, -43, -63, -89, -44, -45, -86,-112, -62, -41, -72, -87, -80, 173,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
 -65,   0, -32, -69, -74,  49, -94,  21, -47, -60, -35, -37, -39, -76, -53, -50, -46, -58, -47,  23,   0, 123,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0};

const short ScoreMatrix::gon120[]={
  59,
   0,   0,
  -1,   0, 144,
 -18,   0, -69,  82,
  -9,   0, -68,  35,  72,
 -48,   0, -26, -87, -78, 102,
  -3,   0, -45, -14, -31, -92,  90,
 -26,   0, -31,  -7,  -6, -14, -37, 110,
 -27,   0, -36, -80, -55,  -3, -87, -48,  72,
 -19,   0, -64,  -8,  11, -67, -34,  -2, -44,  69,
 -30,   0, -39, -82, -57,  15, -82, -42,  28, -44,  66,
 -17,   0, -26, -64, -40,  11, -65, -28,  29, -27,  34,  95,
 -20,   0, -41,  26,  -1, -58,  -7,  14, -55,   5, -61, -46,  80,
  -6,   0, -68, -28, -22, -72, -41, -31, -56, -24, -44, -56, -32, 105,
 -12,   0, -56,   1,  25, -53, -30,  17, -43,  20, -30, -14,   1, -14,  74,
 -23,   0, -45, -27, -10, -68, -30,  -1, -53,  36, -44, -38, -10, -30,  16,  83,
  16,   0,  -7,  -5,  -9, -58,  -6, -14, -44, -10, -47, -29,  10,  -5,  -7, -15,  60,
   2,   0, -21, -13, -15, -47, -35, -14, -17,  -6, -34, -16,   0, -10,  -9, -16,  26,  64,
   0,   0, -11, -65, -38, -17, -65, -47,  42, -39,  13,  10, -50, -42, -36, -44, -28,  -3,  65,
 -68,   0, -29, -96, -78,  27, -66, -28, -46, -68, -29, -31, -68, -89, -49, -30, -57, -67, -59, 166,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
 -48,   0, -20, -53, -56,  55, -74,  26, -31, -44, -20, -22, -28, -59, -38, -37, -35, -42, -33,  33,   0, 111,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0};

const short ScoreMatrix::gon160[]={
  46,
   0,   0,
   3,   0, 135,
 -11,   0, -53,  70,
  -4,   0, -52,  34,  59,
 -38,   0, -18, -70, -62,  91,
   2,   0, -34,  -7, -21, -76,  82,
 -18,   0, -23,  -1,  -1,  -7, -27,  93,
 -18,   0, -25, -62, -43,   3, -70, -37,  59,
 -12,   0, -48,  -1,  13, -53, -24,   2, -35,  55,
 -22,   0, -29, -65, -45,  19, -67, -32,  30, -34,  57,
 -12,   0, -19, -50, -31,  14, -52, -21,  29, -21,  34,  76,
 -12,   0, -31,  26,   5, -47,  -2,  15, -44,   8, -48, -36,  65,
  -1,   0, -52, -19, -14, -58, -30, -22, -43, -16, -35, -42, -22,  96,
  -7,   0, -42,   6,  23, -41, -21,  17, -32,  20, -24, -12,   5,  -8,  56,
 -16,   0, -35, -16,  -3, -53, -21,   3, -41,  35, -35, -29,  -4, -21,  17,  71,
  16,   0,  -2,   0,  -3, -45,  -1,  -8, -33,  -4, -36, -23,  11,   0,  -2,  -9,  44,
   5,   0, -14,  -6,  -8, -36, -24,  -8, -12,  -2, -24, -11,   3,  -4,  -4,  -9,  23,  50,
   1,   0,  -6, -49, -30,  -8, -52, -35,  40, -30,  17,  14, -38, -32, -27, -34, -20,   0,  53,
 -55,   0, -21, -78, -64,  32, -55, -19, -34, -54, -20, -22, -55, -74, -40, -24, -47, -54, -45, 158,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
 -37,   0, -13, -42, -44,  56, -60,  27, -20, -35, -11, -13, -22, -48, -29, -29, -28, -32, -24,  38,   0, 100,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0};

const short ScoreMatrix::gon250[]={ // オリジナルを40倍にしたもの
  24,
   0,   0,
   5,   0, 115,
  -3,   0, -32,  47,
   0,   0, -30,  27,  36,
 -23,   0,  -8, -45, -39,  70,
   5,   0, -20,   1,  -8, -52,  66,
  -8,   0, -13,   4,   4,  -1, -14,  60,
  -8,   0, -11, -38, -27,  10, -45, -22,  40,
  -4,   0, -28,   5,  12, -33, -11,   6, -21,  32,
 -12,   0, -15, -40, -28,  20, -44, -19,  28, -21,  40,
  -7,   0,  -9, -30, -20,  16, -35, -13,  25, -14,  28,  43,
  -3,   0, -18,  22,   9, -31,   4,  12, -28,   8, -30, -22,  38,
   3,   0, -31,  -7,  -5, -38, -16, -11, -26,  -6, -23, -24,  -9,  76,
  -2,   0, -24,   9,  17, -26, -10,  12, -19,  15, -16, -10,   7,  -2,  27,
  -6,   0, -22,  -3,   4, -32, -10,   6, -24,  27, -22, -17,   3,  -9,  15,  47,
  11,   0,   1,   5,   2, -28,   4,  -2, -18,   1, -21, -14,   9,   4,   2,  -2,  22,
   6,   0,  -5,   0,  -1, -22, -11,  -3,  -6,   1, -13,  -6,   5,   1,   0,  -2,  15,  25,
   1,   0,   0, -29, -19,   1, -33, -20,  31, -17,  18,  16, -22, -18, -15, -20, -10,   0,  34,
 -36,   0, -10, -52, -43,  36, -40,  -8, -18, -35,  -7, -10, -36, -50, -27, -16, -33, -35, -26, 142,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
 -22,   0,  -5, -28, -27,  51, -40,  22,  -7, -21,   0,  -2, -14, -31, -17, -18, -19, -19, -11,  41,   0,  78,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0};

const short ScoreMatrix::gon300[]={
  16,
   0,   0,
   5,   0, 104,
  -1,   0, -24,  37,
   1,   0, -23,  23,  27,
 -18,   0,  -5, -37, -31,  60,
   5,   0, -15,   3,  -4, -42,  58,
  -6,   0, -10,   5,   4,   0, -10,  45,
  -6,   0,  -7, -30, -21,  11, -36, -16,  33,
  -2,   0, -21,   6,  11, -26,  -7,   5, -17,  24,
  -9,   0, -10, -32, -22,  19, -36, -14,  25, -17,  33,
  -5,   0,  -6, -24, -16,  15, -28, -10,  22, -11,  24,  31,
  -1,   0, -14,  18,   9, -25,   5,  10, -22,   8, -24, -17,  27,
   3,   0, -23,  -4,  -2, -30, -11,  -8, -20,  -3, -18, -19,  -6,  66,
  -1,   0, -18,   9,  14, -20,  -6,   9, -15,  13, -13,  -8,   7,  -1,  18,
  -4,   0, -17,   0,   5, -25,  -6,   6, -19,  22, -18, -13,   4,  -6,  13,  37,
   8,   0,   1,   5,   3, -22,   4,  -1, -14,   2, -17, -11,   7,   4,   2,   0,  15,
   5,   0,  -3,   1,   1, -17,  -7,  -1,  -4,   2,  -9,  -5,   4,   2,   1,  -1,  11,  17,
   0,   0,   1, -23, -15,   4, -26, -15,  26, -13,  17,  15, -17, -14, -12, -15,  -8,   0,  26,
 -29,   0,  -7, -42, -36,  36, -34,  -5, -13, -28,  -4,  -6, -30, -41, -23, -14, -27, -28, -19, 132,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
 -17,   0,  -3, -22, -22,  46, -33,  18,  -3, -17,   3,   1, -12, -25, -14, -14, -15, -15,  -7,  40,   0,  67,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0};

const short ScoreMatrix::gon350[]={
  10,
   0,   0,
   4,   0,  93,
   0,   0, -19,  29,
   1,   0, -17,  19,  20,
 -14,   0,  -3, -30, -25,  51,
   5,   0, -12,   4,  -2, -35,  51,
  -4,   0,  -8,   5,   4,   1,  -7,  33,
  -4,   0,  -5, -24, -17,  11, -29, -13,  27,
  -1,   0, -16,   6,   9, -21,  -4,   5, -13,  18,
  -7,   0,  -7, -25, -18,  18, -30, -11,  22, -14,  28,
  -4,   0,  -4, -19, -13,  14, -23,  -8,  19,  -9,  21,  23,
   0,   0, -11,  15,   9, -20,   5,   8, -18,   7, -19, -14,  20,
   3,   0, -18,  -2,   0, -25,  -7,  -5, -16,  -2, -15, -14,  -3,  56,
   0,   0, -14,   8,  11, -16,  -4,   7, -11,  10, -11,  -7,   6,   0,  12,
  -2,   0, -13,   2,   6, -20,  -4,   6, -15,  18, -14, -11,   4,  -4,  10,  28,
   6,   0,   1,   5,   3, -18,   5,   0, -11,   2, -13,  -9,   6,   4,   2,   1,  10,
   4,   0,  -2,   2,   1, -13,  -5,  -1,  -3,   2,  -7,  -4,   4,   2,   1,   0,   8,  11,
   0,   0,   2, -18, -12,   5, -21, -11,  22, -10,  16,  14, -13, -11,  -9, -12,  -6,   0,  21,
 -24,   0,  -4, -35, -29,  35, -30,  -3,  -9, -23,  -1,  -3, -24, -34, -19, -12, -22, -23, -14, 124,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
 -14,   0,  -1, -18, -17,  42, -27,  15,  -1, -14,   5,   2, -10, -20, -11, -12, -12, -12,  -4,  39,   0,  57,
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0};

const int ScoreMatrix::DEF_PAM_GAP_BEGIN = -11;
const int ScoreMatrix::DEF_PAM_GAP_EXT = -1;

