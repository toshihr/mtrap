#pragma once
#include <map>
#include <string>
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include "utility.h"

//TODO: DNAと両対応にする
//TODO: 一次元配列で管理し速度重視にする．Eigen3は一次元で管理するのでEigenへ切り替えていく
//TODO: 重み計算済み配列を出力できるようにする
// Eigenのデフォルトはcolumn-major: column1 column2 ... と格納される
// C言語の多次元配列はrow-major: row1 row2 ...
namespace ScoreMatrix {

    enum { RANK = 23 }; // 20 + B,Z,X
// --- MODE ---
    enum MODE { CID, CPAM, CGONNET40, CGONNET80, CGONNET120, CGONNET160, CGONNET250, CGONNET300, CGONNET350, ORIGINAL };

    typedef Eigen::Matrix<double, RANK, RANK> AminoMatrix; /*!< Amino Acid Substitution Matrix */

    typedef char EncodedResidue; // same as the residue index of the matrix

    // Symbol code
    extern const char GAP;				 /*!< Gap */
    extern const double NQ;				/*!< Not a quantity */

    extern EncodedResidue RESIDUE_INDEX[256^sizeof(char)]; // Map: Residues -> EncodedReisdue (index of the matrix), sizeof(char) = 1

    extern const char RESIDUE[RANK+3]; // 23 amino acids + gap + prefix + null
    extern const char clustalOrder[RANK+1]; // 23 amino acids + null

    extern const char geneticCode_AAs[4][65];
    extern const char geneticCode_Starts[4][65];

    extern const int PAM[RANK][RANK]; // ARN... order
    extern const short gon40[];			// go40-350 use clustalOrder
    extern const short gon80[];
    extern const short gon120[];
    extern const short gon160[];
    extern const short gon250[];
    extern const short gon300[];
    extern const short gon350[];

    extern const int DEF_PAM_GAP_BEGIN;
    extern const int DEF_PAM_GAP_EXT;

    extern const double distance_min;
    extern const double distance_max;

    extern MODE mode;										// 現在のスコア行列
    extern MODE rmode;								// Ramachandran Matrix This Value should be changed by setRamachandran()
    extern int TRANS_MODE; // 0: A,B両方あり 1:幅1, B=0  2: 反転のみA=-1, B=0

    extern bool distance;								// 距離モードかを表すフラグ
    extern bool normalized;                             // 距離モードへ変換したのかを表すフラグ．distance==trueだからといってnormalized==trueとは限らない

    extern AminoMatrix amatrix;		/*!< Amino Acid genetic matrix for calculation. */
    extern AminoMatrix rmatrix; // Averaged Ramachandran quantity 将来的にはweighted matrixに移行し，これは廃止
    extern AminoMatrix pmatrix; // used for partition function

    extern bool averageExtendedAA;								 // B,Z,Xを微調整するかを表すフラグ

    // normalized == true の時のみ意味のあるパラメータ
    // TODO: Weighted Matrix mode では各マトリックスでの変換パラメータの平均値とする
    extern double score_min;
    extern double score_max;
    extern double score_ave;
    extern double TRANS_A;	//!< ギャップスコア-距離変換で利用するための傾き
    extern double TRANS_B;	//!< ギャップスコア-距離変換で利用するための切片

    extern double gamma;							// Ramachandran no omomi (1-gamma)*FinalScore + gamma*rama

    void init();
    void initTable();

    void setRamachandran(const std::string& fileName);

    void loadEMBOSSformat(const std::string& fileName, AminoMatrix& mat = ScoreMatrix::amatrix);
    void saveEMBOSSformat(const std::string& fileName, const AminoMatrix& mat, const double& iscale, const std::string& database, const double& mutualEntropy);
    void loadTriangleFormat(const std::string& fileName, AminoMatrix& mat);
    void saveTriangleFormat(const std::string& fileName, const AminoMatrix& mat);
    bool loadMatrix(const std::string& fileName, bool translateToDistance, AminoMatrix& mat);

    void normalize(AminoMatrix& mat, double& score_min, double& score_max, double& score_ave, double& TRANS_A, double& TRANS_B);
    bool testScore(double target, double base);					// スコアモードの時，target > baseでtrueを返す
    double score_to_dist(double score, int pairs);
    double dist_to_score(double distance, int pairs);
    void getNormalizedParam(AminoMatrix& mat, double& score_min, double& score_max, double& score_ave, double& TRANS_A, double& TRANS_B);
    bool isTriangleDistanceMatrixFile(const std::string& fileName);
    char fixAA(const char c);
    ScoreMatrix::MODE StringToMODE(const std::string& name);

    inline const EncodedResidue& encode(unsigned char a) { return RESIDUE_INDEX[a]; }
    inline const char& decode(EncodedResidue a) { return RESIDUE[a]; }
    inline void encode(std::string& s) { for(std::string::iterator i=s.begin(); i!=s.end(); ++i) *i = encode(*i); }
    inline void decode(std::string& s) { for(std::string::iterator i=s.begin(); i!=s.end(); ++i) *i = decode(*i); }
    std::string translate_CDS_to_AA(const std::string& name, const std::string& encoded_cds, const char trans_table[]);

    // --- DEBUG ---
    void showScoreMatrix(AminoMatrix& mat);
}

