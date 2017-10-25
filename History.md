## 20171025 v2.0.0
* Switch repogitory to GitHub.
* Remove windows support.

## 20130129
* Improved: CDS modeにおいてgenetic codeを配列ごとに指定できるように機能を拡張

## 20130124
* Fixed: CDS modeを正しく動くよう修正．入力配列中の-,STOP CODONはすべてXにすることで対応

## 20130116
* Added: アミノ酸, DNA, RNA に存在しないシンボルが入力配列にあるときに，エラーで止まるように修正

## 20121130
* Changed: Family estimationをデフォルトでOffに. 結果がよくないため. そのうちFold,ClassレベルMatrixを使って再挑戦する
* Added: CDS("protein" coding DNA sequence)のサポート. これは売りになるかも.

## 20121126
* Added: Family estimation
* TODO: blastの存在チェック. 存在しないときはfamily estimationをoffにする

## 20121124
* TODO: Blocks14.3でmatrix,TQを学習してみる。もしかしたらVTML200をoutperformするかも

## 20121120 Major fixed point
* Changed: デフォルト値をsmatrix,pmatrix=VTML200I,beta2=1.5と変更(homstradでベストスコア)
* ただしスピードとのバランスを考えると-ph=1.0でもよいかも

## 20121119
* Improved: -pf 0.0の時はpartition functionを計算しないようにし高速化
* Bug fixed: -pmオプションが反応しないのを修正

## 20121113
* pathMatrixをIntegerで計算するようにした. homstradで5sec程度の高速化
* 重大バグの修正！！pathMatrixの計算にてINFの扱いがおかしかったのを修正(分配関数のほうは未修正)
* バグ修正の結果、若干homstradでのスコアが上昇した.
* partition functionをfloatにしたところ20sec近い大幅な高速化. でも1%弱精度は下がった

## 20121112
* Clang, VC, GCCに対応するコードに微調整
* CPackによるインストーラの調整
* ソースコードをパックするスクリプトを作成
* 32bit compile with 64bit environmentをサポート
* TODO: mac環境でもopen mpを

## 20121108  ver.2.0-rc0
* ソースコードをバイナリ公開へ向けて整理
* オプション表記を整理
* -msaprobsオプションを-postprobへ名前を変更
* cmakeによるconfig.hの生成フォルダを調整
* openmpライブラリの存在を自動認識するように変更
* openmpを有効にした際に正しく動くように調整
* VC++10はUTF8BOMでないとC4819警告を発する
* Windows環境でのコンパイルはQTCreator2.5以上で行うことにした
* utility.h内inline SDのバグを修正
* log2,log10,roundはC99で規程されているがVC++はなぜかC99完全準拠ではないためサポートされない
* cmakeはUTF8BOMをサポートしない
* VC++はUTF8BOMしかサポートしない
* gccの最近のバージョンではUTF8のBOMあり、なし両方をサポート
* Linux,Windowsと環境を変えてコンパイルする際はcmakeclean.sh,remove_BOM_from_CMakeLists.shを呼ぶことにする
* Eigenのライセンスには要注意。Eigen自体はヘッダーライブラリだけど、これをincludeするとコピーレフト発生？？

## 20120727
* binary transition quantityをsnappy libraryを用いてサポート。圧縮サイズが1/10くらいになっておどろいた。
* ファイルを探しに行く際に、まずはそのまま開くことをトライするように変更。
* TQ作成ルーチンをブラッシュアップ
* Win環境(VCコンパイラ)にて、レジストリに記録されているインストール先をサーチパスに加える機能を実装
* TODO:RNA検索用パラメータのデフォルト値を設定するフラグを追加

## 20120413
* pmtrap::calc_passMatrixでは壁際のギャップ計算時にTQを足していないことを確認．内部では足している．今のmtrapでは壁際ギャップ計算時にもTQを足している．
* partition functionのFULL MEMORY時のバグを修正
* TODO: pmtrapと同じギャップ計算ルーチンを実装後，consistency計算の検証, iterationは乱数の問題もあり検証は難しい？

## 20120413
* Gap openバグ修正

## 20120411
* Treeバグ修正
* Profileベースのマルチプルアライメントを実装
* Consistency実装

## 20120303
* Profile(TCoffee風，頻度行列ではない)アライメントを実装
* ProgressiveProfileAlignmentを実装
* MTRAPによるProfileは事後確率ほど柔軟ではなく，経路におけるアミノ酸ペアを確率１としている．

## 20111226
* EERによるデータベース検索をサポート for 寺本君修論用

## 20111012
* 複数のアミノ酸置換行列を重ねあわせるモードを作成開始

## 20110713
* ビルドツールをCMakeへ移行．Ubuntu,Cygwin,MinGW上でのコンパイルをサポート．
* また，QT Creatorでのプロジェクトインポートが楽になった．QT Creatorプロジェクトを環境ごとに用意したい！

1. 現状，MinGWでコンパイルするとエラーがおきる．おそらくファイルパスの扱いがWindows風になるためだと思う．

  DATA_DIRをMinGW時だけ変えることで対応できる可能性あり．
2. Cygwin上でのmake installでバグる．mtrap.exeではなくmtrapをコピーしようとするため

## 20110703
* Eigenライブラリの導入．RamachandranMatrixを外部リソース化

## 20110627
* ソースコード整理 ver.1.3へ

## 20110322
* バイナリモードTQを追加．読み込みが超高速化．

## 20091125
* データベース検索機能を追加. now, read all sequences at initialize
* 伊藤君のソースへあとでマージンする

## 20091129
* 出力するアライメント値をSPmetricからmetric(Sum of Transition-quantity)に変更
* データベース検索機能を少し改良した与えられた配列の総当たりの距離を返す機能を追加(Make TCoffee Libraryでの利用を想定)
