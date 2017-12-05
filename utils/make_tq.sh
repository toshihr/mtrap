#!/bin/sh
# === DIRECTORY ===
# 出力先
outdir="./TQ"
# TEMP_DIRにはデータベースが展開される．読み込みの高速な場所あるいは別デバイスにすると計算が速くなる．
TEMP_DIR="./temp"
# データベースの格納されている場所
cd
homeDir="`pwd`"
cd - >/dev/null
DATABASE_DIR="${homeDir}/research/database"
# === MODE ===
# 0: default, 1: skip gap-gap site 2: cancel gap-* site 3: TQ=prob. 4: without lowercase 5: Global TQ, f: 対象化を頻度の段階で行う
# 2,3,4は正しい動作をしているか確認が必要
# モードは複数同時指定可: 123 等
opt="0n"
# 拡張文字を平均化するかどうか。しないほうが若干性能がよさそう。
#ave="-averageExtendedAA"
ave=""
# === FILTER ===
filter="*.ref_fasta"
#filter="BBS*.ref_fasta"
#filter="BB[^S]*.ref_fasta"
# === DATABASE ===
# RV40はBBSデータ(より相同な部分配列)が存在しない！BBデータ(フル)のみ
database=""
#database="$database bali3pdbm"
#database="$database dali_z50_id40"
#database="$database mini"
#database="$database mini2"
#database="$database minipairwise"
#database="$database minimulti_with_core"
#database="$database minimulti_only_full"
#database="$database ox"
#database="$database oxm"
#database="$database oxx"
#database="$database prefab4refm"
#database="$database sabre"
#database="$database sabrem"
#database="$database RV11"
#database="$database RV12"
#database="$database RV20"
#database="$database RV30"
#database="$database RV40"
#database="$database RV50"
#database="$database RV911"
#database="$database RV912"
#database="$database RV913"
#database="$database prefab4_training"
#database="$database prefab4_test"
#database="$database prefab4"
#database="$database homstrad_training"
#database="$database homstrad_test"
database="$database homstrad"
#database="$database homstrad2010Mar1"
#database="$database IRMBaseAll"
#database="$database IRMBaseRef1"
#database="$database IRMBaseRef2"
#database="$database IRMBaseRef3"
#database="$database SABmark_twi"
#database="$database SABmark_sup"
#database="$database SABmark1.63_sup"
#===RNA===
# database="$database bralibase2_SRP"
# database="$database bralibase2_U5"
# database="$database bralibase2_g2intron"
# database="$database bralibase2_rRNA"
# database="$database bralibase2_tRNA"
# database="$database bralibase2_SRP_p"
# database="$database bralibase2_U5_p"
# database="$database bralibase2_g2intron_p"
# database="$database bralibase2_rRNA_p"
# database="$database bralibase2_tRNA_p"
# database="$database bralibase2_SRP_p_test"
# database="$database bralibase2_U5_p_test"
# database="$database bralibase2_g2intron_p_test"
# database="$database bralibase2_rRNA_p_test"
# database="$database bralibase2_tRNA_p_test"
# database="$database bralibase2_SRP_p_training"
# database="$database bralibase2_U5_p_training"
# database="$database bralibase2_g2intron_p_training"
# database="$database bralibase2_rRNA_p_training"
# database="$database bralibase2_tRNA_p_training"

# === TODO ===
# 将来的に重み付きカウントに対応させる
# grep ">" ./*.ref_fasta | cut -d'>' -f2 | sort | uniq -c などを利用して
# 作成した配列頻度リスト(名前->頻度の一覧)を渡す．今は""を渡している
# === BEGIN ===

make_list (){
	for i in $1 ; do
		find $i -type f -name "$filter" >> $2
	done
}

# === UNPACK TARGET DATABASE ===
if [ ! -d "$TEMP_DIR" ] ; then
	mkdir "$TEMP_DIR"
fi
for targetdatabase in $database ; do
	cd "$TEMP_DIR"
  tar xf "${DATABASE_DIR}/${targetdatabase}.tar.gz"
	cd - >/dev/null
done

# === GENERATE TQ ===
if [ ! -d "$outdir" ] ; then
	mkdir "$outdir"
fi

allName="all"
listFile="_list.tmp"
freqFile="_freq.tmp"
printf "" > $listFile
for targetdatabase in $database ; do
	localListFile="_list.tmp2"
	localFreqFile="_freq.tmp2"
	printf "" > "$localListFile"
	make_list "$TEMP_DIR/$targetdatabase/ref" "$localListFile"
	# 各配列の頻度を計算
	cat "$localListFile" | xargs grep ">" | cut -d'>' -f2 | sort | uniq -c | sort -nr | awk '{print $2 "," $1}' > "$localFreqFile"

	mtrap $ave -count "$opt" "$localListFile" "$localFreqFile" "${outdir}/${targetdatabase}_weighted_`date +"%Y%m%d"`"
	mtrap $ave -count "$opt" "$localListFile" "" "${outdir}/${targetdatabase}_`date +"%Y%m%d"`"
	cat "$localListFile" >> $listFile
	# rm "$localListFile"
  #   rm "$localFreqFile"
	allName="${allName}_${targetdatabase}"
done
# ALL DATABASES
#cat "$listFile" | xargs grep ">" | cut -d'>' -f2 | sort | uniq -c | sort -nr | awk '{print $2 "," $1}' > "$freqFile"
#mtrap $ave -count "$opt" "$listFile" "$freqFile" "${outdir}/${allName}_weighted_`date +"%Y%m%d"`"
#mtrap $ave -count "$opt" "$listFile" "" "${outdir}/${allName}_`date +"%Y%m%d"`"
#rm "$listFile"
#rm "$freqFile"

# === REMOVE TEMP. DATA ===
cd "$TEMP_DIR"
for targetdatabase in $database ; do
  rm -r "./${targetdatabase}"
done
cd - >/dev/null
