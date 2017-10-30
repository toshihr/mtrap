#!/bin/sh
# データベース上のすべてのreference alignmentの%IDを計算しソートし
# %IDが偏らない様に半分づつのリストをつくる
# !! PairWise=1 として実行しても出力フォルダは変わらないことに注意
# USE: ./separateDatabaseToTrainingTest.sh [database name] [flag: ONLY_PAIRWISE]
# ex) ./separateDatabaseToTrainingTest.sh prefab4 0
# OUTPUT: database[_pairwise]_list, database[_pairwise]_list_training, database[_pairwise]_list_test と 対応するデータベース
database=$1
ONLY_PAIRWISE=$2

name=`basename $database`
dir=`dirname $database`
if [ $ONLY_PAIRWISE -ne 0 ] ; then
	list="${dir}/${name}_pairwise.list"
else
	list="${dir}/${name}.list"
fi

# 対象ファイル一覧を作成
printf "" > _temp.tmp
find $database/ref -type f -name "*.ref_fasta" | while read target; do
	id=`eerdist -average -id $target`
	num=`cat $target | grep ">" | wc -l`
	if [ $ONLY_PAIRWISE -ne 0 -a $num -eq 2 ] ; then
		echo "$target,$id,$num" >> _temp.tmp
	elif [ $ONLY_PAIRWISE -eq 0 ] ; then
		echo "$target,$id,$num" >> _temp.tmp
	fi
done
cat _temp.tmp | sort -t',' -k2 > $list
rm _temp.tmp

# training testに振り分け
printf "" > ${list}_training
printf "" > ${list}_test

training=0
cat $list | cut -f1 -d',' | while read target; do
	training=`expr \( $training + 1 \) % 2`
#	echo $training
	if [ $training -eq 0 ] ; then
		echo $target >> ${list}_test
	else
		echo $target >> ${list}_training
	fi
done

mkdir ${database}_training
mkdir ${database}_training/ref
mkdir ${database}_training/inputdata
mkdir ${database}_test
mkdir ${database}_test/ref
mkdir ${database}_test/inputdata

cat ${list}_training | while read target; do
#	cp ${target%.*}.* ${database}_training/
	targetname=`basename $target`
	cp ${database}/ref/${targetname%.*}.* ${database}_training/ref/
	cp ${database}/inputdata/${targetname%.*}.* ${database}_training/inputdata/
done

cat ${list}_test | while read target; do
#	cp ${target%.*}.* ${database}_test/
	targetname=`basename $target`
	cp ${database}/ref/${targetname%.*}.* ${database}_test/ref/
	cp ${database}/inputdata/${targetname%.*}.* ${database}_test/inputdata/
done
