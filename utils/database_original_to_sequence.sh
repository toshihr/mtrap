#!/bin/sh
# === DATABASE CONVERTER ===
# requires: emboss package
# 変換元データベース：下記変数にて指定(2012/2/15現在，prefab4のみ正常動作)
# 出力先：変換元データベース_converted
# scriptsフォルダで実行することが前提
# 2012/2/15 prefab4のスクリプトを作り直し

balibase="RV11 RV12 RV20 RV30 RV40 RV50"
prefab1_dir="../original_data/prefab1"
prefab2_dir="../original_data/prefab2"
prefab3_dir="../original_data/prefab3"
prefab4_dir="../original_data/prefab4"
sabmark="SABmark"
dali="./original_data/dali_z50_id40.dat"

do_dali (){
		outdir="dali_z50_id40"
		if [ ! -d $outdir ] ; then
				mkdir $outdir
		fi
		if [ ! -d $outdir/ref ] ; then
				mkdir $outdir/ref
		fi
		if [ ! -d $outdir/inputdata ] ; then
				mkdir $outdir/inputdata
		fi

		cat $dali | while read i ; do
				zscore=`echo ${i} | cut -f1 -d':'`
				id=`echo ${i} | cut -f2 -d':'`
				name1=`echo ${i} | cut -f3 -d':'`
				seq1=`echo ${i} | cut -f4 -d':'`
				name2=`echo ${i} | cut -f5 -d':'`
				seq2=`echo ${i} | cut -f6 -d':'`
				
				outname=./$outdir/ref/${name1}_${name2}.ref_fasta
				echo ">${name1}" >$outname 
				echo $seq1 >>$outname
				echo ">${name2}" >>$outname
				echo $seq2 >>$outname

				degapseq -auto $outname _temp.tmp
				cat _temp.tmp | tr "[a-z]" "[A-Z]" > ./$outdir/inputdata/${name1}_${name2}.fasta
				rm _temp.tmp
		done
}

do_balibase (){
	for database in $balibase; do
		cd ./$database

		if [ ! -d inputdata ] ; then
			mkdir inputdata
		fi

		for i in $( ls *.tfa ); do
			echo TARGET: $i

			seqret -sequence ${i%.*}.msf -outseq ${i%.*}.ref_fasta
			degapseq ${i%.*}.ref_fasta ./inputdata/${i%.*}.fasta
		done

		cd ..
	done
}


# === PREFAB4データベースについて ===
# in/* : 入力配列２本＋２本それぞれと似た配列（PSI BLASTにより取得）２４本が格納されている．
#        QScoreプログラムはtest alignmentにref alignmentにない配列があっても正しく計算する．
#        なので，入力として情報をたくさん与えた場合に，その情報をどれだけ加味できるかを見るのに利用する配列だと思われる．
# inw/* : 上の配列利用率を偏らせたバージョン
#        BAliBASEと似せるため？ 基本的には使わなくてよさそう
# ref/* : pairwise reference alignmentが格納されている．相同部位はUppercase，非相同部位はLowerCase
# 
# PREFAB2 はinとrefalnのファイル数がそもそも違うので何かが変．よくわからない．prefab2multiinputはとりあえず作らない
do_prefab_common (){
	prefab_dir="$1"
	ref_dir="$2"
	useInDir=$3

	if [ ! -d "${prefab_dir}_converted" ] ; then
		mkdir "${prefab_dir}_converted"
	fi
	if [ ! -d "${prefab_dir}_converted/inputdata" ] ; then
		mkdir "${prefab_dir}_converted/inputdata"
	fi
	if [ ! -d "${prefab_dir}_converted/ref" ] ; then
		mkdir "${prefab_dir}_converted/ref"
	fi

	find "${prefab_dir}/${ref_dir}" -type f -name "*" | while read i; do
		inFile="$i"
		outRefFile="${prefab_dir}_converted/ref/${i##*/}.ref_fasta"
		outMsfFile="${prefab_dir}_converted/ref/${i##*/}.msf"
		outInputFile="${prefab_dir}_converted/inputdata/${i##*/}.fasta"
		tmpFile="${prefab_dir}_converted/${i##*/}.temp"
		echo TARGET: $inFile

		# ref/*.ref_fasta
		# ref fastaへ変換(小文字を大文字へ変換はしない)
		seqret -sequence "$inFile" -outseq "$outRefFile" 2>/dev/null

		# ref/*.msf
		# msfフォーマットに変換
		seqret -sequence "$outRefFile" -outseq msf::$tmpFile 2>/dev/null
		sed 's/[~]/./g' $tmpFile > "$outMsfFile"
		rm "$tmpFile"

		# inputdata/*.fasta
		if [ $useInDir -ne 0 ] ; then
			cp "${prefab_dir}/in/${i##*/}" "$outInputFile"
		else
			# gapをとりのぞき入力配列を作成(小文字を大文字へ変換)
			degapseq "$outRefFile" "$tmpFile" 2>/dev/null
			cat "$tmpFile" | awk '/^[>]/ { print $0;} /^[^>]/ { print toupper($0);}' > "$outInputFile"
			rm "$tmpFile"
		fi
	done
}

do_prefab1 (){
	do_prefab_common "$prefab1_dir" "refalns" $1
}

do_prefab2 (){
	do_prefab_common "$prefab2_dir" "refaln" $1
}

do_prefab3 (){
	do_prefab_common "$prefab3_dir" "ref" $1
}

do_prefab4 (){
	do_prefab_common "$prefab4_dir" "ref" $1
}

do_sabmark (){
	dirName="$1"
  outDir="$2"

	if [ ! -d $outDir ] ; then
		mkdir $outDir
	fi

	if [ ! -d $outDir/inputdata ] ; then
		mkdir $outDir/inputdata
	fi

	if [ ! -d $outDir/ref ] ; then
		mkdir $outDir/ref
	fi

  find "$dirName" -type f -name "*.fasta" | while read i ; do
    name=`basename $i`
  	name=${name%.*}
	  if [[ ! -z `echo "$name" | grep "group.*"` ]] ; then
			# 複数本が登録されている入力配列なので対象としない
		  echo CANCEL: $i
		else
			echo TARGET: $i
			seqret -auto $i $outDir/ref/${name}.ref_fasta
			seqret -auto $outDir/ref/${name}.ref_fasta msf:$outDir/ref/${name}.msf
			degapseq -auto $outDir/ref/${name}.ref_fasta $outDir/inputdata/${name}.fasta
		fi
  done

	# 圧縮
	databaseName=`basename $outDir`
	tar cvzf "${databaseName}.tar.gz" "$databaseName"
}

# do_sabmark (){
# 	cd ./$sabmark

# 	for j in twi sup ; do
# 		if [ ! -d ${sabmark}_${j} ] ; then
# 			mkdir ${sabmark}_${j}
# 		fi

# 		cd ${sabmark}_${j}

# 		if [ ! -d inputdata ] ; then
# 			mkdir inputdata
# 		fi

# 		for i in $(ls ../$j/*/reference/*.fasta); do
# 			echo TARGET: $i
# 			filename=${i##*/}
# 			echo OUTPUT: ${filename%.*}.ref_fasta
			
# 			seqret -sequence $i -outseq ${filename%.*}.ref_fasta
# 			degapseq ${filename%.*}.ref_fasta ./inputdata/${filename%.*}.fasta
# 			# msfフォーマットに変換
# 			seqret -sequence ${filename%.*}.ref_fasta -outseq msf::${filename%.*}.temp
# 			sed 's/[~]/./g' ${filename%.*}.temp > ${filename%.*}.msf
# 			rm ${filename%.*}.temp
# 		done

# 		cd ..
# 	done

# 	cd ..
# }

do_homstrad (){
	dirName="$1"
  outDir="$2"

	if [ ! -d $outDir ] ; then
		mkdir $outDir
	fi

	if [ ! -d $outDir/inputdata ] ; then
		mkdir $outDir/inputdata
	fi

	if [ ! -d $outDir/ref ] ; then
		mkdir $outDir/ref
	fi

	find "$dirName" -type f -name "*.ali" | while read i ; do
		name=`basename $i`
		name=${name%.*}
	  echo TARGET: $i
		# 入力データを作成
		sed -e '/C;/d' -e '/structure/d' -e 's/;/_/g' -e 's/[\*-]//g' -e 's/[\/]//g' $i > $outDir/inputdata/${name}.fasta
		# リファレンスアライメントを作成
		sed -e '/C;/d' -e '/structure/d' -e 's/;/_/g' -e 's/[\*]//g' -e 's/[\/]/-/g' $i > $outDir/ref/${name}.ref_fasta
		# msfフォーマットに変換
		seqret -auto -sequence $outDir/ref/${name}.ref_fasta -outseq msf::$outDir/ref/${name}.temp
		sed 's/[~]/./g' $outDir/ref/${name}.temp > $outDir/ref/${name}.msf
		rm $outDir/ref/${name}.temp
	done
}

#do_balibase
#do_homstrad "./homstrad_ali_only_2010_Mar_1" "./homstrad2010Mar1"
# [prefab] input ref.と同じ本数
#do_prefab1 0
#do_prefab2 0
#do_prefab3 0
#do_prefab4 0
# [prefab] input 複数本
#do_prefab1 1
do_prefab2 1
#do_prefab3 1
#do_prefab4 1
#do_sabmark "./SABmark1.63_superfamilies" "./SABmark1.63_sup"
#do_dali
