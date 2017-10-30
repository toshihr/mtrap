#!/bin/sh
# ./ref/*.ref_fastaを元に
# ./ref/*.msfを生成し
# ./inputdata/に入力配列を生成する

ref_fasta_dir="$PWD/ref"
inputdata_dir="$PWD/inputdata"

cd $ref_fasta_dir

for i in $( ls *.ref_fasta ); do
	echo TARGET: $i

	seqret -auto -sequence ${i%.*}.ref_fasta -outseq msf::${i%.*}.msf
	degapseq -auto ${i%.*}.ref_fasta $inputdata_dir/${i%.*}.fasta
done

cd -


