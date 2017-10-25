#!/bin/sh
bin="./src/mtrap"
#bin="mtrap"
#input="./test/sialidase_N"
target="BB11001"
$bin -go -11 -ge -0.3 -m CGONNET250 -tm SABmark1.63_sup_weighted.btq -e 0.8 -o ./result-${target}-GONNET250.fasta ./samples/${target}.fasta
