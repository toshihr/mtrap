#!/usr/bin/ruby
# -*- coding: utf-8 -*-
require 'bio'
serv = Bio::DDBJ::XML::GetEntry.new

# ３文字表記から１文字表記へのハッシュ
A3toA1 = {}
Bio::AminoAcid::NAMES.each_pair {|key, value| A3toA1[value.upcase] = key}

def getSeqFromPDB(pdb)
  s = ""
  pdb.each_residue {|residue| s += A3toA1[residue.resName]}
  s
end

def isUpperCase(c)
  if "A"[0] <= c && c <= "Z"[0] then
    true
  else
    false
  end
end

# PDB SEQUENCES をすべてメモリへ読み込み(40MB程度)
#puts "Loading..."
#SEQUENCES = {}
#open("/Users/keru/research/databases/PDB/pdb.fasta") do |file|
#  header = file.gets
#  SEQUENCES[header] = file.gets
#end
#puts "done."

# DALIIDからPDBID+CHAINIDを抽出
def getPDBIDfromDALIID(s)
  if s.size == 4 then
    return s.upcase
  else
    # chainidに対しては大文字小文字変換を行わない
    return s[0..s.size-2].upcase #+ s[s.size-1..s.size-1]
  end
end

# ID seq1 seq2 zscore
# seqはPDBID+(大文字１文字)
TABLE_ALIGN = []
TABLE_ALIGN.push(["3916304","1tu9A","1jpbA","16.7"])
#TABLE_ALIGN.push(["2","101m","102m","30.6"])
#open("../dalidb_dccp") do |file|
#  while l = file.gets
#    TABLE_ALIGN.push(l.split(" "))
#  end
#end

# ID -> [[start1,start2,blocklength]]
# 複数のセグメントに対応
TABLE_SEGMENT = {}
TABLE_SEGMENT["3916304"] = []
TABLE_SEGMENT["3916304"].push(["1","5","16"])
TABLE_SEGMENT["3916304"].push(["18","23","24"])
TABLE_SEGMENT["3916304"].push(["42","60","23"])
TABLE_SEGMENT["3916304"].push(["66","83","12"])
TABLE_SEGMENT["3916304"].push(["79","95","27"])
TABLE_SEGMENT["3916304"].push(["106","123","26"])
#open("../dalidb_segments") do |file|
#  while l = file.gets
#    v = l.split(" ")
#    TABLE_ALIGN[v[0]].push(v[1..3])
#  end
#end


TABLE_ALIGN.each do |v|
  align_id = v[0]

  # 配列を取得 (TODO: エラー処理を追加すべき)
  id1 = getPDBIDfromDALIID(v[1])
  id2 = getPDBIDfromDALIID(v[2])
  pdb1 = Bio::PDB.new(serv.getPDBEntry(id1))
  pdb2 = Bio::PDB.new(serv.getPDBEntry(id2))
  seq1 = getSeqFromPDB(pdb1)
  seq2 = getSeqFromPDB(pdb2)

  puts pdb1.entry_id + ":" + seq1
  puts pdb2.entry_id + ":" + seq2
  puts ""

  # アライメントを構築
  len1 = seq1.size
  len2 = seq2.size
  align1 = ""
  align2 = ""

  # C++イテレータと同じ．最後のインデックス＋１を指す
  preSegEnd1 = 1
  preSegEnd2 = 1
  # 各セグメント
  TABLE_SEGMENT[align_id].each do |seg|
    # align1: preSegEnd1 から seg[0]-1 までギャップあるいはブロック外のアミノ酸を挿入
    # align2: preSegEnd2 から seg[1]-1 までギャップあるいはブロック外のアミノ酸を挿入
    nonSegLen1 = seg[0].to_i - preSegEnd1
    nonSegLen2 = seg[1].to_i - preSegEnd2
    # nonBlockでの対応はわからないため左詰でアミノ酸を出力してしまう
    # nonBlockなアミノ酸を出力
#    if preSegEnd1 > 1
      align1 += seq1[preSegEnd1-1, nonSegLen1].downcase
#    end
#    if preSegEnd2 > 1
      align2 += seq2[preSegEnd2-1, nonSegLen2].downcase
#    end
    # アライメント長を調整 gapLen > 0 なら gapLen分 align2にギャップを挿入． < 0 ならその逆
    gapLen = align1.size - align2.size
    if gapLen > 0
      align2 += "-" * gapLen
    elsif gapLen < 0
      align1 += "-" * (-gapLen)
    end
    # Blockを出力
    align1 += seq1[seg[0].to_i-1,seg[2].to_i]
    align2 += seq2[seg[1].to_i-1,seg[2].to_i]
    # 更新
    preSegEnd1 = seg[0].to_i + seg[2].to_i
    preSegEnd2 = seg[1].to_i + seg[2].to_i
  end

  # 最後に長さがそろっていなかったらnonBlockが残っていると考え，余ったアミノ酸，ギャップを挿入
  if preSegEnd1 <= seq1.size
    align1 += seq1[preSegEnd1-1, seq1.size - preSegEnd1 + 1].downcase
  end
  if preSegEnd2 <= seq2.size
    align2 += seq2[preSegEnd2-1, seq2.size - preSegEnd2 + 1].downcase
  end

  gapLen = align1.size - align2.size
  if gapLen > 0
    align2 += "-" * gapLen
  elsif gapLen < 0
    align1 += "-" * (-gapLen)
  end

  puts pdb1.entry_id + ":" + align1 + ":" + pdb2.entry_id + ":" + align2
end
