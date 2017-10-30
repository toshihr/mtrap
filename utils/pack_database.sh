#!/bin/sh
# ./pack_database.sh [folder name]
# folder nameに./はいらない．単純にフォルダ名のみ
name=`basename ${1}`
tar cvzf "${name}.tar.gz" "$1"
#cd $1
#tar czf "inputdata.tar.gz" "inputdata"
#tar czf "ref.tar.gz" "ref"
#cd - >/dev/null

# 中身確認は以下のコマンド．ファイル名のリストが出力される
# tar tf ***.tar.gz

# 指定ファイルの解凍は
# tar xf ***.tar.gz filename
