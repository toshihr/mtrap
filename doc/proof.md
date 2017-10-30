```
A: A1,A2,...,AnumRow
B: B1,B2,...,BnumCol
--- MTRAP ---
reference: MTRAP
PM[0] <-> D2 (from diagonal)
PM[1] <-> D1 (from upper) 論文ではPk=(ik,jk)と表現しているがx,y座標ではなくik行，jk列を表すと思えばよい
PM[2] <-> D3 (from left)
                                                                                                          The difference between MTRAP and gotoh's affine gap equation
D[r,c] = min {PM[0][r,c] = min { PM[0][r-1,c-1] + td((Ar-1,Bc-1) -> (Ar,Bc))   (r-2,c-2)D(r-1,c-1)D(r,c)  }} Dmn = min {Dm-1n-1 + score
                              { PM[1][r-1,c-1] + td((Ar-1,*)    -> (Ar,Bc))   (r-2,c-1)U(r-1,c-1)D(r,c)  }             This line is not defined
                              { PM[2][r-1,c-1] + td((*   ,Bc-1) -> (Ar,Bc))   (r-1,c-2)L(r-1,c-1)D(r,c)  }             This line is not defined
            {PM[1][r,c] = min { PM[0][r-1,c]   + td((Ar-1,Bc)   -> (Ar,* ))   (r-2,c-1)D(r-1,c  )U(r,c)  }            {Pmn = min { Dm-1n + Open
                              { PM[1][r-1,c]   + td((Ar-1,*)    -> (Ar,* ))   (r-2,c  )U(r-1,c  )U(r,c)  }                       { Pm-1n + Ext
                              { PM[2][r-1,c]   + td((*   ,Bc)   -> (Ar,* ))   (r-1,c-1)L(r-1,c  )U(r,c)  }                       { transition LU is not defined
            {PM[2][r,c] = min { PM[0][r,c-1]   + td((Ar,Bc-1)   -> (* ,Bc))   (r-1,c-2)D(r  ,c-1)L(r,c)  }            {Qmn = min { Dmn-1 + Open
                              { PM[1][r,c-1]   + td((Ar,*)      -> (* ,Bc))   (r  ,c-2)U(r  ,c-1)L(r,c)  }                       { transition UL is not defined
                              { PM[2][r,c-1]   + td((* ,Bc-1)   -> (* ,Bc))   (r-1,c-1)L(r  ,c-1)L(r,c)  }                       { Qmn-1 + Ext
            for r=1,...,numRow
                c=1,...,numCol
td((Ar-1,Bc-1) -> (Ar,Bc)) := tq + s(Ar,Bc)
td((Ar-1,Bc  ) -> (Ar,* )) := tq + gap open
td((Ar  ,Bc-1) -> (* ,Bc)) := tq + gap open
td((Ar-1,*  ) -> (Ar,* )) := tq + gap ext
td((*  ,Bc-1) -> (* ,Bc)) := tq + gap ext
where tq = 0 when r = 1 or c = 1
--- boundary conditions ---
                 PM[0][0,0] =   0
                 PM[1][0,0] =   +inf (-inf for score mode)
                 PM[2][0,0] =   +inf (-inf for score mode)
                 PM[0][1,0] =   +inf (-inf for score mode)
                 PM[1][1,0] =   PM[0][0,0]     + td((@   ,*)    -> (A1,* ))
                 PM[2][1,0] =   +inf (-inf for score mode)
                 PM[0][r,0] =   +inf (-inf for score mode)
                 PM[1][r,0] =   PM[1][r-1,0]   + td((Ar-1,*)    -> (Ar,* ))
                 PM[2][r,0] =   +inf (-inf for score mode)
                   for r=2,...,numRow
                 PM[0][0,1] =   +inf (-inf for score mode)
                 PM[1][0,1] =   +inf (-inf for score mode)
                 PM[2][0,1] =   PM[0][0,0]     + td((* ,@   )   -> (* ,B1))
                 PM[0][0,c] =   +inf (-inf for score mode)
                 PM[1][0,c] =   +inf (-inf for score mode)
                 PM[2][0,c] =   PM[2][0,c-1]   + td((* ,Bc-1)   -> (* ,Bc))
                   for c=1,...,numCol
--- MTRAP流定義が，SmithWaterman流再帰によるdistance matrixの定義を含んだものであるかの証明を epsilon = 0 の場合を用いて行う---
STEP1: tq -> s
D[r,c] = min {PM[0][r,c] = min { PM[0][r-1,c-1] + s(Ar,Bc)   (r-2,c-2)D(r-1,c-1)D(r,c)  }} Dmn = min {Dm-1n-1 + score
                              { PM[1][r-1,c-1] + s(Ar,Bc)   (r-2,c-1)U(r-1,c-1)D(r,c)  }             This line is not defined
                              { PM[2][r-1,c-1] + s(Ar,Bc)   (r-1,c-2)L(r-1,c-1)D(r,c)  }             This line is not defined
            {PM[1][r,c] = min { PM[0][r-1,c]   + GOPEN      (r-2,c-1)D(r-1,c  )U(r,c)  }            {Pmn = min { Dm-1n + Open
                              { PM[1][r-1,c]   + GEXT       (r-2,c  )U(r-1,c  )U(r,c)  }                       { Pm-1n + Ext
                              { PM[2][r-1,c]   + GOPEN      (r-1,c-1)L(r-1,c  )U(r,c)  }                       { transition LU is not defined
            {PM[2][r,c] = min { PM[0][r,c-1]   + GOPEN      (r-1,c-2)D(r  ,c-1)L(r,c)  }            {Qmn = min { Dmn-1 + Open
                              { PM[1][r,c-1]   + GOPEN      (r  ,c-2)U(r  ,c-1)L(r,c)  }                       { transition UL is not defined
                              { PM[2][r,c-1]   + GEXT       (r-1,c-1)L(r  ,c-1)L(r,c)  }                       { Qmn-1 + Ext
STEP2: place out s
D[r,c] = min {PM[0][r,c] = min { PM[0][r-1,c-1] } + s(Ar,Bc)
                              { PM[1][r-1,c-1] }           
                              { PM[2][r-1,c-1] }           
            {PM[1][r,c] = min { PM[0][r-1,c]   + GOPEN     
                              { PM[1][r-1,c]   + GEXT      
                              { PM[2][r-1,c]   + GOPEN     
            {PM[2][r,c] = min { PM[0][r,c-1]   + GOPEN     
                              { PM[1][r,c-1]   + GOPEN     
                              { PM[2][r,c-1]   + GEXT      
STEP3: replace D[r-1,c-1]
D[r,c] = min {PM[0][r,c] = D[r-1,c-1] + s(Ar,Bc)
            {PM[1][r,c] = min { PM[0][r-1,c]   + GOPEN
                              { PM[1][r-1,c]   + GEXT
                              { PM[2][r-1,c]   + GOPEN
            {PM[2][r,c] = min { PM[0][r,c-1]   + GOPEN
                              { PM[1][r,c-1]   + GOPEN
                              { PM[2][r,c-1]   + GEXT
STEP4: omit PM[0]
D[r,c] = min {D[r-1,c-1] + s(Ar,Bc)
            {PM[1][r,c] = min { D[r-2,c-1] + s(Ar-1,Bc)   + GOPEN  ここで D[r-2,c-1] + s(Ar-1,Bc) <==> D[r-1,c]=min{D[r-2,c-1]+s(Ar-1,Bc), PM[1][r-1,c], PM[2][r-1,c]}
                              { PM[1][r-1,c]   + GEXT                               gotohではこれが等しい場合を考えている
                              { PM[2][r-1,c]   + GOPEN
            {PM[2][r,c] = min { D[r-1,c-2] + s(Ar,Bc-1)   + GOPEN  ここで D[r-1,c-2] + s(Ar,Bc-1) <==> D[r,c-1]=min{D[r-1,c-2]+s(Ar,Bc-1), PM[1][r,c-1], PM[2][r,c-1]}
                              { PM[1][r,c-1]   + GOPEN
                              { PM[2][r,c-1]   + GEXT
STEP5: ここで，おきかえをした場合を考えてみる
D[r,c] = min {D[r-1,c-1] + s(Ar,Bc)
            {PM[1][r,c] = min { min{ (1) D[r-2,c-1] + s(Ar-1,Bc), (2) PM[1][r-1,c], (3) PM[2][r-1,c]}  + GOPEN
                              { (4) PM[1][r-1,c]   + GEXT
                              { (5) PM[2][r-1,c]   + GOPEN
            {PM[2][r,c] = min { min{D[r-1,c-2] + s(Ar,Bc-1), PM[1][r,c-1], PM[2][r,c-1]}  + GOPEN
                              { PM[1][r,c-1]   + GOPEN
                              { PM[2][r,c-1]   + GEXT
ここで，min{(1),(2),(3)} = (2) となってしまう場合を考えてみると，GOPEN > GEXT (convex gap function) より min{(1),(2),(3)} + GOPEN > (4) となり，この場合は考えなくて良い
ここで，min{(1),(2),(3)} = (3) となってしまう場合を考えてみると，これは (5) と同値であり，この場合は考えなくてよい．逆に言えば(5)は定義に含めなくてよい．
結局，min{(1),(2),(3)} = (1) と置き換えをしてしまって問題はない．また，置き換える場合において (5) は最初の定義からいらないことが分かる．
STEP6: ここまでを整理すると
D[r,c] = min {D[r-1,c-1] + s(Ar,Bc)
            {PM[1][r,c] = min { D[r-1,c]  + GOPEN
                              { PM[1][r-1,c]   + GEXT
            {PM[2][r,c] = min { D[r,c-1]  + GOPEN
                              { PM[2][r,c-1]   + GEXT
よって，gotoh's affine gap algorithmをepsilon = 0 の場合として含むことが示せた．
※一見定義から(5)をはぶくとクロスギャップはゆるしていないように見えるが，STEP5にてmin{(1),(2),(3)} = (3)となる場合はクロスギャップとなり，実はクロスギャップは認めている
※STEP3における式はmsaprobsの式と対応する．ただし，msaprobsの式はクロスギャップの項がたりない．D[]でリプレースを行った段階で省略はできる．でも論文ではリプレースを行なっていない．
 hasegawa論文の式(32)ではリプレースを行なっている．そのため問題ない．
--- gotoh's affine gap ---
reference: Gotoh, An Improved Algorithm for Matching Biological Sequences, J.Mol.Biol., 162, 705-708 (1982)
PM[0] <-> D (from diagonal)
PM[1] <-> P (from upper)
PM[2] <-> Q (from left)
(r,c) <-> (m,n)
PM[0][r,c] = min( PM[0][r-1,c-1] + d(Ar, Bc), PM[1][r,c], PM[2][r,c] ) for r=1,2,...,numRow, c=1,2,...,NumCol
PM[1][r,c] = min( PM[0][r-1,c] + GAP_OPEN, PM[1][r-1,c] + GAP_EXT ) for r=1,2,...,numRow, c=1,2,...,NumCol
PM[2][r,c] = min( PM[0][r,c-1] + GAP_OPEN, PM[2][r,c-1] + GAP_EXT ) for r=1,2,...,numRow, c=1,2,...,NumCol
--- boundary conditions for global alignment ---
PM[0][r,0] = PM[1][r,0] = GAP_OPEN + GAP_EXT*(r-1) for r=1,...,numRow
PM[0][0,c] = PM[2][0,c] = GAP_OPEN + GAP_EXT*(r-1) for c=1,...,numCol
PM[2][r,0] = undefined for r=1,...,numRow
PM[1][0,c] = undefined for c=1,...,numCol
PM[0][0,0] = PM[1][0,0] = PM[2][0,0] = undefined
```
