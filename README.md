MTRAP is a sequence alignment program by a new measure based on transition
probability between two consecutive pairs of residues.

## Build & Install

  ```
  $ cmake .
  $ make
  $ sudo make install
  ```

Trouble shoots:
* if make failed, try `./clean.sh` to clean up unnecessary files.

## Usage

```
mtrap [OPTIONS] inputfile outputfile

===   OPTIONS  ===
= GENERAL =
-I                  [...]  add search path e.g. -I ~/mymatrix
-filelist           [...]
-i                  [...]  input file (FASTA format)
-o                  [...]  output file
-nosort                    the results are not sorted
-noestimation                 do not estimate a family
-cds      [GENETICCODES (comma separate)]  "protein" coding DNA sequence
             GENETICCODE 0: The Standard Code (transl_table=1)
                         1: The Vertebrate Mitochondrial Code (transl_table=2)
                         2: The Invertebrate Mitochondrial Code (transl_table=5)
                         3: The Bacterial, Archaeal and Plant Plastid Code (transl_table=11)

= GAP =
-go                 [-11]  gap open [score]
-ge                [-0.3]  gap extension [score]

= MULTIPLE ALIGNMENT =
-itr                 [10]  set the number of iterations
-c                    [2]  set the number of consistency transformations

= SUBSTITUTION MATRIX =
-m             [VTML200I]  substitution matrix e.g. CGONNET250, EBLOSUM62
-rm              [MRAMA1]  ramachandran matrix e.g. MRAMA1
-averageExtendedAA                 average extended amino acids (B,Z and X) substitutions

= TRANSITION QUANTITY MATRIX =
-tm       [SABmark1.63_sup_weighted.btq]  transition-quantity matrix
-e                [0.775]  epsilon, weight for transition quantity

= RAMACHANDRAN QUANTITY =
-gamma                [0]  gamma, weight for ramachandran matrix

= PARTITION FUNCTION POSTERIOR PROBABILITY =
-pf                 [0.7]  the degree of partition function
-pm            [VTML200I]  substitution matrix used for partition function
-pgo                [-22]  gap open [score] for partition function
-pge                 [-1]  gap extension [score] for partition function
-beta               [0.2]  beta, weight for partition function
-beta2              [1.5]  beta2, weight for a part of TQ of partition function

= GENERATING MODE =
-primarylibrary           [ID]  output T-Coffee primary library, ID: sequence identity (same as T-Coffee)

-count    [MODES FILELIST WEIGHTLIST]  count the frequency of transition
           MODES 0: normal 1:skip gap-gap site 2:set zero gap site
                 3: direct prob 4: without lower case 5: global TQ
                 f: symmetrize the frequencies not the tq
                 n: do not round the substitution matrix
-outbintm      [FILENAME]  output binary format transition-quantity matrix
```

## Version 2
MTRAP version 2 uses the following state-of-the-art approaches in addition to Transition-quantity.

* Partition function posterior probability
* Iterations & consistency transformations

These are introduced by [ProbCons](http://genome.cshlp.org/content/15/2/330.abstract)

As a result, profile matrix is calculated as below:
1. Calculate profile := pf * (MTRAP profile) + (1 - pf) * (posterior probability profile),
  where MTRAP profile is the paper version profile, i.e. e * Transition-quantity + (1 - e) * Score-matrix.

2. Smooth the profile by using iterations & consistency transformations.

## Documents
* [Making transition quantity](./doc/make_transition_quantity.md)
* [Details of the dynamic programming](./doc/proof.md)
* [Resources](./doc/resource.md)
* [Utility scripts](./doc/utils.md)

## Links
* Original site: [Ohya Lab.](http://www.rs.noda.tus.ac.jp/~ohya-m/) (Ohya Lab. was closed)
* Paper: Hara T., Sato K., Ohya M., "MTRAP: Pairwise sequence alignment algorithm by a new measure based on transition probability between two consecutive pairs of residues", BMC Bioinformatics 2010, 11:235. [link](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-235)

## License
MTRAP is released under the terms of the MIT license.

Note that the previous versions (before version 2) was licensed under BSD.

#### MTRAP uses the following open source libraries
* [Eigen](https://eigen.tuxfamily.org/) 3.1.2 (MPL2 License)
* [Snappy](https://github.com/google/snappy) 1.0.5 (BSD License)
