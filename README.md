Copyright (C) 2012-2020, Marx Gomes van der Linden

# HmmPred 1.1

HmmPred implements an algorithm based on a Hidden Markov Model for
prediction of structural features of proteins from amino acid
sequences. It is designed as generalization of a method developed by
Crooks and Brenner (2004) for secondary structure prediction, with
the main purpose of performing prediction of burial levels of
residues in proteins.

## Prerequisites

HmmPred was compiled with GCC (g++) and tested in a GNU/Linux
environment, but it should probably work just fine in any common
operating system with a standard C++ compiler. To compile HmmPred,
you will need the following libraries:

- The "program_options" library from Boost
- The GNU Scientific Library (GSL)
- The CppTest unit testing framework

All of them are easily available as free software.

## Compiling

The easiest way to compile HmmPred is simply to run:

    g++ src/*.cpp -oHmmPred -O3 -lboost_program_options -lgsl -lgslcblas -lcpptest


## The Datasets

The distribution comes with several datasets, each with training and
evaluation files. In each group, the alphabet of the feature to be
predicted has a different size and meaning:

**`burial<CA|CB>/<N>layers`**

  Burial layer of alpha (CA) or beta (CB) carbons, split into N
  layers.

  The alphabet is "01" for 2 layers, "012" for 3 layers, etc..
  '0' is the inner layer.

  When using these datasets, add options `-slayers -lN` to HmmPred,
  where N is the number of layers.

**`burial<CA|CB>/+secstruct`**

  Burial layer for 2 layers, plus information about secondary
  structure.

  The alphabet is "abcdef" (use `-schars -l6`) where:

-  a: Inner layer, sheet
-  b: Inner layer, helix
-  c: Inner layer, loop
-  d: Outer layer, sheet
-  e: Outer layer, helix
-  f: Outer layer, loop

  If you want the output to contain only information about burials,
  add option `-r000111`. If you want to predict only secondary
  structures, add `-rEHLEHL`.

**`burial<CA|CB>/+secstruct9`**

  Burial layer for 3 layers, plus information about secondary
  structure.

  The alphabet is "abcdefghi" (use `-schars -l9`) where:

-  a:  Inner layer, sheet
-  b:  Inner layer, helix
-  c:  Inner layer,  loop
-  d: Middle layer, sheet
-  e: Middle layer, helix
-  f: Middle layer, loop
-  g:  Outer layer, sheet
-  h:  Outer layer, helix
-  i:  Outer layer, loop

  If you want the output to contain only information about burials,
  add option `-r000111222`. If you want to predict only secondary
  structures, add `-rEHLEHLEHL`.

**`burial<CAdir|CBdir>/<N>layers`**

  Burial layer of the alpha carbon, plus the orientation of the
  beta carbon (CAdir) or the other way around (CBdir).

  The alphabet for 2 layers is "0^v1" (`-scadir`), where:

-  0: CA in inner layer, CB directed inwards
-  ^: CA in inner layer, CB directed outwards
-  v: CA in outer layer, CB directed inwards
-  1: CA in outer layer, CB directed outwards

  To output only the position of CA, add `-r0011`.

 **Others:**

  For 3-6 layers, the alphabet is in the form "abcd(...)" (`-schars
  -lN`, where N is twice the number of layers), following the same
  pattern. For example, for 4 layers:

-  a:  CA in inner layer, CB directed inwards
-  b:  CA in inner layer, CB directed outwards
-  c: CA in second layer, CB directed inwards
-  d: CA in second layer, CB directed outwards
-  e:  CA in third layer, CB directed inwards
- f:  CA in third layer, CB directed outwards
-  g:  CA in outer layer, CB directed inwards
-  h:  CA in outer layer, CB directed outwards

  To output only the position of CA, use:
-  3 layers: "-schars -l6  -r001122"
-  4 layers: "-schars -l8  -r00112233"
-  5 layers: "-schars -l10 -r0011223344"
-  6 layers: "-schars -l12 -r001122334455"


  The same options are valid for beta carbon ouput in the CBdir
  directory.

**`secstruct-stride`**

  Secondary structure, using the same database employed by Crooks
  and Brenner (2004). The source alphabet is "HGIEBbTSC _-L" but
  it should be reduced to "EHL" (sheet, helix, loop) by using
  either the `-sehl` or `-sck` options.

In addition to the above, all alphabets also include an 'X' symbol
to represent unknown values.

## Examples

Run HmmPred without any parameters to see the usage instructions.
You have to provide a training file (with known associations between
primary structure and whatever you are trying to predict) and an
evaluation file (which is just a list of sequences to predict), as
well as filename to output the results. Other options describe what
kind of data is to be expected in the input files, as well as the
algorithm parameters.

Predicts alpha carbon burial level into 3 layers, with a
  window size of 4:

    ./HmmPred datasets/burial-CA/3layers/lista.training \
            datasets/burial-CA/3layers/lista.eval \
            output_file \
            -slayers -l3 -w4


Same as above, but performs full bootstraping with 30
         replicas and provides detailed statistics.
         (Note the `.eval+` file, instead of `.eval`):

    ./HmmPred datasets/burial-CA/3layers/lista.training \
            datasets/burial-CA/3layers/lista.eval+ \
            output_file \
            -slayers -l3 -w4 --stats -n30


Predicts alpha carbon burial and the relative orientation
           of the beta carbon, with a window size of 5:       

    ./HmmPred datasets/burial-CAdir/2layers/lista.training \
            datasets/burial-CAdir/2layers/lista.eval \
            output_file \
            -scadir -w5

Same as above, but represents only the position of CA in
           the output:

    ./HmmPred datasets/burial-CAdir/2layers/lista.training \
            datasets/burial-CAdir/2layers/lista.eval \
            output_file \
            -scadir -w5 -r0011


## References

HmmPred is used in the research published on the following papers:

- ROCHA, J. R. ; VAN DER LINDEN, M. G. ; FERREIRA, D. C. ; AZEVEDO, P. H. ; PEREIRA DE ARAUJO, A. F. . Information-theoretic analysis and prediction of protein atomic burials: on the search for an informational intermediate between sequence and structure. Bioinformatics (Oxford. Print), v. 28, p. 2755-2762, 2012. 

- VAN DER LINDEN, MARX GOMES; FERREIRA, DIOGO CÉSAR ; DE OLIVEIRA, LEANDRO CRISTANTE ; ONUCHIC, JOSÉ N. ; DE ARAÚJO, ANTÔNIO F. PEREIRA . Ab initio protein folding simulations using atomic burials as informational intermediates between sequence and structure. Proteins (Print), v. 82, p. n/a-n/a, 2013. 

- FERREIRA, DIOGO C. ; VAN DER LINDEN, MARX G. ; DE OLIVEIRA, LEANDRO C. ; ONUCHIC, JOSÉ N. ; PEREIRA DE ARAÚJO, ANTÔNIO F. . Information and redundancy in the burial folding code of globular proteins within a wide range of shapes and sizes. Proteins (Print), v. 84, p. 515-531, 2016. 

The algorithm implemented in HmmPred is based on:

- Crooks Ge, Brenner Se. 2004. Protein Secondary Structure: Entropy,
Correlations, And Prediction. Bioinformatics 20:1603-1611
