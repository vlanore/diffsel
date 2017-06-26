# diffsel

[![Build Status](https://travis-ci.org/vlanore/diffsel.svg?branch=master)](https://travis-ci.org/vlanore/diffsel) 
[![codecov](https://codecov.io/gh/vlanore/diffsel/branch/master/graph/badge.svg)](https://codecov.io/gh/vlanore/diffsel)

## How to download and build

To get diffsel from a machine connected to the internet, type in a terminal:

```bash
git clone https://github.com/vlanore/diffsel.git
```

This should create a folder called `diffsel` (the diffsel root folder). You must go there before compiling diffsel:

```bash
cd diffsel
```

Then, to build diffsel simply run:

```bash
make
```

To check that everything ran well, look into the `_build`  folder to see if executables are present:

```bash
ls _build
```

This should display a list of files. If you see the following files, then diffsel was correctly built:
* `diffsel`
* `singleomega`
* `multigenesingleomega`
* `readdiffsel`

## How to format your data

Diffsel requires a tree and an alignment file to run.

Your alignment file must follow the Phylip format. In addition, the number of bases should be a multiple of 3 (as it will be interpreted as codons). For example, the following file would be a valid alignment for diffsel:

```phylip
8 6
S0      TCCTGA
S1      AATAGT
S2      GGATTT
S3      AATTCA
S4      CGAAGG
S5      AACGCT
S6      ACGAGT
S7      AATATT
```

Your tree file must follow the newick format. In addition, nodes should be labelled with "conditions". Conditions are integer labels comprised between `0` and `P-1` where `P` is the total number of conditions. The tree does not need to have branch lenghts. In addition, the leaves of the tree should have the same names as the sequences in your alignment file. For example, the following file would be a valid tree file for diffsel matching the alignment file above (with `P=2`):

```newick
((((((((S0:0,S1:1):0,(S2:0,S3:1):0):0,(S4:0,S5:1):0,(S6:0,S7:1):0):0):0,(S8:0,S9:1):0,(S10:0,S11:1):0):0,(S12:0,S13:1):0,(S14:0,S15:1):0):0):0):0
```

The `datà` folder in the diffsel root folder contains examples of data files usable with diffsel.

## How to run diffsel

Basic usage for diffsel is (from the diffsel root folder):

```
_build/diffsel -t <path to tree file> -d <path to alignment file> [options...] <name of run>
```
The name of the run is a string which is used to create files to store the results of the run. Note that the name of the run must be at the very end of the command.

Useful options include:
* `-ncond <number of conditions>̀` to specify the number of conditions `P` (see above);
* `-x <integer e> <integer u>` to tell diffsel to write to disc every `e` datapoints and to stop after `u` cycles of `e` datapoints.

For example, to run a diffsel chain called "myrun" on example data (from the `data` folder), with 3 conditions until 10k points have been written to disc (from the diffsel root folder):

```bash
_build/diffsel -t data/samhd1.tree -d data/samhd1.ali -ncond 3 -x 1 10000 myrun
```
