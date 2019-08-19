# diffsel

[![Build Status](https://travis-ci.org/vlanore/diffsel.svg?branch=master)](https://travis-ci.org/vlanore/diffsel) 
[![codecov](https://codecov.io/gh/vlanore/diffsel/branch/master/graph/badge.svg)](https://codecov.io/gh/vlanore/diffsel)
[![licence CeCILL](https://img.shields.io/badge/license-CeCILL--C-blue.svg)](http://www.cecill.info/licences.en.html)

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

## Check MCMC convergence

To ensure that analysis results are correct, one should check that the MCMC chain has converged.
You can do that either by using your MCMC convergence tool of choice using the `<name of run>.trace` files,
or by using the provided `tracecomp` program.

Tracecomp requires two runs with identical parameters.
You can invoke it by typing `_build/tracecomp -x <burn-in> <name of run 1> <name of run2>` from the diffsel root folder (assuming you ran your analysis from the diffsel root folder). If you don't know what burn-in to specify, 20% of the total number of iterations should be a decent default value.

For example, with the example analysis command above, assuming a second identical run named `myrun2` was done, the command could be (with a burn-in at 20% of 10000 iterations):

```bash
_build/tracecomp -x 2000 myrun1 myrun2
```

and it would output something like this (example output from a very small run):

```
setting upper limit to : 11
tmp.trace	 burnin : 5	sample size : 6
tmp2.trace	 burnin : 5	sample size : 6
name                effsize	rel_diff

#logprior           6		1.83775
lnL                 4		1.03644
length              5		1.25434
globent             5		0.632021
selvar1             6		1.11341
statent             5		1.03834
rrent               5		0.550811
diag                6		0
```

If all lines have a large effsize (>300) and a small rel_diff (<0.1) then your chain has converged well.
In the above example, the chain *has not* converged (but it only had 11 iterations).
For test runs, these conditions can be slightly relaxed while still giving decently reliable results (e.g., effsize>100 and rel_diff<0.15).

## How to run post-analysis after a diffsel run

When it runs, diffsel produces a series of files with raw data about the run. These files are named `<name of run>.<something>`, for example `myrun.run`, `myrun.param`, `myrun.chain`, etc...

The diffsel repository contains an executable, called `readdiffsel`, which can analyze these raw files to see if convergent sites were detected. To run this exectuable, run the command `_build/readdiffsel <name of run>` from the diffsel root folder (assuming you ran your analysis from the diffsel root folder). For example, if your analysis was run using the example command above, the command for post-analysis is:

```bash
_build/readdiffsel myrun
```

This command produces two files per condition after the first. For example, in the `myrun` example with 3 conditions, it should have produced the files `myrun_1.meandiffsel`, `myrun_1.signdiffsel`, `myrun_2.meandiffsel` and `myrun_2.signdiffsel`.

Among the produced files, the interesting one is the `*.signdiffsel` file for the convergent condition.
It contains the list of detected convergent effects.
It looks like this:

```tsv
0       17      0.909091        1.19692
1       3       0.0909091       -0.549043
1       10      0       -0.992412
1       12      0.909091        0.493234
1       17      1       1.29196
2       10      0.909091        1.19423
```

The columns are: site number, aminoacid number, "convergence probability", and effect strength.
Effects are already filtered to keep only those with probability above 0.9 or below 0.1, i.e., significant effects.

"Convergence probability" measures how consistently the corresponding aminoacid has a fitness that is better than the other aminoacids at the same position. See the original diffsel paper for more details.

If the `*.signdiffsel` file is empty, it means that no effects were detected.
You can lower the detection threshold in the `src/ReadDiffsel.cpp` file. Just change the lign that says `double cutoff = 0.9;` to another value, recompile and rerun `readdiffsel`.

## Another way to use diffsel: convergence pipeline

If you want to use diffsel but don't need to see the low-level MCMC-related files, you might prefer using our convergence detection pipeline.
This pipeline is capable of running diffsel, along with other convergence detection methods, and produces a standardized output for all methods (a score between 0 and 1 for each method and for each site).
This output might be simpler to understand and use than the post-analysis described above.

The pipeline and instructions are located at: https://gitlab.in2p3.fr/pveber/reviewphiltrans
