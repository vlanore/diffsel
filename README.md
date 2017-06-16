# diffsel

[![Build Status](https://travis-ci.org/vlanore/diffsel.svg?branch=master)](https://travis-ci.org/vlanore/diffsel) 
[![codecov](https://codecov.io/gh/vlanore/diffsel/branch/master/graph/badge.svg)](https://codecov.io/gh/vlanore/diffsel)

## How to download and build

To get diffsel from a machine connected to the internet, type in a terminal:

```bash
git clone https://github.com/vlanore/diffsel.git
```

This should create a folder called `diffsel` (the diffsel root folder). This folder contains a `src` subfolder in which you must go before compiling diffsel:

```bash
cd diffsel/src
```

Then, to build diffsel simply run:

```bash
make
```

To check that everything ran well, you can go back to the diffsel root folder and see if executables were correctly created:
```bash
cd ..
ls
```

This should display a list of files. If you see the following files, then diffsel was correctly built:
* `flatdiffsel_bin`
* `flatgeneglobom_bin`
* `flatglobom_bin`
* `readflatdiffsel_bin`

## How to format your data

Diffsel requires a tree and an alignment file to run.

Your alignment file must follow the Phylip format. In addition, the number of bases should be a multiple of 3 (as it will be interpreted as codons). For example, the following file would be a valid alignment for diffsel:

```
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
