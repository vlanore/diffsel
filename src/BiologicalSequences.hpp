/*Copyright or © or Copr. Centre National de la Recherche Scientifique (CNRS) (2017-06-14).
Contributors:
* Nicolas LARTILLOT - nicolas.lartillot@univ-lyon1.fr

This software is a computer program whose purpose is to detect convergent evolution using Bayesian
phylogenetic codon models.

This software is governed by the CeCILL-C license under French law and abiding by the rules of
distribution of free software. You can use, modify and/ or redistribute the software under the terms
of the CeCILL-C license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and rights to copy, modify and redistribute
granted by the license, users are provided only with a limited warranty and the software's author,
the holder of the economic rights, and the successive licensors have only limited liability.

In this respect, the user's attention is drawn to the risks associated with loading, using,
modifying and/or developing or reproducing the software by the user in light of its specific status
of free software, that may mean that it is complicated to manipulate, and that also therefore means
that it is reserved for developers and experienced professionals having in-depth computer knowledge.
Users are therefore encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or data to be ensured and,
more generally, to use and operate it in the same conditions as regards security.

The fact that you are presently reading this means that you have had knowledge of the CeCILL-C
license and that you accept its terms.*/

#ifndef BIOLOGICALSEQUENCES_H
#define BIOLOGICALSEQUENCES_H

#include <iostream>
#include <vector>

const int AAN = 49;
const int DNAN = 38;
const int RNAN = 38;
const int unknown = -1;
const int Naa = 20;
const int Nnuc = 4;
const int Ncodon = 64;
const int precision = 10000;
const int MtMamNStopCodons = 4;
const int MtInvNStopCodons = 2;
const int UniNStopCodons = 3;

static const std::string Path = "";

extern const std::vector<char> AAset;
extern const std::vector<char> DNAset;
extern const std::vector<char> RNAset;

// amino acids
extern const std::string Alphabet;
extern const std::vector<char> AminoAcids;
extern const std::vector<char> aminoacids;
extern const std::vector<char> DNAletters;
extern const std::vector<char> dnaletters;
extern const std::vector<char> RNAletters;
extern const std::vector<char> rnaletters;

extern const std::vector<std::string> Codons;
extern const int codonpos[][64];

enum GeneticCodeType { Universal = 0, MtMam = 1, MtInv = 2, MtProt = 3, MtEch = 4 };

extern const std::vector<int> UniStopCodons;
extern const std::vector<int> UniCodonCode;
extern const std::vector<int> UniStopPos1;
extern const std::vector<int> UniStopPos2;
extern const std::vector<int> UniStopPos3;

extern const std::vector<int> MtInvStopCodons;
extern const std::vector<int> MtInvCodonCode;
extern const std::vector<int> MtInvStopPos1;
extern const std::vector<int> MtInvStopPos2;
extern const std::vector<int> MtInvStopPos3;

// mammal mitochondrial genetic code
extern const std::vector<int> MtMamStopCodons;
extern const std::vector<int> MtMamCodonCode;
extern const std::vector<int> MtMamStopPos1;
extern const std::vector<int> MtMamStopPos2;
extern const std::vector<int> MtMamStopPos3;

#endif  // BIOLOGICALSEQUENCES_H
