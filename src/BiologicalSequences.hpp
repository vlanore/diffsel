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
