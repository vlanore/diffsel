#ifndef BIOLOGICALSEQUENCES_H
#define BIOLOGICALSEQUENCES_H

#include <iostream>
#include <vector>

static const std::string Path = "";

enum DataType { DNA = 0, RNA = 1, Protein = 2, Other = 3 };

const int RYN = 2;
const std::vector<char> RYset{'R', 'Y'};

const std::vector<char> AAset{'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P',
                              'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'a', 'c', 'd', 'e', 'f', 'g',
                              'h', 'i', 'k', 'l', 'm', 'n', 'p', 'q', 'r', 's', 't', 'v', 'w',
                              'y', '-', '?', '$', '.', 'B', 'Z', '*', 'X', 'x'};
const int AAN = 49;
const int DNAN = 38;
const int RNAN = 38;
const std::vector<char> DNAset{'A', 'C', 'G', 'T', 'a', 'c', 'g', 't', 'B', 'D', 'H', 'K', 'M',
                               'N', 'R', 'S', 'V', 'W', 'Y', 'b', 'd', 'h', 'k', 'm', 'n', 'r',
                               's', 'v', 'w', 'y', '-', '?', '$', '.', '*', 'X', 'x'};
const std::vector<char> RNAset{'A', 'C', 'G', 'U', 'a', 'c', 'g', 'u', 'B', 'D', 'H', 'K', 'M',
                               'N', 'R', 'S', 'V', 'W', 'Y', 'b', 'd', 'h', 'k', 'm', 'n', 'r',
                               's', 'v', 'w', 'y', '-', '?', '$', '.', '*', 'X', 'x'};

// amino acids

const int precision = 10000;
const std::string Alphabet = "Amino_Acids";
const std::vector<char> AminoAcids{'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M',
                                   'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '-'};
const std::vector<char> aminoacids{'a', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'k', 'l', 'm',
                                   'n', 'p', 'q', 'r', 's', 't', 'v', 'w', 'y', '-'};
const std::vector<char> RYletters{'R', 'Y'};
const std::vector<char> DNAletters{'A', 'C', 'G', 'T'};
const std::vector<char> dnaletters{'a', 'c', 'g', 't'};
const std::vector<char> RNAletters{'A', 'C', 'G', 'U'};
const std::vector<char> rnaletters{'a', 'c', 'g', 'u'};

const std::vector<int> Dayhoff6Table{3, 5, 2, 2, 4, 3, 1, 0, 1, 0, 0, 2, 3, 2, 1, 3, 3, 0, 4, 4};
const std::vector<int> Dayhoff4Table{3, -1, 2, 2, 0, 3, 1, 0, 1, 0, 0, 2, 3, 2, 1, 3, 3, 0, 2, 2};

const int unknown = -1;

const int Naa = 20;
const int Nnuc = 4;
const int Ncodon = 64;
const std::vector<std::string> Codons{
    "TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG", "TAT", "TAC", "TAA", "TAG", "TGT",
    "TGC", "TGA", "TGG", "CTT", "CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG", "CAT", "CAC",
    "CAA", "CAG", "CGT", "CGC", "CGA", "CGG", "ATT", "ATC", "ATA", "ATG", "ACT", "ACC", "ACA",
    "ACG", "AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA", "AGG", "GTT", "GTC", "GTA", "GTG",
    "GCT", "GCC", "GCA", "GCG", "GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG"};
const int codonpos[][64]{{3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1,
                          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                         {3, 3, 3, 3, 1, 1, 1, 1, 0, 0, 0, 0, 2, 2, 2, 2, 3, 3, 3, 3, 1, 1,
                          1, 1, 0, 0, 0, 0, 2, 2, 2, 2, 3, 3, 3, 3, 1, 1, 1, 1, 0, 0, 0, 0,
                          2, 2, 2, 2, 3, 3, 3, 3, 1, 1, 1, 1, 0, 0, 0, 0, 2, 2, 2, 2},
                         {3, 1, 0, 2, 3, 1, 0, 2, 3, 1, 0, 2, 3, 1, 0, 2, 3, 1, 0, 2, 3, 1,
                          0, 2, 3, 1, 0, 2, 3, 1, 0, 2, 3, 1, 0, 2, 3, 1, 0, 2, 3, 1, 0, 2,
                          3, 1, 0, 2, 3, 1, 0, 2, 3, 1, 0, 2, 3, 1, 0, 2, 3, 1, 0, 2}};

enum GeneticCodeType { Universal = 0, MtMam = 1, MtInv = 2, MtProt = 3, MtEch = 4 };
// universal genetic code
// const std::string UniStopCodons[]  {"TAA","TAG","TGA"};
const int UniNStopCodons = 3;
const std::vector<int> UniStopCodons{10, 11, 14};
const std::vector<int> UniCodonCode{4,  4,  9,  9,  15, 15, 15, 15, 19, 19, -1, -1, 1,  1,  -1, 18,
                                    9,  9,  9,  9,  12, 12, 12, 12, 6,  6,  13, 13, 14, 14, 14, 14,
                                    7,  7,  7,  10, 16, 16, 16, 16, 11, 11, 8,  8,  15, 15, 14, 14,
                                    17, 17, 17, 17, 0,  0,  0,  0,  2,  2,  3,  3,  5,  5,  5,  5};
const std::vector<int> UniStopPos1{3, 3, 3};
const std::vector<int> UniStopPos2{0, 0, 2};
const std::vector<int> UniStopPos3{0, 2, 0};

const int MtInvNStopCodons = 2;
const std::vector<int> MtInvStopCodons{10, 11};
const std::vector<int> MtInvCodonCode{
    4,  4,  9,  9,  15, 15, 15, 15, 19, 19, -1, -1, 1,  1,  18, 18, 9,  9,  9,  9,  12, 12,
    12, 12, 6,  6,  13, 13, 14, 14, 14, 14, 7,  7,  10, 10, 16, 16, 16, 16, 11, 11, 8,  8,
    15, 15, 15, 15, 17, 17, 17, 17, 0,  0,  0,  0,  2,  2,  3,  3,  5,  5,  5,  5};
const std::vector<int> MtInvStopPos1{3, 3};
const std::vector<int> MtInvStopPos2{0, 0};
const std::vector<int> MtInvStopPos3{0, 2};

// mammal mitochondrial genetic code
const int MtMamNStopCodons = 4;
const std::vector<int> MtMamStopCodons{10, 11, 46, 47};
const std::vector<int> MtMamCodonCode{
    4,  4,  9,  9,  15, 15, 15, 15, 19, 19, -1, -1, 1,  1,  18, 18, 9,  9,  9,  9,  12, 12,
    12, 12, 6,  6,  13, 13, 14, 14, 14, 14, 7,  7,  10, 10, 16, 16, 16, 16, 11, 11, 8,  8,
    15, 15, -1, -1, 17, 17, 17, 17, 0,  0,  0,  0,  2,  2,  3,  3,  5,  5,  5,  5};
const std::vector<int> MtMamStopPos1{3, 3, 0, 0};
const std::vector<int> MtMamStopPos2{0, 0, 2, 2};
const std::vector<int> MtMamStopPos3{0, 2, 0, 2};

inline std::istream &operator>>(std::istream &is, GeneticCodeType &type) {
    std::string t;
    is >> t;
    if (t == "Universal") {
        type = Universal;
    } else if (t == "MtMam") {
        type = MtMam;
    } else if (t == "MtInv") {
        type = MtInv;
    } else if (t == "MtProt") {
        type = MtProt;
    } else if (t == "MtEch") {
        type = MtEch;
    } else {
        std::cerr << "error in std::istream genetic code type\n";
        std::cerr << type << '\n';
        exit(1);
    }
    return is;
}

inline std::ostream &operator<<(std::ostream &os, GeneticCodeType type) {
    if (type == Universal) {
        os << "Universal\n";
    } else if (type == MtMam) {
        os << "MtMam\n";
    } else if (type == MtInv) {
        os << "MtInv\n";
    } else if (type == MtProt) {
        os << "MtProt\n";
    } else if (type == MtEch) {
        os << "MtEch\n";
    } else {
        std::cerr << "error in std::ostream genetic code type\n";
        std::cerr << (int)type << '\n';
        exit(1);
    }
    return os;
}

#endif  // BIOLOGICALSEQUENCES_H
