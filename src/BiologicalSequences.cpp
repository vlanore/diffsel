#include <string>
#include <vector>

extern const std::vector<char> AAset{'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P',
                              'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'a', 'c', 'd', 'e', 'f', 'g',
                              'h', 'i', 'k', 'l', 'm', 'n', 'p', 'q', 'r', 's', 't', 'v', 'w',
                              'y', '-', '?', '$', '.', 'B', 'Z', '*', 'X', 'x'};
extern const std::vector<char> DNAset{'A', 'C', 'G', 'T', 'a', 'c', 'g', 't', 'B', 'D', 'H', 'K', 'M',
                               'N', 'R', 'S', 'V', 'W', 'Y', 'b', 'd', 'h', 'k', 'm', 'n', 'r',
                               's', 'v', 'w', 'y', '-', '?', '$', '.', '*', 'X', 'x'};
extern const std::vector<char> RNAset{'A', 'C', 'G', 'U', 'a', 'c', 'g', 'u', 'B', 'D', 'H', 'K', 'M',
                               'N', 'R', 'S', 'V', 'W', 'Y', 'b', 'd', 'h', 'k', 'm', 'n', 'r',
                               's', 'v', 'w', 'y', '-', '?', '$', '.', '*', 'X', 'x'};

// amino acids
extern const std::string Alphabet = "Amino_Acids";
extern const std::vector<char> AminoAcids{'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M',
                                   'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '-'};
extern const std::vector<char> aminoacids{'a', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'k', 'l', 'm',
                                   'n', 'p', 'q', 'r', 's', 't', 'v', 'w', 'y', '-'};
extern const std::vector<char> DNAletters{'A', 'C', 'G', 'T'};
extern const std::vector<char> dnaletters{'a', 'c', 'g', 't'};
extern const std::vector<char> RNAletters{'A', 'C', 'G', 'U'};
extern const std::vector<char> rnaletters{'a', 'c', 'g', 'u'};

extern const std::vector<std::string> Codons{
    "TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG", "TAT", "TAC", "TAA", "TAG", "TGT",
    "TGC", "TGA", "TGG", "CTT", "CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG", "CAT", "CAC",
    "CAA", "CAG", "CGT", "CGC", "CGA", "CGG", "ATT", "ATC", "ATA", "ATG", "ACT", "ACC", "ACA",
    "ACG", "AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA", "AGG", "GTT", "GTC", "GTA", "GTG",
    "GCT", "GCC", "GCA", "GCG", "GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG"};
extern const int codonpos[][64]{{3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1,
                          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                         {3, 3, 3, 3, 1, 1, 1, 1, 0, 0, 0, 0, 2, 2, 2, 2, 3, 3, 3, 3, 1, 1,
                          1, 1, 0, 0, 0, 0, 2, 2, 2, 2, 3, 3, 3, 3, 1, 1, 1, 1, 0, 0, 0, 0,
                          2, 2, 2, 2, 3, 3, 3, 3, 1, 1, 1, 1, 0, 0, 0, 0, 2, 2, 2, 2},
                         {3, 1, 0, 2, 3, 1, 0, 2, 3, 1, 0, 2, 3, 1, 0, 2, 3, 1, 0, 2, 3, 1,
                          0, 2, 3, 1, 0, 2, 3, 1, 0, 2, 3, 1, 0, 2, 3, 1, 0, 2, 3, 1, 0, 2,
                          3, 1, 0, 2, 3, 1, 0, 2, 3, 1, 0, 2, 3, 1, 0, 2, 3, 1, 0, 2}};

extern const std::vector<int> UniStopCodons{10, 11, 14};
extern const std::vector<int> UniCodonCode{4,  4,  9,  9,  15, 15, 15, 15, 19, 19, -1, -1, 1,  1,  -1, 18,
                                    9,  9,  9,  9,  12, 12, 12, 12, 6,  6,  13, 13, 14, 14, 14, 14,
                                    7,  7,  7,  10, 16, 16, 16, 16, 11, 11, 8,  8,  15, 15, 14, 14,
                                    17, 17, 17, 17, 0,  0,  0,  0,  2,  2,  3,  3,  5,  5,  5,  5};
extern const std::vector<int> UniStopPos1{3, 3, 3};
extern const std::vector<int> UniStopPos2{0, 0, 2};
extern const std::vector<int> UniStopPos3{0, 2, 0};

extern const std::vector<int> MtInvStopCodons{10, 11};
extern const std::vector<int> MtInvCodonCode{
    4,  4,  9,  9,  15, 15, 15, 15, 19, 19, -1, -1, 1,  1,  18, 18, 9,  9,  9,  9,  12, 12,
    12, 12, 6,  6,  13, 13, 14, 14, 14, 14, 7,  7,  10, 10, 16, 16, 16, 16, 11, 11, 8,  8,
    15, 15, 15, 15, 17, 17, 17, 17, 0,  0,  0,  0,  2,  2,  3,  3,  5,  5,  5,  5};
extern const std::vector<int> MtInvStopPos1{3, 3};
extern const std::vector<int> MtInvStopPos2{0, 0};
extern const std::vector<int> MtInvStopPos3{0, 2};

// mammal mitochondrial genetic code
extern const std::vector<int> MtMamStopCodons{10, 11, 46, 47};
extern const std::vector<int> MtMamCodonCode{
    4,  4,  9,  9,  15, 15, 15, 15, 19, 19, -1, -1, 1,  1,  18, 18, 9,  9,  9,  9,  12, 12,
    12, 12, 6,  6,  13, 13, 14, 14, 14, 14, 7,  7,  10, 10, 16, 16, 16, 16, 11, 11, 8,  8,
    15, 15, -1, -1, 17, 17, 17, 17, 0,  0,  0,  0,  2,  2,  3,  3,  5,  5,  5,  5};
extern const std::vector<int> MtMamStopPos1{3, 3, 0, 0};
extern const std::vector<int> MtMamStopPos2{0, 0, 2, 2};
extern const std::vector<int> MtMamStopPos3{0, 2, 0, 2};
