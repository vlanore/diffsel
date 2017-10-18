/*Copyright or © or Copr. Centre National de la Recherche Scientifique (CNRS) (2017-06-14).
Contributors:
* Nicolas LARTILLOT - nicolas.lartillot@univ-lyon1.fr
* Vincent LANORE - vincent.lanore@univ-lyon1.fr

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

#include <string>
#include <vector>

extern const std::vector<char> AAset{
    'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T',
    'V', 'W', 'Y', 'a', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'k', 'l', 'm', 'n', 'p', 'q',
    'r', 's', 't', 'v', 'w', 'y', '-', '?', '$', '.', 'B', 'Z', '*', 'X', 'x'};
extern const std::vector<char> DNAset{
    'A', 'C', 'G', 'T', 'a', 'c', 'g', 't', 'B', 'D', 'H', 'K', 'M', 'N', 'R', 'S', 'V', 'W', 'Y',
    'b', 'd', 'h', 'k', 'm', 'n', 'r', 's', 'v', 'w', 'y', '-', '?', '$', '.', '*', 'X', 'x'};
extern const std::vector<char> RNAset{
    'A', 'C', 'G', 'U', 'a', 'c', 'g', 'u', 'B', 'D', 'H', 'K', 'M', 'N', 'R', 'S', 'V', 'W', 'Y',
    'b', 'd', 'h', 'k', 'm', 'n', 'r', 's', 'v', 'w', 'y', '-', '?', '$', '.', '*', 'X', 'x'};

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
extern const std::vector<int> UniCodonCode{
    4,  4,  9,  9,  15, 15, 15, 15, 19, 19, -1, -1, 1, 1,  -1, 18, 9,  9,  9,  9,  12, 12,
    12, 12, 6,  6,  13, 13, 14, 14, 14, 14, 7,  7,  7, 10, 16, 16, 16, 16, 11, 11, 8,  8,
    15, 15, 14, 14, 17, 17, 17, 17, 0,  0,  0,  0,  2, 2,  3,  3,  5,  5,  5,  5};
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
