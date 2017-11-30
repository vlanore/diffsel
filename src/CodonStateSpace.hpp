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

#ifndef CODONSTATESPACE_H
#define CODONSTATESPACE_H

#include <map>
#include "StateSpace.hpp"

class CodonStateSpace : public StateSpace {
    using Codon = int;
    using Base = int;
    using AminoAcid = int;

  public:
    static const int Npos = 3;

    // by default, codons always exclude stops
    // if a method takes or returns a codon stops INCLUDED, then this is made
    // explicit in the method's name
    explicit CodonStateSpace(GeneticCodeType type);

    ~CodonStateSpace() override;

    // -----
    // generic methods
    // exist for any state space, and should have a consistent meaning throughout

    int GetNstate() const override { return Nstate; }

    // gets a three letter code, returns codon (if stop exits with error message)
    Codon GetState(std::string word) const override;

    // give a codon (stops excluded), returns a three letter code
    std::string GetState(Codon codon) const override;

    // -----
    // codon specific methods

    const DNAStateSpace *GetDNAStateSpace() const { return nucstatespace; }

    const ProteinStateSpace *GetProteinStateSpace() const { return protstatespace; }

    // returns a codon based on three letters
    // returns -1 (== unknown) if at least one of the positions is unknown
    // if stop exits with error message...
    Codon GetCodonFromDNA(Base pos1, Base pos2, Base pos3) const;

    // 2 codons excluding stops are compared
    // method returns -1 if identical
    // returns 3 is codons differ at more than one position
    // otherwise, returns the position at which codons differ (i.e. returns 0,1 or
    // 2 if the codons
    // differ at position 1,2 or 3)
    int GetDifferingPosition(Codon i, Codon j) const {
        // identical
        if ((GetCodonPosition(0, i) == GetCodonPosition(0, j)) &&
            (GetCodonPosition(1, i) == GetCodonPosition(1, j)) &&
            (GetCodonPosition(2, i) == GetCodonPosition(2, j))) {
            return -1;
        }
        if (GetCodonPosition(0, i) != GetCodonPosition(0, j)) {
            if ((GetCodonPosition(1, i) == GetCodonPosition(1, j)) &&
                (GetCodonPosition(2, i) == GetCodonPosition(2, j))) {
                return 0;
            }
            return 3;
        }
        if (GetCodonPosition(1, i) != GetCodonPosition(1, j)) {
            if ((GetCodonPosition(0, i) == GetCodonPosition(0, j)) &&
                (GetCodonPosition(2, i) == GetCodonPosition(2, j))) {
                return 1;
            }
            return 3;
        }
        if (GetCodonPosition(2, i) != GetCodonPosition(2, j)) {
            if ((GetCodonPosition(1, i) == GetCodonPosition(1, j)) &&
                (GetCodonPosition(0, i) == GetCodonPosition(0, j))) {
                return 2;
            }
            return 3;
        }
        return 3;
    }

    // return the integer encoding for the base at requested position
    // stops excluded
    Base GetCodonPosition(int pos, Codon codon) const {
        if ((pos < 0) || (pos >= Npos)) {
            std::cerr << "GetCodonPosition: pos out of bound\n";
            std::cerr << pos << '\n';
            exit(1);
        }
        if (codon == -1) {
            return -1;
        }
        if ((codon < 0) || (codon >= Nstate)) {
            std::cerr << "GetCodonPosition: codon out of bound\n";
            std::cerr << codon << '\n';
            exit(1);
        }
        return CodonPos[pos][codon];
    }

    int IsNonCTNearest(int a, int b) const;

    // translation stops excluded
    int Translation(Codon codon) const { return CodonCode[codon]; }

    // stops excluded
    bool Synonymous(Codon codon1, Codon codon2) const {
        return (CodonCode[codon1] == CodonCode[codon2]);
    }

    // returns -1 if stop codon
    // otherwise returns integer in [0,19] standing for an amino-acid (one letter
    // code, alphabetical
    // order)
    AminoAcid TranslationWithStops(Codon codon) const { return CodonCodeWithStops[codon]; }

    bool CheckStop(int pos1, int pos2, int pos3) const;

    int GetDegeneracy(int codon) const;

    int GetNstop() const { return Nstop; }

    const int *GetStopPos1() const { return StopPos1; }

    const int *GetStopPos2() const { return StopPos2; }

    const int *GetStopPos3() const { return StopPos3; }

  private:
    void MakeDegeneracyMap() const;

    GeneticCodeType code;
    const DNAStateSpace *nucstatespace;
    const ProteinStateSpace *protstatespace;
    // number of codons, not including stops (61 in general)
    int Nstate;

    // and array of size Ncodon = 64
    // whose entries are between -1 and 19
    // -1 : stop codon
    // 0..19 : amino acid encoded (1 letter code, alphabetical order)
    int *CodonCodeWithStops;
    int *CodonCode;
    int **CodonPos;
    int *StopCodons;
    int Nstop;
    int *StopPos1;
    int *StopPos2;
    int *StopPos3;

    std::map<int, int> degeneracy;
};

#endif
