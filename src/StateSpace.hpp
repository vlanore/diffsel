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

#ifndef STATESPACE_H
#define STATESPACE_H

#include <doctest.h>
#include <algorithm>
#include <string>
#include <vector>
#include "BiologicalSequences.hpp"

class StateSpace {
  public:
    virtual ~StateSpace() = default;

    virtual std::string GetState(int state) const = 0;
    virtual int GetNstate() const = 0;
    virtual int GetState(std::string from) const = 0;

    virtual bool isCompatible(int state1, int state2) const {
        return ((state1 == unknown) || (state2 == unknown) || (state1 == state2));
    }
};

// simple state space: assumes that states are referred to using a one-letter code
class SimpleStateSpace : public StateSpace {
  public:
    SimpleStateSpace(const std::vector<char> &alphabet, const std::vector<char> &alphabetSet)
        : Alphabet(alphabet), AlphabetSet(alphabetSet) {}

    int GetState(std::string from) const override {
        if (from.length() != 1) {
            std::cerr << "error in SimpleStateSpace\n";
            exit(1);
        }
        char c = from[0];
        auto pos = std::find(AlphabetSet.begin(), AlphabetSet.end(), c);
        if (pos == AlphabetSet.end()) {
            std::cerr << "error: does not recognise character " << c << '\n';
            exit(1);
        }
        int index = pos - AlphabetSet.begin();
        if (index >= 2 * GetNstate()) {  // if input is outside of Alphabet
            return unknown;
        } else {
            return index % GetNstate();  // return correct index even if lowercase input
        }
    }

    int GetNstate() const override { return Alphabet.size(); }

    std::string GetState(int state) const override {
        if (state == unknown) {
            return "-";
        } else {
            return std::string(1, Alphabet.at(state));
        }
    }

  private:
    const std::vector<char> &Alphabet;
    const std::vector<char> &AlphabetSet;
};

class DNAStateSpace : public SimpleStateSpace {
  public:
    DNAStateSpace(const std::vector<char> &alphabet = DNAletters,
                  const std::vector<char> &alphabetSet = DNAset)
        : SimpleStateSpace(alphabet, alphabetSet) {}
};

class RNAStateSpace : public SimpleStateSpace {
  public:
    RNAStateSpace(const std::vector<char> &alphabet = RNAletters,
                  const std::vector<char> &alphabetSet = RNAset)
        : SimpleStateSpace(alphabet, alphabetSet) {}
};

class ProteinStateSpace : public SimpleStateSpace {
  public:
    ProteinStateSpace(const std::vector<char> &alphabet = AminoAcids,
                      const std::vector<char> &alphabetSet = AAset)
        : SimpleStateSpace(alphabet, alphabetSet) {}
};

/*
#### TESTS #########################################################################################
*/
TEST_CASE("StateSpace tests") {
    DNAStateSpace myDnaSS;
    CHECK(myDnaSS.GetState(0) == "A");
    CHECK(myDnaSS.GetState(1) == "C");
    CHECK(myDnaSS.GetState(2) == "G");
    CHECK(myDnaSS.GetState(3) == "T");
    CHECK(myDnaSS.GetNstate() == 4);
    CHECK(myDnaSS.GetState("A") == 0);
    CHECK(myDnaSS.GetState("C") == 1);
    CHECK(myDnaSS.GetState("G") == 2);
    CHECK(myDnaSS.GetState("T") == 3);
    CHECK(myDnaSS.GetState("a") == 0);
    CHECK(myDnaSS.GetState("c") == 1);
    CHECK(myDnaSS.GetState("g") == 2);
    CHECK(myDnaSS.GetState("t") == 3);
    CHECK(myDnaSS.GetState("*") == unknown);
}

#endif  // STATESPACE_H
