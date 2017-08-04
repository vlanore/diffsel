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

#ifndef SEQUENCEALIGNMENT_H
#define SEQUENCEALIGNMENT_H

#include <vector>
#include "StateSpace.hpp"
#include "TaxonSet.hpp"

// this class works like an interface
// it does not do any job
class SequenceAlignment {
  public:
    SequenceAlignment() : Data(nullptr) {}
    SequenceAlignment(const SequenceAlignment &) = delete;
    virtual ~SequenceAlignment() = default;

    // the set of characters (A,C,G,T for nucleotides, etc..)
    const StateSpace *GetStateSpace() const { return statespace; }
    const TaxonSet *GetTaxonSet() const { return taxset; }  // the list of taxa
    int GetNstate() const { return statespace->GetNstate(); }
    int GetNsite() const { return Nsite; }
    int GetNtaxa() const { return taxset->GetNtaxa(); }
    void SetState(int taxon, int site, int state) { Data[taxon][site] = state; }
    int GetState(int taxon, int site) const { return Data[taxon][site]; }

    bool isMissing(int taxon, int site) const { return Data[taxon][site] == -1; }

    bool NoMissingColumn(int site) const {
        bool ret = true;
        int tax = 0;
        while ((tax < GetNtaxa()) && ret) {
            ret &= static_cast<int>(Data[tax][site] != unknown);
            tax++;
        }
        return ret;
    }

    void ToStream(std::ostream &os) const;

    // data fields
    int Ntaxa{0};
    int Nsite{0};

  protected:
    const TaxonSet *taxset{nullptr};
    const StateSpace *statespace{nullptr};
    int **Data{nullptr};
};

class FileSequenceAlignment : public SequenceAlignment {
  public:
    explicit FileSequenceAlignment(std::istream &is);
    explicit FileSequenceAlignment(std::string filename);

  private:
    void FailFormat();
    int ReadDataFromFile(std::string filespec, int forceinterleaved = 0);
    int TestPhylipSequential(std::string filespec);
    void ReadPhylipSequential(std::string filespec);
    int TestPhylip(std::string filespec, int repeattaxa);
    void ReadPhylip(std::string filespec, int repeattaxa);

    std::vector<std::string> SpeciesNames;
};

#endif  // SEQUENCEALIGNMENT_H
