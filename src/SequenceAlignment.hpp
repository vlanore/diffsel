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
    int ReadDataFromFile(std::string filespec, int forceinterleaved = 0);
    int TestPhylipSequential(std::string filespec);
    void ReadPhylipSequential(std::string filespec);
    int TestPhylip(std::string filespec, int repeattaxa);
    void ReadPhylip(std::string filespec, int repeattaxa);

    std::vector<std::string> SpeciesNames;
};

#endif  // SEQUENCEALIGNMENT_H
