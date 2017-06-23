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

    void RegisterWith(const TaxonSet *intaxset) {
        std::cerr << "register\n";
        if (taxset->GetNtaxa() != intaxset->GetNtaxa()) {
            std::cerr << "error in register seq\n";
            exit(1);
        }
        auto tmp = new int *[GetNtaxa()];
        for (int i = 0; i < GetNtaxa(); i++) {
            int k = 0;
            while ((k < GetNtaxa()) && (taxset->GetTaxon(k) != intaxset->GetTaxon(i))) {
                k++;
            }
            if (k == GetNtaxa()) {
                std::cerr << "error in register seq : overflow\n";
                exit(1);
            }
            tmp[i] = Data[k];
        }
        for (int i = 0; i < GetNtaxa(); i++) {
            Data[i] = tmp[i];
        }
        delete taxset;
        taxset = intaxset;
        delete[] tmp;
        std::cerr << "register ok\n";
    }

    void Unclamp() {
        for (int i = 0; i < Ntaxa; i++) {
            for (int j = 0; j < Nsite; j++) {
                Data[i][j] = unknown;
            }
        }
    }

    int GetNstate() const { return statespace->GetNstate(); }

    // the list of taxa
    const TaxonSet *GetTaxonSet() const { return taxset; }

    int GetNsite() const { return Nsite; }

    int GetNtaxa() const { return taxset->GetNtaxa(); }

    bool isMissing(int taxon, int site) const { return Data[taxon][site] == -1; }

    bool AllMissingTaxon(int tax) const {
        bool ret = true;
        int site = 0;
        while ((site < GetNsite()) && ret) {
            ret &= static_cast<int>(Data[tax][site] == unknown);
            site++;
        }
        return ret;
    }

    bool AllMissingTaxon(std::string taxname) const {
        int index = taxset->GetTaxonIndex(taxname);
        if (index == -1) {
            std::cerr << "error in all missing taxon: did not recognize " << taxname << '\n';
            exit(1);
        }
        return AllMissingTaxon(index);
    }

    bool AllMissingColumn(int site) const {
        bool ret = true;
        int tax = 0;
        while ((tax < GetNtaxa()) && ret) {
            ret &= static_cast<int>(Data[tax][site] == unknown);
            tax++;
        }
        return ret;
    }

    bool NoMissingColumn(int site) const {
        bool ret = true;
        int tax = 0;
        while ((tax < GetNtaxa()) && ret) {
            ret &= static_cast<int>(Data[tax][site] != unknown);
            tax++;
        }
        return ret;
    }

    bool ConstantColumn(int site) const {
        bool ret = true;
        int tax = 0;
        while ((tax < GetNtaxa()) && (Data[tax][site] == unknown)) {
            tax++;
        }

        if (tax < GetNtaxa()) {
            int refstate = Data[tax][site];

            while ((tax < GetNtaxa()) && ret) {
                if (Data[tax][site] != -1) {
                    ret &= static_cast<int>(Data[tax][site] == refstate);
                }
                tax++;
            }
        }
        return ret;
    }

    void SetState(int taxon, int site, int state) { Data[taxon][site] = state; }

    int GetState(int taxon, int site) const { return Data[taxon][site]; }

    void GetEmpiricalFreq(double *in) const;

    void GetSiteEmpiricalFreq(double **in, double pseudocount = 0) const;

    void ToStream(std::ostream &os) const;
    void ToStreamTriplet(std::ostream &os) const;
    int GetNonMissingTriplet() const;
    void ToFasta(std::ostream &os) const;

    void DeleteConstantSites() {
        int i = 0;
        int j = 0;
        int Eliminated = 0;
        while (i < Nsite) {
            int k = 0;
            while ((k < Ntaxa) && (Data[k][i] == unknown)) {
                k++;
            }
            if (k < Ntaxa) {
                int a = Data[k][i];
                k++;
                while ((k < Ntaxa) && ((Data[k][i] == unknown) || (Data[k][i] == a))) {
                    k++;
                }
                if (k == Ntaxa) {
                    Eliminated++;
                } else {
                    for (int k = 0; k < Ntaxa; k++) {
                        Data[k][j] = Data[k][i];
                    }
                    j++;
                }
            }
            i++;
        }
        Nsite -= Eliminated;
        std::cout << "number of positions eliminated : " << Eliminated << '\n';
    }

    void GetMissingCellsFromTemplate(SequenceAlignment *from) {
        int fromnsite = from->GetNsite();

        for (int i = 0; i < Nsite; i++) {
            for (int j = 0; j < Ntaxa; j++) {
                int ii = i % fromnsite;
                int state = from->GetState(j, ii);
                if (state == unknown) {
                    Data[j][i] = unknown;
                }
            }
        }
    }

    double GetMeanDiversity() const {
        int Nstate = GetNstate();
        std::vector<int> found(Nstate);

        double mean = 0;
        for (int i = 0; i < Nsite; i++) {
            for (int k = 0; k < Nstate; k++) {
                found[k] = 0;
            }
            for (int j = 0; j < Ntaxa; j++) {
                int state = GetState(j, i);
                if (state != unknown) {
                    found[state] = 1;
                }
            }
            double div = 0;
            for (int k = 0; k < Nstate; k++) {
                div += found[k];
            }
            mean += div;
        }
        mean /= Nsite;
        return mean;
    }

    void GetTaxEmpiricalFreq(double **taxfreq) const {
        int Nstate = GetNstate();
        for (int j = 0; j < Ntaxa; j++) {
            for (int k = 0; k < Nstate; k++) {
                taxfreq[j][k] = 0;
            }
        }

        for (int i = 0; i < Nsite; i++) {
            for (int j = 0; j < Ntaxa; j++) {
                int state = GetState(j, i);
                if (state != unknown) {
                    taxfreq[j][state]++;
                }
            }
        }

        for (int j = 0; j < Ntaxa; j++) {
            double total = 0;
            for (int k = 0; k < Nstate; k++) {
                total += taxfreq[j][k];
            }
            for (int k = 0; k < Nstate; k++) {
                taxfreq[j][k] /= total;
            }
        }
    }

    double CompositionalHeterogeneity(std::ostream *os) const {
        int Nstate = GetNstate();
        auto taxfreq = new double *[Ntaxa];
        for (int j = 0; j < Ntaxa; j++) {
            taxfreq[j] = new double[Nstate];
            for (int k = 0; k < Nstate; k++) {
                taxfreq[j][k] = 0;
            }
        }

        for (int i = 0; i < Nsite; i++) {
            for (int j = 0; j < Ntaxa; j++) {
                int state = GetState(j, i);
                if (state != unknown) {
                    taxfreq[j][state]++;
                }
            }
        }

        // make global freqs out of tax-specific freqs
        auto globalfreq = new double[Nstate];
        for (int k = 0; k < Nstate; k++) {
            globalfreq[k] = 0;
            for (int j = 0; j < Ntaxa; j++) {
                globalfreq[k] += taxfreq[j][k];
            }
        }

        // normalise
        double total = 0;
        for (int k = 0; k < Nstate; k++) {
            total += globalfreq[k];
        }
        for (int k = 0; k < Nstate; k++) {
            globalfreq[k] /= total;
        }
        for (int j = 0; j < Ntaxa; j++) {
            double total = 0;
            for (int k = 0; k < Nstate; k++) {
                total += taxfreq[j][k];
            }
            for (int k = 0; k < Nstate; k++) {
                taxfreq[j][k] /= total;
                if (os != nullptr) {
                    (*os) << taxfreq[j][k] << '\t';
                }
            }
            if (os != nullptr) {
                (*os) << '\n';
            }
        }
        if (os != nullptr) {
            (*os) << '\n';
        }

        // compute max distance
        double maxdist = 0;
        for (int j = 0; j < Ntaxa; j++) {
            double dist = 0;
            for (int k = 0; k < Nstate; k++) {
                double tmp = (taxfreq[j][k] - globalfreq[k]);
                dist += tmp * tmp;
            }
            if (maxdist < dist) {
                maxdist = dist;
            }
        }

        delete[] globalfreq;
        for (int j = 0; j < Ntaxa; j++) {
            delete[] taxfreq[j];
        }
        delete[] taxfreq;

        return maxdist;
    }

    // data fields

    int Ntaxa{0};
    int Nsite{0};
    const TaxonSet *taxset{nullptr};
    const StateSpace *statespace{nullptr};
    int **Data{nullptr};
};

class FileSequenceAlignment : public SequenceAlignment {
  public:
    explicit FileSequenceAlignment(std::istream &is);
    FileSequenceAlignment(std::string filename, int fullline = 0);

  private:
    int ReadDataFromFile(std::string filespec, int forceinterleaved = 0);
    int ReadNexus(std::string filespec);
    int ReadSpecial(std::string filename);
    int TestPhylipSequential(std::string filespec);
    void ReadPhylipSequential(std::string filespec);
    int TestPhylip(std::string filespec, int repeattaxa);
    void ReadPhylip(std::string filespec, int repeattaxa);

    std::vector<std::string> SpeciesNames;
};

#endif  // SEQUENCEALIGNMENT_H
