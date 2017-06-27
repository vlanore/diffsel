#ifndef TAXONSET_H
#define TAXONSET_H

#include <iostream>
#include <map>
#include <string>
#include <vector>

class Tree;
class Link;

class TaxonSet {
  public:
    TaxonSet(const std::vector<std::string> &names, int ntaxa) {
        Ntaxa = ntaxa;
        taxlist.resize(Ntaxa);
        for (int i = 0; i < ntaxa; i++) {
            if (taxmap[names[i]] != 0) {
                std::cerr << "found several taxa with same name : " << names[i] << '\n';
                exit(1);
            }
            taxlist[i] = names[i];
            taxmap[names[i]] = i + 1;
        }
    }

    int GetNtaxa() const { return Ntaxa; }

    std::string GetTaxon(int index) const { return taxlist[index]; }

    int GetTaxonIndex(std::string intaxon) const { return taxmap[intaxon] - 1; }

    void ToStream(std::ostream &os) const {
        os << Ntaxa << '\n';
        for (auto e : taxlist) {
            os << e << '\n';
        }
    }

  private:
    int Ntaxa;
    mutable std::map<std::string, int> taxmap;
    std::vector<std::string> taxlist;
};

#endif  // TAXONSET_H
