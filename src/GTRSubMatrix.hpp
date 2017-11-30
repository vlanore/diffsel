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

#ifndef GTRSUBMATRIX_H
#define GTRSUBMATRIX_H

#include "BiologicalSequences.hpp"  //FIXME only used for Naa (const int)
#include "SubMatrix.hpp"

class GTRSubMatrix : public virtual SubMatrix {
  public:
    GTRSubMatrix(int inNstate, const double *rr, const double *stat, bool innormalise = false);
    ~GTRSubMatrix() = default;

    int GetNRelativeRate() const { return Nrr; }
    double RelativeRate(int i, int j) const { return mRelativeRate[rrindex(i, j, GetNstate())]; }

    static int rrindex(int i, int j, int nstate) {
        return (i < j) ? (2 * nstate - i - 1) * i / 2 + j - i - 1
                       : (2 * nstate - j - 1) * j / 2 + i - j - 1;
    }

    // make a copy of the entries (not of the pointer)
    void CopyStationary(const double *instat) const;

  protected:
    void ComputeArray(int i) const override;
    void ComputeStationary() const override {}

    // data members
    const double *mRelativeRate;
    int Nrr;
};

#endif
