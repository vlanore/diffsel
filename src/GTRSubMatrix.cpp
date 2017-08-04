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

#include "GTRSubMatrix.hpp"
#include <iostream>
using namespace std;

// ---------------------------------------------------------------------------
//     GTRSubMatrix
// ---------------------------------------------------------------------------

GTRSubMatrix::GTRSubMatrix(int inNstate, const double *rr, const double *stat, bool innormalise)
    : SubMatrix(inNstate, innormalise) {
    Nrr = Nstate * (Nstate - 1) / 2;
    mRelativeRate = rr;
    if (stat != nullptr) {
        CopyStationary(stat);
    }
}

void GTRSubMatrix::CopyStationary(const double *instat) const {
    for (int k = 0; k < Nstate; k++) {
        mStationary[k] = instat[k];
    }
}

// ---------------------------------------------------------------------------
//     ComputeArray
// ---------------------------------------------------------------------------

void GTRSubMatrix::ComputeArray(int i) const {
    if (mRelativeRate != nullptr) {
        double total = 0;
        for (int j = 0; j < Nstate; j++) {
            if (i != j) {
                Q(i, j) = RelativeRate(i, j) * mStationary[j];
                total += Q(i, j);
            }
        }

        // should always ensure that the diagonal entry of the matrix Q(i, i) is
        // such that
        // the sum over all entries of the row is equal to 0
        Q(i, i) = -total;
    } else {
        for (int j = 0; j < Nstate; j++) {
            Q(i, j) = mStationary[j];
        }
        Q(i, i) -= 1;
    }
}
