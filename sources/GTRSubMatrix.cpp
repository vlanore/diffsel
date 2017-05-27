#include "GTRSubMatrix.hpp"
#include <iostream>
using namespace std;

// ---------------------------------------------------------------------------
//     GTRSubMatrix
// ---------------------------------------------------------------------------

GTRSubMatrix::GTRSubMatrix(int inNstate, const std::vector<double>& rr, const std::vector<double>& stat, bool innormalise)
    : SubMatrix(inNstate, innormalise), mRelativeRate(rr)  {
    Nrr = Nstate * (Nstate - 1) / 2;
        CopyStationary(stat);
}

void GTRSubMatrix::CopyStationary(const std::vector<double>& instat) {
    for (int k = 0; k < Nstate; k++) {
        mStationary[k] = instat[k];
    }
}

// ---------------------------------------------------------------------------
//     ComputeArray
// ---------------------------------------------------------------------------

void GTRSubMatrix::ComputeArray(int i) const {
	double total = 0;
	for (int j = 0; j < Nstate; j++) {
	    if (i != j) {
		Q[i][j] = RelativeRate(i, j) * mStationary[j];
		total += Q[i][j];
	    }
	}

	// should always ensure that the diagonal entry of the matrix Q[i][i] is
	// such that
	// the sum over all entries of the row is equal to 0
	Q[i][i] = -total;
}
