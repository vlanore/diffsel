#include "MSCodonSubMatrix.hpp"
using namespace std;

void MGSRFitnessSubMatrix::ComputeStationary() const {
    // compute stationary probabilities
    double total = 0;
    for (int i = 0; i < GetNstate(); i++) {
        mStationary[i] = NucMatrix->Stationary(GetCodonPosition(0, i)) *
                         NucMatrix->Stationary(GetCodonPosition(1, i)) *
                         NucMatrix->Stationary(GetCodonPosition(2, i)) *
                         GetFitness(GetCodonStateSpace()->Translation(i));
        total += mStationary[i];
    }

    // renormalize stationary probabilities
    for (int i = 0; i < GetNstate(); i++) {
        mStationary[i] /= total;
    }
}

void MGSRFitnessSubMatrix::ComputeArray(int i) const {
    double total = 0;
    for (int j = 0; j < GetNstate(); j++) {
        if (i != j) {
            // grep the position (0, 1 or 2) at which the two codons differ
            // if the codons differ at more than one position then pos = -1 or 3 (and the rate
            // should be 0)
            int pos = GetDifferingPosition(i, j);
            if ((pos != -1) && (pos != 3)) {
                // get the initial and final nucleotides
                int a = GetCodonPosition(pos, i);
                int b = GetCodonPosition(pos, j);
                Q[i][j] = (*NucMatrix)(a, b);

                // if non synonymous, then multiply by the fixation bias
                if (!Synonymous(i, j)) {
                    Q[i][j] *= sqrt(((GetFitness(GetCodonStateSpace()->Translation(j))) /
                                     (GetFitness(GetCodonStateSpace()->Translation(i)))));
                }
            } else {
                // non-nearest neighbor coodns: rate is 0
                Q[i][j] = 0;
            }
            total += Q[i][j];
        }
    }

    Q[i][i] = -total;
    if (total < 0) {
        cerr << "negative rate away\n";
        exit(1);
    }
}

void MGMSFitnessSubMatrix::ComputeStationary() const {
    // compute stationary probabilities
    double total = 0;
    for (int i = 0; i < GetNstate(); i++) {
        mStationary[i] = NucMatrix->Stationary(GetCodonPosition(0, i)) *
                         NucMatrix->Stationary(GetCodonPosition(1, i)) *
                         NucMatrix->Stationary(GetCodonPosition(2, i)) *
                         GetFitness(GetCodonStateSpace()->Translation(i));
        total += mStationary[i];
    }

    // renormalize stationary probabilities
    // double min = 1;
    for (int i = 0; i < GetNstate(); i++) {
        mStationary[i] /= total;
    }
}

void MGMSFitnessSubMatrix::ComputeArray(int i) const {
    double total = 0;
    for (int j = 0; j < GetNstate(); j++) {
        if (i != j) {
            // grep the position (0, 1 or 2) at which the two codons differ
            // if the codons differ at more than one position then pos = -1 or 3 (and the rate
            // should be 0)
            int pos = GetDifferingPosition(i, j);
            if ((pos != -1) && (pos != 3)) {
                // get the initial and final nucleotides
                int a = GetCodonPosition(pos, i);
                int b = GetCodonPosition(pos, j);

                Q[i][j] = (*NucMatrix)(a, b);

                // calculate the scaled selection coefficient
                double S = 0;
                if (!Synonymous(i, j)) {
                    S = (log(GetFitness(GetCodonStateSpace()->Translation(j))) -
                         log(GetFitness(GetCodonStateSpace()->Translation(i))));
                }

                // first order development when S << 1
                if ((fabs(S)) < 1e-30) {
                    Q[i][j] *= 1 + S / 2;
                }
                // asymptotic value for very large S
                else if (S > 50) {
                    Q[i][j] *= S;
                }
                // if S << -1 : fixation probability is essentially 0
                else if (S < -50) {
                    Q[i][j] = 0;
                }
                // otherwise, standard formula for fixation probability as a function of S
                else {
                    Q[i][j] *= S / (1.0 - exp(-S));
                }

            } else {
                // non-nearest neighbor coodns: rate is 0
                Q[i][j] = 0;
            }
            total += Q[i][j];

            if (std::isinf(Q[i][j])) {
                cerr << "Q matrix infinite: " << Q[i][j] << '\n';
                exit(1);
            }

            if (Q[i][j] < 0) {
                cerr << "Q matrix negative: " << Q[i][j] << '\n';
                exit(1);
            }
        }
    }

    Q[i][i] = -total;

    if (total < 0) {
        cerr << "negative rate away\n";
        exit(1);
    }
}
