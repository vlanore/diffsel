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

#include "CodonSubMatrix.hpp"
using namespace std;

void MGCodonSubMatrix::ComputeArray(int i) const {
    double total = 0;
    for (int j = 0; j < GetNstate(); j++) {
        if (i != j) {
            int pos = GetDifferingPosition(i, j);
            if ((pos != -1) && (pos != 3)) {
                int a = GetCodonPosition(pos, i);
                int b = GetCodonPosition(pos, j);
                Q(i, j) = (*NucMatrix)(a, b);
                total += Q(i, j);
            } else {
                Q(i, j) = 0;
            }
        }
    }
    Q(i, i) = -total;
}

void MGCodonSubMatrix::ComputeStationary() const {
    // compute stationary probabilities
    double total = 0;
    for (int i = 0; i < GetNstate(); i++) {
        mStationary[i] = NucMatrix->Stationary(GetCodonPosition(0, i)) *
                         NucMatrix->Stationary(GetCodonPosition(1, i)) *
                         NucMatrix->Stationary(GetCodonPosition(2, i));
        total += mStationary[i];
    }

    // renormalize stationary probabilities
    for (int i = 0; i < GetNstate(); i++) {
        mStationary[i] /= total;
    }
}

/*
void MGCodonSubMatrix::ComputeNucArrays() {
    cerr << "error: in mg codon sub matrix compute nuc array\n";
    exit(1);
    for (int i = 0; i < Nnuc; i++) {
        for (int j = 0; j < Nnuc; j++) {
            if (i != j) {
                synnucarray[i][j] = (*NucMatrix)(i, j);
                nonsynnucarray[i][j] = (*NucMatrix)(i, j);
            }
        }
    }
    const double *stat = NucMatrix->GetStationary();
    const int *stoppos1 = GetCodonStateSpace()->GetStopPos1();
    const int *stoppos2 = GetCodonStateSpace()->GetStopPos2();
    const int *stoppos3 = GetCodonStateSpace()->GetStopPos3();
    stopstat = 0;
    for (int i = 0; i < GetCodonStateSpace()->GetNstop(); i++) {
        stopstat += stat[stoppos1[i]] * stat[stoppos2[i]] * stat[stoppos3[i]];
    }
}
*/

void MGOmegaCodonSubMatrix::ComputeArray(int i) const {
    double total = 0;
    for (int j = 0; j < GetNstate(); j++) {
        if (i != j) {
            int pos = GetDifferingPosition(i, j);
            if ((pos != -1) && (pos != 3)) {
                int a = GetCodonPosition(pos, i);
                int b = GetCodonPosition(pos, j);
                if (a == b) {
                    cerr << GetCodonStateSpace()->GetState(i) << '\t'
                         << GetCodonStateSpace()->GetState(j) << '\n';
                    cerr << pos << '\n';
                    exit(1);
                }
                Q(i, j) = (*NucMatrix)(a, b);
                if (!Synonymous(i, j)) {
                    Q(i, j) *= GetOmega();
                }
            } else {
                Q(i, j) = 0;
            }
            total += Q(i, j);
        }
    }
    Q(i, i) = -total;
    if (total < 0) {
        cerr << "negative rate away\n";
        exit(1);
    }
}

/*
void MGOmegaCodonSubMatrix::ComputeNucArrays() {
    for (int i = 0; i < Nnuc; i++) {
        for (int j = 0; j < Nnuc; j++) {
            if (i != j) {
                synnucarray[i][j] = (*NucMatrix)(i, j);
                nonsynnucarray[i][j] = omega * (*NucMatrix)(i, j);
            }
        }
    }
    const double *stat = NucMatrix->GetStationary();
    const int *stoppos1 = GetCodonStateSpace()->GetStopPos1();
    const int *stoppos2 = GetCodonStateSpace()->GetStopPos2();
    const int *stoppos3 = GetCodonStateSpace()->GetStopPos3();
    stopstat = 0;
    for (int i = 0; i < GetCodonStateSpace()->GetNstop(); i++) {
        stopstat += stat[stoppos1[i]] * stat[stoppos2[i]] * stat[stoppos3[i]];
    }
}
*/

/*
void AminoAcidReducedCodonSubMatrix::ComputeStationary() {
    // Stat[a] = \sum_{i|a} CodonStat[i]

    for (int a = 0; a < GetNstate(); a++) {
        mStationary[a] = 0;
    }

    for (int i = 0; i < GetCodonStateSpace()->GetNstate(); i++) {
        mStationary[GetCodonStateSpace()->Translation(i)] += GetCodonSubMatrix()->Stationary(i);
    }
}

void AminoAcidReducedCodonSubMatrix::ComputeArray(int a) {
    // Q[a][[b] = [ \sum _{i|a, j|b} CodonQ(i, j) ] / [ \sum_{i|a} CodonStat[i] ]

    for (int b = 0; b < GetNstate(); b++) {
        Q[a][b] = 0;
    }

    for (int i = 0; i < GetCodonStateSpace()->GetNstate(); i++) {
        if (GetCodonStateSpace()->Translation(i) == a) {
            for (int j = 0; j < GetCodonStateSpace()->GetNstate(); j++) {
                int b = GetCodonStateSpace()->Translation(j);
                if (b != a) {
                    Q[a][b] += GetCodonSubMatrix()->Stationary(i) * (*GetCodonSubMatrix())(i, j);
                }
            }
        }
    }

    for (int b = 0; b < GetNstate(); b++) {
        if (b != a) {
            Q[a][b] /= mStationary[a];
        }
    }

    double total = 0;
    for (int b = 0; b < GetNstate(); b++) {
        if (b != a) {
            total += Q[a][b];
        }
    }
    Q[a][a] = -total;
}
*/
