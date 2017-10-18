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

#include "SubMatrix.hpp"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
// #include "linalg.hpp"
using namespace std;

int SubMatrix::nuni = 0;
int SubMatrix::nunimax = 0;
int SubMatrix::nunisubcount = 0;
int SubMatrix::diagcount = 0;

double SubMatrix::nz = 0;
double SubMatrix::meanz = 0;
double SubMatrix::maxz = 0;

void SubMatrix::Create() {
    // Q = new double *[Nstate];
    // for (int i = 0; i < Nstate; i++) {
    //     Q[i] = new double[Nstate];
    // }
    Q = EMatrix(Nstate, Nstate);

    u = EMatrix(Nstate, Nstate);
    // u = new double *[Nstate];
    // for (int i = 0; i < Nstate; i++) {
    //     u[i] = new double[Nstate];
    // }

    invu = EMatrix(Nstate, Nstate);
    // invu = new double *[Nstate];
    // for (int i = 0; i < Nstate; i++) {
    //     invu[i] = new double[Nstate];
    // }

    v = EVector(Nstate);
    vi = EVector(Nstate);
    // v = new double[Nstate];
    // vi = new double[Nstate];

    mStationary = EVector(Nstate);
    oldStationary = new double[Nstate];

    UniMu = 1;
    mPow = new double**[UniSubNmax];
    for (int n = 0; n < UniSubNmax; n++) {
        mPow[n] = nullptr;
    }

    flagarray = new bool[Nstate];
    diagflag = false;
    statflag = false;
    for (int i = 0; i < Nstate; i++) {
        flagarray[i] = false;
    }
    powflag = false;

    // aux = new double *[Nstate];
    // for (int i = 0; i < Nstate; i++) {
    //     aux[i] = new double[Nstate];
    // }
}

// ---------------------------------------------------------------------------
//     ~SubMatrix()
// ---------------------------------------------------------------------------

SubMatrix::~SubMatrix() {
    for (int i = 0; i < Nstate; i++) {
        // delete[] Q[i];
        // delete[] u[i];
        // delete[] invu[i];
        // delete[] aux[i];
    }
    // delete[] Q;
    // delete[] u;
    // delete[] invu;

    if (mPow != nullptr) {
        for (int n = 0; n < UniSubNmax; n++) {
            if (mPow[n] != nullptr) {
                for (int i = 0; i < Nstate; i++) {
                    delete[] mPow[n][i];
                }
                delete[] mPow[n];
            }
        }
        delete[] mPow;
    }
    // delete[] mStationary;
    delete[] oldStationary;
    delete[] flagarray;
    // delete[] v;
    // delete[] vi;

    // delete[] aux;
}

// ---------------------------------------------------------------------------
//     void ScalarMul()
// ---------------------------------------------------------------------------
void SubMatrix::ScalarMul(double e) {
    for (int i = 0; i < Nstate; i++) {
        v[i] *= e;
        vi[i] *= e;
    }
    UniMu *= e;
}

// ---------------------------------------------------------------------------
//     Diagonalise()
// ---------------------------------------------------------------------------
int SubMatrix::Diagonalise() const {
    if (!ArrayUpdated()) {
        UpdateMatrix();
    }

    diagcount++;
    auto& stat = GetStationary();

    EMatrix a(Nstate, Nstate);

    for (int i = 0; i < Nstate; i++) {
        for (int j = 0; j < Nstate; j++) {
            a(i, j) = Q(i, j) * sqrt(stat[i] / stat[j]);
        }
    }

    solver.compute(a);
    v = solver.eigenvalues().real();
    vi = solver.eigenvalues().imag();
    u = solver.eigenvectors().real();

    for (int i = 0; i < Nstate; i++) {
        for (int j = 0; j < Nstate; j++) {
            invu(i, j) = u(j, i) * sqrt(stat[j]);
        }
    }
    for (int i = 0; i < Nstate; i++) {
        for (int j = 0; j < Nstate; j++) {
            u(i, j) /= sqrt(stat[i]);
        }
    }

    diagflag = true;
    return 0;
}

// ---------------------------------------------------------------------------
//     ComputeRate()
// ---------------------------------------------------------------------------
double SubMatrix::GetRate() const {
    if (!ArrayUpdated()) {
        UpdateStationary();
        for (int k = 0; k < Nstate; k++) {
            ComputeArray(k);
        }
    }
    for (int k = 0; k < Nstate; k++) {
        flagarray[k] = true;
    }
    double norm = 0;
    for (int i = 0; i < Nstate - 1; i++) {
        for (int j = i + 1; j < Nstate; j++) {
            norm += mStationary[i] * Q(i, j);
        }
    }
    return 2 * norm;
}

const EVector& SubMatrix::GetEigenVal() const {
    if (!diagflag) {
        Diagonalise();
    }
    return v;
}

const EMatrix& SubMatrix::GetEigenVect() const {
    if (!diagflag) {
        Diagonalise();
    }
    return u;
}

const EMatrix& SubMatrix::GetInvEigenVect() const {
    if (!diagflag) {
        Diagonalise();
    }
    return invu;
}

// ---------------------------------------------------------------------------
//     Update()
// ---------------------------------------------------------------------------

void SubMatrix::UpdateMatrix() const {
    UpdateStationary();
    for (int k = 0; k < Nstate; k++) {
        ComputeArray(k);
    }
    for (int k = 0; k < Nstate; k++) {
        flagarray[k] = true;
    }
    if (isNormalised()) {
        Normalise();
    }
    // CheckReversibility();
}

// ---------------------------------------------------------------------------
//     Normalise()
// ---------------------------------------------------------------------------

void SubMatrix::Normalise() const {
    double norm = GetRate();
    for (int i = 0; i < Nstate; i++) {
        for (int j = 0; j < Nstate; j++) {
            Q(i, j) /= norm;
        }
    }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//     Powers
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
void SubMatrix::ActivatePowers() const {
    if (!powflag) {
        if (!ArrayUpdated()) {
            UpdateMatrix();
        }

        UniMu = 0;
        for (int i = 0; i < Nstate; i++) {
            if (UniMu < fabs(Q(i, i))) {
                UniMu = fabs(Q(i, i));
            }
        }

        CreatePowers(0);
        for (int i = 0; i < Nstate; i++) {
            for (int j = 0; j < Nstate; j++) {
                mPow[0][i][j] = 0;
            }
        }
        for (int i = 0; i < Nstate; i++) {
            mPow[0][i][i] = 1;
        }
        for (int i = 0; i < Nstate; i++) {
            for (int j = 0; j < Nstate; j++) {
                mPow[0][i][j] += Q(i, j) / UniMu;
                if (mPow[0][i][j] < 0) {
                    cerr << "error in SubMatrix::ComputePowers: negative prob : ";
                    cerr << i << '\t' << j << '\t' << mPow[0][i][j] << '\n';
                    cerr << "Nstate : " << Nstate << '\n';
                    exit(1);
                }
            }
        }
        npow = 1;
        powflag = true;
    }
}

void SubMatrix::InactivatePowers() const {
    if (powflag) {
        for (int n = 0; n < UniSubNmax; n++) {
            if (mPow[n] != nullptr) {
                for (int i = 0; i < Nstate; i++) {
                    delete[] mPow[n][i];
                }
                delete[] mPow[n];
                mPow[n] = nullptr;
            }
        }
        nunimax += npow;
        nuni++;

        npow = 0;
        powflag = false;
    }
}

void SubMatrix::CreatePowers(int n) const {
    if (mPow[n] == nullptr) {
        mPow[n] = new double*[Nstate];
        for (int i = 0; i < Nstate; i++) {
            mPow[n][i] = new double[Nstate];
        }
    }
}

double SubMatrix::GetUniformizationMu() const {
    if (!powflag) {
        ActivatePowers();
    }
    return UniMu;
}

double SubMatrix::Power(int n, int i, int j) const {
    if (!powflag) {
        ActivatePowers();
    }
    if (n == 0) {
        return static_cast<double>(i == j);
    }
    if (n > UniSubNmax) {
        return Stationary(j);
    }
    if (n > npow) {
        ComputePowers(n);
    }
    return mPow[n - 1][i][j];
}

void SubMatrix::ComputePowers(int N) const {
    if (!powflag) {
        ActivatePowers();
    }
    if (N > npow) {
        for (int n = npow; n < N; n++) {
            CreatePowers(n);
            for (int i = 0; i < Nstate; i++) {
                for (int j = 0; j < Nstate; j++) {
                    double& t = mPow[n][i][j];
                    t = 0;
                    for (int k = 0; k < Nstate; k++) {
                        t += mPow[n - 1][i][k] * mPow[0][k][j];
                    }
                }
            }
        }
        npow = N;
    }
}

void SubMatrix::ToStream(ostream& os) const {
    os << GetNstate() << '\n';
    os << "stationaries:\n";
    for (int i = 0; i < GetNstate(); i++) {
        os << Stationary(i) << '\t';
    }
    os << '\n';

    os << "rate matrix:\n";
    for (int i = 0; i < GetNstate(); i++) {
        for (int j = 0; j < GetNstate(); j++) {
            os << Q(i, j) << '\t';
        }
        os << '\n';
    }

    os << '\n';
    for (int i = 0; i < GetNstate(); i++) {
        os << v[i] << '\t';
    }
    os << '\n';
}

double SubMatrix::SuffStatLogProb(PathSuffStat* suffstat) {
    double total = 0;
    auto& stat = GetStationary();
    for (auto i : suffstat->rootcount) {
        total += i.second * log(stat[i.first]);
    }
    for (auto i : suffstat->waitingtime) {
        total += i.second * (*this)(i.first, i.first);
    }
    for (auto i : suffstat->paircount) {
        total += i.second * log((*this)(i.first.first, i.first.second));
    }
    return total;
}
