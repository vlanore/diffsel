#include "SubMatrix.hpp"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "linalg.hpp"
// #include "linalg2.hpp"
using namespace std;

int SubMatrix::nuni = 0;
int SubMatrix::nunimax = 0;
int SubMatrix::nunisubcount = 0;

double SubMatrix::nz = 0;
double SubMatrix::meanz = 0;
double SubMatrix::maxz = 0;

// ---------------------------------------------------------------------------
//     SubMatrix()
// ---------------------------------------------------------------------------
SubMatrix::SubMatrix(int inNstate, bool innormalise) : Nstate(inNstate), normalise(innormalise) {
    ndiagfailed = 0;
    // discn = 10;
    Create();
}

void SubMatrix::Create() {
    Q = new double *[Nstate];
    for (int i = 0; i < Nstate; i++) {
        Q[i] = new double[Nstate];
    }

    u = new double *[Nstate];
    for (int i = 0; i < Nstate; i++) {
        u[i] = new double[Nstate];
    }

    invu = new double *[Nstate];
    for (int i = 0; i < Nstate; i++) {
        invu[i] = new double[Nstate];
    }

    v = new double[Nstate];
    vi = new double[Nstate];

    mStationary = new double[Nstate];

    UniMu = 1;
    mPow = new double **[UniSubNmax];
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

    aux = new double *[Nstate];
    for (int i = 0; i < Nstate; i++) {
        aux[i] = new double[Nstate];
    }
}

// ---------------------------------------------------------------------------
//     ~SubMatrix()
// ---------------------------------------------------------------------------

SubMatrix::~SubMatrix() {
    for (int i = 0; i < Nstate; i++) {
        delete[] Q[i];
        delete[] u[i];
        delete[] invu[i];
        delete[] aux[i];
    }
    delete[] Q;
    delete[] u;
    delete[] invu;

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
    delete[] mStationary;
    delete[] flagarray;
    delete[] v;
    delete[] vi;

    delete[] aux;
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

    int nmax = 1000;
    double epsilon = 1e-20;
    int n = LinAlg::DiagonalizeRateMatrix(Q, mStationary, Nstate, v, u, invu, nmax, epsilon);
    bool failed = (n == nmax);
    if (failed) {
        cerr << "in submatrix: diag failed\n";
        cerr << "n : " << n << '\t' << nmax << '\n';
        ofstream os("mat");
        ToStream(os);
        exit(1);
    }
    diagflag = true;
    return static_cast<int>(failed);
}

/*
int SubMatrix::Diagonalise() const {
    if (!ArrayUpdated()) {
        UpdateMatrix();
    }

        // copy Q into aux
        for (int i=0; i<Nstate; i++)	{
                for (int j=0; j<Nstate; j++)	{
                        aux[i][j] = Q[i][j];
                }
        }

        double * w = new double[Nstate];
        int* iw = new int[Nstate];

        // diagonalise a into v and u
    cerr << "diag\n";
        int success = EigenRealGeneral(Nstate, aux, v, vi, u, iw, w);
    cerr << "diag ok\n";
    cerr << success << '\n';

        // copy u into aux
        for (int i=0; i<Nstate; i++)	{
                for (int j=0; j<Nstate; j++)	{
                        aux[i][j] = u[i][j];
                }
        }

        // invert a into invu
        InvertMatrix(aux, Nstate, w, iw, invu);

        // CheckDiag();

        delete[] w;
        delete[] iw;
        return success;
}
*/

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
            norm += mStationary[i] * Q[i][j];
        }
    }
    return 2 * norm;
}

double *SubMatrix::GetEigenVal() const {
    if (!diagflag) {
        Diagonalise();
    }
    return v;
}

double **SubMatrix::GetEigenVect() const {
    if (!diagflag) {
        Diagonalise();
    }
    return u;
}

double **SubMatrix::GetInvEigenVect() const {
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
            Q[i][j] /= norm;
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
            if (UniMu < fabs(Q[i][i])) {
                UniMu = fabs(Q[i][i]);
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
                mPow[0][i][j] += Q[i][j] / UniMu;
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
        mPow[n] = new double *[Nstate];
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
                    double &t = mPow[n][i][j];
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

void SubMatrix::ToStream(ostream &os) const {
    os << GetNstate() << '\n';
    os << "stationaries: \n";
    for (int i = 0; i < GetNstate(); i++) {
        os << Stationary(i) << '\t';
    }
    os << '\n';

    os << "rate matrix\n";
    for (int i = 0; i < GetNstate(); i++) {
        for (int j = 0; j < GetNstate(); j++) {
            os << Q[i][j] << '\t';
        }
        os << '\n';
    }

    os << '\n';
    for (int i = 0; i < GetNstate(); i++) {
        os << v[i] << '\t';
    }
    os << '\n';
}

void SubMatrix::CheckReversibility() const {
    double max = 0;
    int imax = 0;
    int jmax = 0;
    for (int i = 0; i < GetNstate(); i++) {
        for (int j = i + 1; j < GetNstate(); j++) {
            double tmp = fabs(Stationary(i) * Q[i][j] - Stationary(j) * Q[j][i]);
            if (max < tmp) {
                max = tmp;
                imax = i;
                jmax = j;
            }
        }
    }
    if (max > 1e-6) {
        cerr << "max irreversibility: " << max << '\n';
        cerr << imax << '\t' << jmax << '\t' << Stationary(imax) << '\t' << Q[imax][jmax] << '\t'
             << Stationary(jmax) << '\t' << Q[jmax][imax] << '\n';
        exit(1);
    }
}

void SubMatrix::ComputeExponential(double range, double **expo) const {
    Diagonalise();

    for (int i = 0; i < Nstate; i++) {
        double tot = 0;
        for (int j = 0; j < Nstate; j++) {
            tot += Q[i][j];
        }
        if (fabs(tot) > 1e-8) {
            cerr << "error : row does not sum to 0\n";
            exit(1);
        }

        tot = 0;
        for (int j = 0; j < Nstate; j++) {
            double tmp = 0;
            for (int k = 0; k < Nstate; k++) {
                tmp += u[i][k] * invu[k][j] * exp(v[k] * range);
            }
            expo[i][j] = tmp;
            tot += tmp;
        }
        if (fabs(tot - 1) > 1e-8) {
            cerr << "error : row does not sum to 1\n";
            cerr << "range : " << range << endl;
            // exit(1);
            throw 1;
        }
    }
}

void SubMatrix::ApproachExponential(double range, double **expo, int /*unused*/) const {
    for (int i = 0; i < Nstate; i++) {
        if (!flagarray[i]) {
            UpdateRow(i);
        }
    }

    for (int i = 0; i < GetNstate(); i++) {
        double tot = 0;
        for (int j = 0; j < GetNstate(); j++) {
            tot += Q[i][j];
        }
        if (fabs(tot) > 10e-6) {
            cerr << "Error in approachExp : row does not sum to 0 " << endl;
            cerr << "Sum : " << tot;
            cerr << Q[i][0] << " " << Q[i][1] << " " << Q[i][2] << " " << Q[i][3] << endl;
            exit(1);
        }
    }

    int z = 1;

    double maxdiag = 0;  // Max of diagonal coefficients
    for (int i = 0; i < GetNstate(); i++) {
        if (maxdiag < fabs(Q[i][i]) * range) {
            maxdiag = fabs(Q[i][i]) * range;
        }
    }

    if (maxdiag > 300) {
        const double *stat = GetStationary();
        // exp(range Q) = stat(Q)
        for (int i = 0; i < GetNstate(); i++) {
            for (int j = 0; j < GetNstate(); j++) {
                expo[i][j] = stat[j];
            }
        }
    }

    else {
        double precision = 0.01;

        while (z < 1024 && maxdiag > z * precision) {
            z *= 2;
        }

        nz++;
        meanz += z;
        if (maxz < z) {
            maxz = z;
        }

        for (int i = 0; i < GetNstate(); i++) {
            for (int j = 0; j < GetNstate(); j++) {
                expo[i][j] = range * Q[i][j] / z + (i == j ? 1 : 0);
            }
        }
        // y is now an approximation of exp(Q/z)

        this->PowerOf2(expo, z);
    }

    for (int i = 0; i < GetNstate(); i++) {
        double tot = 0;
        for (int j = 0; j < GetNstate(); j++) {
            tot += expo[i][j];
        }
        if (fabs(tot - 1) > 10e-6) {
            cerr << "Error in approachExp : row does not sum to 1 " << endl;
            cerr << "Sum : " << tot << endl;
            for (int a = 0; a < GetNstate(); a++) {
                for (int b = 0; b < GetNstate(); b++) {
                    cerr << Q[a][b] << '\t';
                }
                cerr << '\n';
            }
            cerr << '\n';
            cerr << maxdiag << '\t' << range << '\t' << z << '\n';
            cerr << '\n';
            for (int a = 0; a < GetNstate(); a++) {
                for (int b = 0; b < GetNstate(); b++) {
                    cerr << range * Q[a][b] / z << '\t';
                }
                cerr << '\n';
            }
            cerr << '\n';

            exit(1);
        }
    }
}

void SubMatrix::PowerOf2(double **y, int z) const {
    if (z == 1) {
        return;
    }
    for (int i = 0; i < GetNstate(); i++) {
        for (int j = 0; j < GetNstate(); j++) {
            aux[i][j] = 0;
            for (int k = 0; k < GetNstate(); k++) {
                aux[i][j] += y[i][k] * y[k][j];
            }
        }
    }
    for (int i = 0; i < GetNstate(); i++) {
        for (int j = 0; j < GetNstate(); j++) {
            y[i][j] = aux[i][j];
        }
    }
    PowerOf2(y, z / 2);
}

double SubMatrix::SuffStatLogProb(PathSuffStat *suffstat) {
    double total = 0;
    const double *stat = GetStationary();
    for (map<int, int>::iterator i = suffstat->rootcount.begin(); i != suffstat->rootcount.end();
         i++) {
        total += i->second * log(stat[i->first]);
    }
    for (map<int, double>::iterator i = suffstat->waitingtime.begin();
         i != suffstat->waitingtime.end(); i++) {
        total += i->second * (*this)(i->first, i->first);
    }
    for (map<pair<int, int>, int>::iterator i = suffstat->paircount.begin();
         i != suffstat->paircount.end(); i++) {
        total += i->second * log((*this)(i->first.first, i->first.second));
    }
    return total;
}