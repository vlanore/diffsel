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

// SubMatrix:
// this class implements
// the instant rate matrix of a substitution process
// but only the mathematical aspects of it
//
// If you need to define a new substitution process
// - derive a new class deriving (directly or indirectly) from SubMatrix
// in this class, implements the ComputeArray and ComputeStationary functions
// (these are the functions that construct the instant rates and the stationary
// probabilities of the
// matrix)
//
// - derive a new class deriving both from RandomSubMatrix, and from your new
// class
// in this class, implement SetParameters
// (this function is responsible for updating the internal parameters that the
// SubMatrix uses in
// ComputeArray And ComputeStationary,
// based on the values stored by the parent nodes)

#ifndef SUBMATRIX_H
#define SUBMATRIX_H

#include <doctest.h>
#include <Eigen/Dense>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include "PathSuffStat.hpp"
#include "Random.hpp"

using Vector = double *;
using ConstVect = const double *;
using Matrix = double **;
using EMatrix = Eigen::MatrixXd;
using EVector = Eigen::VectorXd;

class SubMatrix {
  protected:
    // these 2 pure virtual functions are the most essential component of the
    // SubMatrix class see GTRSubMatrix.cpp and CodonSubMatrix.cpp for examples

    // ComputeArray(int state) is in charge of computing the row of the rate
    // matrix corresponding to all possible rates of substitution AWAY from state
    virtual void ComputeArray(int state) const = 0;

    // ComputeStationary() is in charge of computing the vector of stationary
    // probabilities (equilibirum frequencies) of the substitution process
    virtual void ComputeStationary() const = 0;

    static const int UniSubNmax = 500;
    static int nunisubcount;
    static int nuni;
    static int nunimax;
    static int diagcount;

    static double meanz;
    static double maxz;
    static double nz;

    mutable bool powflag;
    mutable bool diagflag;
    mutable bool statflag;
    mutable bool *flagarray;

    int Nstate;
    mutable int npow;
    mutable double UniMu;
    bool normalise;

    Matrix *mPow;
    mutable EMatrix Q;            // Q : the infinitesimal generator matrix
    mutable EVector mStationary;  // the stationary probabilities of the matrix
    double *oldStationary;
    // mutable Matrix aux;          // an auxiliary matrix

    mutable Eigen::EigenSolver<EMatrix> solver;

    mutable EMatrix u;     // u : the matrix of eigen vectors
    mutable EMatrix invu;  // invu : the inverse of u
    mutable EVector v;     // v : eigenvalues
    mutable EVector vi;    // vi : imaginary part

    mutable int ndiagfailed;

  public:
    SubMatrix(int inNstate, bool innormalise = false)
        : Nstate(inNstate), normalise(innormalise), ndiagfailed(0) {
        Create();
    }

    SubMatrix(const SubMatrix &) = delete;
    virtual ~SubMatrix();

    void Create();  // manual allocation of internal resources (requires Nstate and UniSubNmax)

    // static getters
    static int GetUniSubCount() { return nunisubcount; }
    static int GetDiagCount() { return diagcount; }
    static double GetMeanUni() { return ((double)nunimax) / nuni; }

    // getters
    int GetNstate() const { return Nstate; }
    const EVector &GetStationary() const;
    double *OldGetStationary() const {
        if (!statflag) {
            UpdateStationary();
        }
        return oldStationary;
    }
    double GetRate() const;
    const EVector &GetEigenVal() const;
    const EMatrix &GetEigenVect() const;
    const EMatrix &GetInvEigenVect() const;
    int GetDiagStat() const { return ndiagfailed; }
    double GetUniformizationMu() const;

    static void ResetDiagCount() { diagcount = 0; }

    double Stationary(int i) const;
    double operator()(int /*i*/, int /*j*/) const;
    void ScalarMul(double e);

    bool isNormalised() const { return normalise; }
    void Normalise() const;

    virtual void CorruptMatrix();
    void UpdateMatrix() const;

    void ActivatePowers() const;
    void InactivatePowers() const;
    double Power(int n, int i, int j) const;

    virtual void ToStream(std::ostream &os) const;

    void BackwardPropagate(const double *up, double *down, double length) const;
    void ForwardPropagate(const double *down, double *up, double length) const;

    // uniformization resampling methods
    // CPU level 1
    void GetFiniteTimeTransitionProb(int state, double *down, double efflength) const;
    double GetFiniteTimeTransitionProb(int stateup, int statedown, double efflength) const;
    int DrawUniformizedTransition(int state, int statedown, int n) const;
    int DrawUniformizedSubstitutionNumber(int stateup, int statedown, double efflength) const;

    // used by accept-reject resampling method
    // CPU level 1
    int DrawOneStep(int state) const;
    double DrawWaitingTime(int state) const;
    int DrawFromStationary() const;

    double SuffStatLogProb(PathSuffStat *suffstat);

  protected:
    EVector GetRow(int i) const;
    void UpdateRow(int state) const;
    void UpdateStationary() const;

    void ComputePowers(int N) const;
    void CreatePowers(int n) const;

    bool ArrayUpdated() const;

    int Diagonalise() const;
};


//-------------------------------------------------------------------------
//	* Inline definitions
//-------------------------------------------------------------------------
inline double SubMatrix::operator()(int i, int j) const {
    if (!flagarray[i]) {
        UpdateRow(i);
    }
    return Q(i, j);
}

inline EVector SubMatrix::GetRow(int i) const {
    if (!flagarray[i]) {
        UpdateRow(i);
    }
    return Q.row(i);
}

inline const EVector &SubMatrix::GetStationary() const {
    if (!statflag) {
        UpdateStationary();
    }
    return mStationary;
}

inline double SubMatrix::Stationary(int i) const {
    if (!statflag) {
        UpdateStationary();
    }
    return mStationary[i];
}

inline void SubMatrix::CorruptMatrix() {
    diagflag = false;
    statflag = false;
    for (int k = 0; k < Nstate; k++) {
        flagarray[k] = false;
    }
    InactivatePowers();
}

inline bool SubMatrix::ArrayUpdated() const {
    bool qflag = true;
    for (int k = 0; k < Nstate; k++) {
        qflag &= static_cast<int>(flagarray[k]);
    }
    return qflag;
}

inline void SubMatrix::UpdateStationary() const {
    ComputeStationary();
    for (int i = 0; i < Nstate; i++) {
        oldStationary[i] = mStationary[i];
    }
    statflag = true;
}

inline void SubMatrix::UpdateRow(int state) const {
    if (isNormalised()) {
        UpdateMatrix();
    } else {
        if (!statflag) {
            UpdateStationary();
        }
        ComputeArray(state);
        flagarray[state] = true;
    }
}

inline void SubMatrix::BackwardPropagate(const double *up, double *down, double length) const {
    auto &eigenvect = GetEigenVect();
    auto &inveigenvect = GetInvEigenVect();
    auto &eigenval = GetEigenVal();

    int matSize = GetNstate();

    auto aux = new double[GetNstate()];

    for (int i = 0; i < GetNstate(); i++) {
        aux[i] = 0;
    }
    for (int i = 0; i < GetNstate(); i++) {
        for (int j = 0; j < GetNstate(); j++) {
            aux[i] += inveigenvect(i, j) * up[j];
        }
    }

    for (int i = 0; i < GetNstate(); i++) {
        aux[i] *= exp(length * eigenval[i]);
    }

    for (int i = 0; i < GetNstate(); i++) {
        down[i] = 0;
    }

    for (int i = 0; i < GetNstate(); i++) {
        for (int j = 0; j < GetNstate(); j++) {
            down[i] += eigenvect(i, j) * aux[j];
        }
    }

    for (int i = 0; i < GetNstate(); i++) {
        if (std::isnan(down[i])) {
            std::cerr << "error in back prop\n";
            for (int j = 0; j < GetNstate(); j++) {
                std::cerr << up[j] << '\t' << down[j] << '\t' << Stationary(j) << '\n';
            }
            exit(1);
        }
    }
    double maxup = 0;
    for (int k = 0; k < matSize; k++) {
        if (up[k] < 0) {
            std::cerr << "error in backward propagate: negative prob : " << up[k] << "\n";
            // down[k] = 0;
        }
        if (maxup < up[k]) {
            maxup = up[k];
        }
    }
    double max = 0;
    for (int k = 0; k < matSize; k++) {
        if (down[k] < 0) {
            down[k] = 0;
        }
        if (max < down[k]) {
            max = down[k];
        }
    }
    if (maxup == 0) {
        std::cerr << "error in backward propagate: null up array\n";
        exit(1);
    }
    if (max == 0) {
        std::cerr << "error in backward propagate: null array\n";
        for (int k = 0; k < matSize; k++) {
            std::cerr << up[k] << '\t' << down[k] << '\n';
        }
        std::cerr << '\n';
        exit(1);
    }
    down[matSize] = up[matSize];

    delete[] aux;
}

inline void SubMatrix::ForwardPropagate(const double *down, double *up, double length) const {
    auto &eigenvect = GetEigenVect();
    auto &inveigenvect = GetInvEigenVect();
    auto &eigenval = GetEigenVal();

    auto aux = new double[GetNstate()];

    for (int i = 0; i < GetNstate(); i++) {
        aux[i] = 0;
    }

    for (int i = 0; i < GetNstate(); i++) {
        for (int j = 0; j < GetNstate(); j++) {
            aux[i] += down[j] * eigenvect(j, i);
        }
    }

    for (int i = 0; i < GetNstate(); i++) {
        aux[i] *= exp(length * eigenval[i]);
    }

    for (int i = 0; i < GetNstate(); i++) {
        up[i] = 0;
    }

    for (int i = 0; i < GetNstate(); i++) {
        for (int j = 0; j < GetNstate(); j++) {
            up[i] += aux[j] * inveigenvect(j, i);
        }
    }

    delete[] aux;
}

inline double SubMatrix::GetFiniteTimeTransitionProb(int stateup, int statedown,
                                                     double efflength) const {
    auto &invp = GetInvEigenVect();
    auto &p = GetEigenVect();
    auto &l = GetEigenVal();
    double tot = 0;
    for (int i = 0; i < GetNstate(); i++) {
        tot += p(stateup, i) * exp(efflength * l[i]) * invp(i, statedown);
    }
    return tot;
}

inline void SubMatrix::GetFiniteTimeTransitionProb(int state, double *p, double efflength) const {
    double *p1 = new double[GetNstate()];
    for (int k = 0; k < GetNstate(); k++) {
        p1[k] = 0;
    }
    p1[state] = 1;
    ForwardPropagate(p1, p, efflength);
    double tot = 0;
    for (int k = 0; k < GetNstate(); k++) {
        tot += p[k];
    }
    if (fabs(1 - tot) > 1e-5) {
        std::cerr << "error in forward propagate: normalization : " << tot << '\t' << fabs(1 - tot)
                  << '\n';
        std::cerr << "eff length : " << efflength << '\n';
        ToStream(std::cerr);
        exit(1);
    }
    delete[] p1;
}

inline int SubMatrix::DrawUniformizedTransition(int state, int statedown, int n) const {
    double *p = new double[GetNstate()];
    double tot = 0;
    for (int l = 0; l < GetNstate(); l++) {
        tot += Power(1, state, l) * Power(n, l, statedown);
        p[l] = tot;
    }

    double s = tot * Random::Uniform();
    int k = 0;
    while ((k < GetNstate()) && (s > p[k])) {
        k++;
    }
    delete[] p;
    if (k == GetNstate()) {
        std::cerr << "error in DrawUniformizedTransition: overflow\n";
        throw;
    }
    return k;
}


inline int SubMatrix::DrawUniformizedSubstitutionNumber(int stateup, int statedown,
                                                        double efflength) const {
    double mu = GetUniformizationMu();
    double fact = exp(-efflength * mu);
    int m = 0;
    double total = (stateup == statedown) * fact;
    double Z = GetFiniteTimeTransitionProb(stateup, statedown, efflength);
    double q = Random::Uniform() * Z;

    while ((m < UniSubNmax) && (total < q)) {
        m++;
        fact *= mu * efflength / m;
        total += Power(m, stateup, statedown) * fact;
        if ((total - Z) > 1e-12) {
            std::cerr << "error in DrawUniformizedSubstitutionNumber: normalising constant\n";
            std::cerr << total << '\t' << Z << '\n';
            std::cerr << mu << '\n';
            std::cerr << m << '\n';
            std::cerr << stateup << '\t' << statedown << '\n';

            ToStream(std::cerr);
            throw;
        }
    }
    if (m == UniSubNmax) {
        nunisubcount++;
    }
    return m;
}

inline double SubMatrix::DrawWaitingTime(int state) const {
    EVector row = GetRow(state);
    double t = Random::sExpo() / (-row[state]);
    return t;
}

inline int SubMatrix::DrawOneStep(int state) const {
    EVector row = GetRow(state);
    double p = -row[state] * Random::Uniform();
    int k = -1;
    double tot = 0;
    do {
        k++;
        if (k != state) {
            tot += row[k];
        }
    } while ((k < GetNstate()) && (tot < p));
    if (k == GetNstate()) {
        std::cerr << "error in DrawOneStep\n";
        std::cerr << GetNstate() << '\n';
        for (int k = 0; k < GetNstate(); k++) {
            std::cerr << row[k] << '\n';
        }
        exit(1);
    }
    return k;
}

/*
#### TESTS #########################################################################################
*/
TEST_CASE("SubMatrix tests") {
    class MyMatrix : public SubMatrix {
      protected:
        void ComputeArray(int) const override {}

        void ComputeStationary() const override {
            for (int i = 0; i < 4; ++i) {
                mStationary[i] = 0.25;
            }
        }

      public:
        MyMatrix() : SubMatrix(4) {}

        double &ref(int i, int j) { return Q(i, j); }

        using SubMatrix::Diagonalise;
    };

    MyMatrix myMatrix{};
    CHECK(myMatrix.GetNstate() == 4);

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            if (i == j) {
                myMatrix.ref(i, j) = -0.75;
            } else {
                myMatrix.ref(i, j) = 0.25;
            }
        }
    }

    myMatrix.Diagonalise();

    std::stringstream ss;
    myMatrix.ToStream(ss);
    CHECK(ss.str() ==
          "4\n\
stationaries:\n\
0.25\t0.25\t0.25\t0.25\t\n\
rate matrix:\n\
-0.75\t0.25\t0.25\t0.25\t\n\
0.25\t-0.75\t0.25\t0.25\t\n\
0.25\t0.25\t-0.75\t0.25\t\n\
0.25\t0.25\t0.25\t-0.75\t\n\
\n\
-1\t0\t-1\t-1\t\n");
}

#endif  // SUBMATRIX_H
