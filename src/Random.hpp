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

#ifndef RANDOM_H
#define RANDOM_H

#include <vector>

#define MT_LEN 624  // (VL) required for magic

static const double Pi = 3.1415926535897932384626;

class Random {
  public:
    static const double INFPROB;

    Random(int seed = -1);

    static void InitRandom(int seed = -1);

    static int GetSeed();

    static double Uniform();
    static int ApproxBinomial(int N, double p);
    static int Poisson(double mu);
    static double Gamma(double alpha, double beta);
    static double sNormal();
    static double sExpo();
    static double sGamma(double);
    static double sGammanew(double);

    static int Choose(int);
    static int FiniteDiscrete(int n, const double *probarray);
    static void DrawFromUrn(std::vector<int> &, int n, int N);
    static int DrawFromDiscreteDistribution(const double *prob, int nstate);

    static double logGamma(double alpha);

    static double logMultivariateGamma(double a, int p);

    static double ProfileProposeMove(double *profile, int dim, double tuning, int n);
    static double RealVectorProposeMove(double *x, int dim, double tuning, int n);

  private:
    static int Seed;
    static int mt_index;
    static unsigned long long mt_buffer[MT_LEN];
};

#endif  // RANDOM_H
