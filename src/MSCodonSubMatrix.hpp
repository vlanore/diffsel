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

#ifndef MSCODONSUBMATRIX_H
#define MSCODONSUBMATRIX_H

#include <iostream>
#include "CodonSubMatrix.hpp"

struct FitnessProxy {
    virtual double GetFitness(int) const = 0;
};

class SparseFitness : public FitnessProxy {
    Eigen::Ref<Eigen::VectorXd> fitness0;
    Eigen::Ref<Eigen::VectorXd> fitness;
    Eigen::Ref<Eigen::Matrix<bool, Eigen::Dynamic, 1>> ind_conv;

  public:
    SparseFitness(Eigen::Ref<Eigen::VectorXd> fitness0, Eigen::Ref<Eigen::VectorXd> fitness,
                  Eigen::Ref<Eigen::Matrix<bool, Eigen::Dynamic, 1>> ind_conv)
        : fitness0(fitness0), fitness(fitness), ind_conv(ind_conv) {}

    double GetFitness(int aastate) const final {
        return (ind_conv[aastate] ? fitness[aastate] : fitness0[aastate]) + 1e-6;
    }
};

// square root
class MGSRFitnessSubMatrix : public MGCodonSubMatrix {
  public:
    MGSRFitnessSubMatrix(const CodonStateSpace* instatespace, const SubMatrix* inNucMatrix,
                         FitnessProxy& proxy, bool innormalise = false)
        : SubMatrix(instatespace->GetNstate(), innormalise),
          CodonSubMatrix(instatespace, innormalise),
          MGCodonSubMatrix(instatespace, inNucMatrix, innormalise),
          proxy(proxy) {}

  private:
    double GetFitness(int aastate) const { return proxy.GetFitness(aastate); }

    void ComputeArray(int i) const override;
    void ComputeStationary() const override;
    FitnessProxy& proxy;
};

// mutation selection
class MGMSFitnessSubMatrix : public MGCodonSubMatrix {
  public:
    MGMSFitnessSubMatrix(const CodonStateSpace* instatespace, const SubMatrix* inNucMatrix,
                         FitnessProxy& proxy, bool innormalise = false)
        : SubMatrix(instatespace->GetNstate(), innormalise),
          CodonSubMatrix(instatespace, innormalise),
          MGCodonSubMatrix(instatespace, inNucMatrix, innormalise),
          proxy(proxy) {}

  private:
    double GetFitness(int aastate) const { return proxy.GetFitness(aastate); }

    void ComputeArray(int i) const override;
    void ComputeStationary() const override;

    FitnessProxy& proxy;
};

#endif
