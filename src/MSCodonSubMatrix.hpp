#ifndef MSCODONSUBMATRIX_H
#define MSCODONSUBMATRIX_H

#include <iostream>
#include "CodonSubMatrix.hpp"

// square root
class MGSRFitnessSubMatrix : public MGCodonSubMatrix {
  public:
    MGSRFitnessSubMatrix(const CodonStateSpace *instatespace, const SubMatrix *inNucMatrix,
                         const double *infitness, bool innormalise = false)
        : SubMatrix(instatespace->GetNstate(), innormalise),
          CodonSubMatrix(instatespace, innormalise),
          MGCodonSubMatrix(instatespace, inNucMatrix, innormalise),
          fitness(infitness) {}

    double GetFitness(int aastate) const { return fitness[aastate] + 1e-6; }

  protected:
    // look at how ComputeArray and ComputeStationary are implemented in
    // CodonSubMatrix.cpp

    void ComputeArray(int i) const override;
    void ComputeStationary() const override;
    // data members
    const double *fitness;
};

// mutation selection
class MGMSFitnessSubMatrix : public MGCodonSubMatrix {
  public:
    MGMSFitnessSubMatrix(const CodonStateSpace *instatespace, const SubMatrix *inNucMatrix,
                         const double *infitness, bool innormalise = false)
        : SubMatrix(instatespace->GetNstate(), innormalise),
          CodonSubMatrix(instatespace, innormalise),
          MGCodonSubMatrix(instatespace, inNucMatrix, innormalise),
          fitness(infitness) {}

    double GetFitness(int aastate) const { return fitness[aastate] + 1e-6; }

  protected:
    // look at how ComputeArray and ComputeStationary are implemented in
    // CodonSubMatrix.cpp

    void ComputeArray(int i) const override;
    void ComputeStationary() const override;

    // data members
    const double *fitness;
};

#endif
