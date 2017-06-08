#ifndef MSCODONSUBMATRIX_H
#define MSCODONSUBMATRIX_H

#include <iostream>
#include "CodonSubMatrix.hpp"

// square root
class MGSRFitnessCodonUsageSubMatrix : public MGCodonSubMatrix {
  public:
    MGSRFitnessCodonUsageSubMatrix(const CodonStateSpace *instatespace,
                                   const SubMatrix *inNucMatrix, const double *infitness,
                                   const double *incodonusageselection, bool innormalise = false)
        : SubMatrix(instatespace->GetNstate(), innormalise),
          CodonSubMatrix(instatespace, innormalise),
          MGCodonSubMatrix(instatespace, inNucMatrix, innormalise),
          fitness(infitness),
          codonusageselection(incodonusageselection) {}

    double GetFitness(int aastate) const { return fitness[aastate] + 1e-6; }

    double GetCodonUsageSelection(int codonstate) const { return codonusageselection[codonstate]; }

  protected:
    // look at how ComputeArray and ComputeStationary are implemented in
    // CodonSubMatrix.cpp

    void ComputeArray(int i) const override;
    void ComputeStationary() const override;
    // data members
    const double *fitness;
    const double *codonusageselection;
};

// mutation selection
class MGMSFitnessCodonUsageSubMatrix : public MGCodonSubMatrix {
  public:
    MGMSFitnessCodonUsageSubMatrix(const CodonStateSpace *instatespace,
                                   const SubMatrix *inNucMatrix, const double *infitness,
                                   const double *incodonusageselection, bool innormalise = false)
        : SubMatrix(instatespace->GetNstate(), innormalise),
          CodonSubMatrix(instatespace, innormalise),
          MGCodonSubMatrix(instatespace, inNucMatrix, innormalise),
          fitness(infitness),
          codonusageselection(incodonusageselection) {}

    double GetFitness(int aastate) const { return fitness[aastate] + 1e-6; }

    double GetCodonUsageSelection(int codonstate) const { return codonusageselection[codonstate]; }

  protected:
    // look at how ComputeArray and ComputeStationary are implemented in
    // CodonSubMatrix.cpp

    void ComputeArray(int i) const override;
    void ComputeStationary() const override;

    // data members
    const double *fitness;
    const double *codonusageselection;
};

#endif
