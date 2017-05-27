#ifndef GTRSUBMATRIX_H
#define GTRSUBMATRIX_H

#include "BiologicalSequences.hpp"  //FIXME only used for Naa (const int)
#include "SubMatrix.hpp"

class GTRSubMatrix : public virtual SubMatrix {
  public:
    GTRSubMatrix(int inNstate, const std::vector<double>& rr, const std::vector<double>& stat, bool innormalise = false);
    ~GTRSubMatrix() override = default;

    int GetNRelativeRate() const { return Nrr; }
    double RelativeRate(int i, int j) const { return mRelativeRate[rrindex(i, j, GetNstate())]; }

    static int rrindex(int i, int j, int nstate) {
        return (i < j) ? (2 * nstate - i - 1) * i / 2 + j - i - 1
                       : (2 * nstate - j - 1) * j / 2 + i - j - 1;
    }

    // make a copy of the entries (not of the pointer)
    void CopyStationary(const std::vector<double>& instat);

  protected:

  protected:
    void ComputeArray(int i) const override;
    void ComputeStationary() const override {}

    const std::vector<double>&  mRelativeRate;
    int Nrr;
};

#endif
