#ifndef CODONSEQUENCEALIGNMENT_H
#define CODONSEQUENCEALIGNMENT_H

#include <cmath>
#include "CodonStateSpace.hpp"
#include "SequenceAlignment.hpp"

class CodonSequenceAlignment : public SequenceAlignment {
  public:
    CodonSequenceAlignment(SequenceAlignment *from, bool force_stops = false,
                           GeneticCodeType type = Universal);
    ~CodonSequenceAlignment() override = default;

    CodonStateSpace *GetCodonStateSpace() {
        // return static_cast<CodonStateSpace*>(statespace);
        return (CodonStateSpace *)(statespace);
    }

    void ToStream(std::ostream &os);
    void ToStream(std::ostream &os, int pos);
    void ToStreamFourFoldThird(std::ostream &os);
    void ToStreamFourFoldTriplet(std::ostream &os);
    void ToStreamFourFoldThirdwoCpG(std::ostream &os);

    void ToStreamRandomJackknife(std::ostream &os, double p);

  private:
    SequenceAlignment *DNAsource;
};

#endif
