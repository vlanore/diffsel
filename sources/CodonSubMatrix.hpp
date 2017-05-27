#ifndef CODONSUBMATRIX_H
#define CODONSUBMATRIX_H

#include "CodonStateSpace.hpp"
#include "SubMatrix.hpp"

const double omegamin = 1e-10;

// a general class representing all codon matrices
// this is still an abstract class
class CodonSubMatrix : public virtual SubMatrix {
  public:
    CodonSubMatrix(const CodonStateSpace *instatespace, bool innormalise)
        : SubMatrix(instatespace->GetNstate(), innormalise), statespace(instatespace) {}

    const CodonStateSpace *GetCodonStateSpace() const { return statespace; }

    bool Synonymous(int codon1, int codon2) const { return statespace->Synonymous(codon1, codon2); }
    int GetCodonPosition(int pos, int codon) const { return statespace->GetCodonPosition(pos, codon); }
    int GetDifferingPosition(int codon1, int codon2) const {
        return statespace->GetDifferingPosition(codon1, codon2);
    }

  protected:
    const CodonStateSpace *statespace;
};

// most codon matrices rely on a mutation process at the level of nucleotides
// thus we create a NucCodonSubMatrix, which takes a substitution matrix
// representing the nucleotide
// 4x4 mutation process
// this is still an abstract class
class NucCodonSubMatrix : public virtual CodonSubMatrix {
  public:
    NucCodonSubMatrix(const CodonStateSpace *instatespace, const SubMatrix *inNucMatrix, bool innormalise)
        : SubMatrix(instatespace->GetNstate(), innormalise),
          CodonSubMatrix(instatespace, innormalise) {
        SetNucMatrix(inNucMatrix);
    }

    const SubMatrix *GetNucMatrix() const { return NucMatrix; }

  protected:
    void SetNucMatrix(const SubMatrix *inmatrix) {
        NucMatrix = inmatrix;
        if (NucMatrix->GetNstate() != Nnuc) {
            std::cerr << "error in CodonSubMatrix: underyling mutation process "
                         "should be a 4x4 matrix\n";
            throw;
        }
    }

    const SubMatrix *NucMatrix;
};

// The Muse and Gaut codon substitution process
// The simplest codon model based on a pure nucleotide mutation process (with
// stops excluded)
// look at how ComputeArray and ComputeStationary are implemented in
// CodonSubMatrix.cpp
class MGCodonSubMatrix : public NucCodonSubMatrix {
  public:
    MGCodonSubMatrix(const CodonStateSpace *instatespace, const SubMatrix *inNucMatrix,
                     bool innormalise = false)
        : SubMatrix(instatespace->GetNstate(), innormalise),
          CodonSubMatrix(instatespace, innormalise),
          NucCodonSubMatrix(instatespace, inNucMatrix, innormalise) {
	/*
        synnucarray = new double *[Nnuc];
        nonsynnucarray = new double *[Nnuc];
        for (int i = 0; i < Nnuc; i++) {
            synnucarray[i] = new double[Nnuc];
            nonsynnucarray[i] = new double[Nnuc];
            for (int j = 0; j < Nnuc; j++) {
                synnucarray[i][j] = 0;
                nonsynnucarray[i][j] = 0;
            }
        }
        nucflag = false;
	*/
    }

    ~MGCodonSubMatrix() override {
	/*
        for (int i = 0; i < Nnuc; i++) {
            delete[] synnucarray[i];
            delete[] nonsynnucarray[i];
        }
	*/
    }

    /*
    const double **GetSynNucArray() const {
        if (!nucflag) {
            ComputeNucArrays();
            nucflag = true;
        }
        return synnucarray;
    }

    const double **GetNonSynNucArray() const {
        if (!nucflag) {
            ComputeNucArrays();
            nucflag = true;
        }
        return nonsynnucarray;
    }

    double GetNormStat() {
        if (!nucflag) {
            ComputeNucArrays();
            nucflag = true;
        }
        return 1 - stopstat;
    }
    */

    void CorruptMatrix() override {
        // nucflag = false;
        SubMatrix::CorruptMatrix();
    }

  protected:
    // look at how ComputeArray and ComputeStationary are implemented in
    // CodonSubMatrix.cpp
    void ComputeArray(int i) const override;
    void ComputeStationary() const override;

    /*
    virtual void ComputeNucArrays();

    double **synnucarray;
    double **nonsynnucarray;
    double stopstat;
    bool nucflag;
    */
};

// The Muse and Gaut codon substitution process
// with an omgea = dN/dS parameter
// look at how ComputeArray and ComputeStationary are implemented in
// CodonSubMatrix.cpp
class MGOmegaCodonSubMatrix : public MGCodonSubMatrix {
  public:
    MGOmegaCodonSubMatrix(const CodonStateSpace *instatespace, const SubMatrix *inNucMatrix, double inomega,
                          bool innormalise = false)
        : SubMatrix(instatespace->GetNstate(), innormalise),
          CodonSubMatrix(instatespace, innormalise),
          MGCodonSubMatrix(instatespace, inNucMatrix, innormalise),
          omega(inomega) {}

    double GetOmega() const { return omega + omegamin; }

    void SetOmega(double inomega) { omega = inomega; CorruptMatrix();}

  protected:
    // look at how ComputeArray and ComputeStationary are implemented in
    // CodonSubMatrix.cpp
    void ComputeArray(int i) const override;

    // void ComputeNucArrays() const override;

    void ToStream(std::ostream &os) const override {
        os << "Omega : " << omega << '\n';
        os << "nuc matrix\n";
        GetNucMatrix()->ToStream(os);
        os << '\n';
        SubMatrix::ToStream(os);
    }

    // data members

    double omega;
};

// this class implements the projection of a 61x61 codon substitution process
// onto a 20x20 amino-acid replacement process
// according to the formula:
//
// R_{ab} = [ \sum \pi_i Q_{ij} ] [ \sum_i \pi_i ]
//
// where a,b runs over amino-acids
// and i (resp. j) over all codons encoding for amino acid a (resp. b)
//
/*
class AminoAcidReducedCodonSubMatrix : public virtual SubMatrix {
  public:
    AminoAcidReducedCodonSubMatrix(CodonSubMatrix *incodonmatrix, bool innormalise = false)
        : SubMatrix(Naa, innormalise), codonmatrix(incodonmatrix) {
        aastatespace = new ProteinStateSpace();
    }

    CodonSubMatrix *GetCodonMatrix() { return codonmatrix; }
    CodonSubMatrix *GetCodonSubMatrix() { return codonmatrix; }
    CodonStateSpace *GetCodonStateSpace() { return GetCodonMatrix()->GetCodonStateSpace(); }
    ProteinStateSpace *GetProteinStateSpace() { return aastatespace; }

  protected:
    void SetCodonMatrix(CodonSubMatrix *incodonmatrix) { codonmatrix = incodonmatrix; }

    void ComputeArray(int a) override;
    void ComputeStationary() override;

    CodonSubMatrix *codonmatrix;
    ProteinStateSpace *aastatespace;
};
*/

#endif
