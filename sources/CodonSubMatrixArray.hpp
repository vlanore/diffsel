
#ifndef CODONMATRIXARRAY_H
#define CODONMATRIXARRAY_H

#include "GTRSubMatrix.hpp"
#include "CodonSubMatrix.hpp"
#include "Array.hpp"

class MGOmegaHeterogeneousCodonSubMatrixArray : public Array<SubMatrix>	{

	public:
	MGOmegaHeterogeneousCodonSubMatrixArray(const CodonStateSpace* incodonstatespace, const GTRSubMatrix* innucmatrix, const Array<double>* inomegaarray) : Array<SubMatrix>(inomegaarray->GetSize()), codonstatespace(incodonstatespace), nucmatrix(innucmatrix), omegaarray(inomegaarray), matrixarray(inomegaarray->GetSize())   {
        Create();
	}

	~MGOmegaHeterogeneousCodonSubMatrixArray()	{
        Delete();
	}
		
	MGOmegaCodonSubMatrix* GetMGOmegaCodonSubMatrix(int site) const {
		return dynamic_cast<MGOmegaCodonSubMatrix*>(matrixarray[site]);
	}

    const SubMatrix& GetVal(int i) const {return *matrixarray[i];}
    SubMatrix& operator[](int i) {return *matrixarray[i];}

	const GTRSubMatrix& GetNucMatrix() const {return *nucmatrix;}

    void UpdateCodonMatrices()  {
		for (int i=0; i<GetSize(); i++)	{
            GetMGOmegaCodonSubMatrix(i)->SetOmega(omegaarray->GetVal(i));
            GetMGOmegaCodonSubMatrix(i)->CorruptMatrix();
		}
    }

	private:

    void Create()   {
		for (int i=0; i<GetSize(); i++)	{
			matrixarray[i] = new MGOmegaCodonSubMatrix(codonstatespace,nucmatrix,omegaarray->GetVal(i));
		}
    }

    void Delete()   {
		for (int i=0; i<GetSize(); i++)	{
			delete matrixarray[i];
		}
    }
        

	const CodonStateSpace* codonstatespace;
	const GTRSubMatrix* nucmatrix;
	const Array<double>* omegaarray;
    vector<MGOmegaCodonSubMatrix*> matrixarray;
};

/*
class MGOmegaHeterogeneousCodonSubMatrixArray : public SimpleArray<SubMatrix*>	{

	public:
	MGOmegaHeterogeneousCodonSubMatrixArray(const CodonStateSpace* incodonstatespace, const GTRSubMatrix* innucmatrix, const Array<double>* inomegaarray) : SimpleArray<SubMatrix*>(inomegaarray->GetSize()), codonstatespace(incodonstatespace), nucmatrix(innucmatrix), omegaarray(inomegaarray)   {
        Create();
	}

	~MGOmegaHeterogeneousCodonSubMatrixArray()	{
        Delete();
	}
		
	MGOmegaCodonSubMatrix* GetMGOmegaCodonSubMatrix(int site) const {
		return dynamic_cast<MGOmegaCodonSubMatrix*>(array[site]);
	}

	const GTRSubMatrix& GetNucMatrix() const {return *nucmatrix;}

    void UpdateCodonMatrices()  {
		for (int i=0; i<GetSize(); i++)	{
            GetMGOmegaCodonSubMatrix(i)->SetOmega(omegaarray->GetVal(i));
            GetMGOmegaCodonSubMatrix(i)->CorruptMatrix();
		}
    }

	private:

    void Create()   {
		for (int i=0; i<GetSize(); i++)	{
			array[i] = new MGOmegaCodonSubMatrix(codonstatespace,nucmatrix,omegaarray->GetVal(i));
		}
    }

    void Delete()   {
		for (int i=0; i<GetSize(); i++)	{
			delete array[i];
		}
    }
        

	const CodonStateSpace* codonstatespace;
	const GTRSubMatrix* nucmatrix;
	const Array<double>* omegaarray;
};
*/

#endif
