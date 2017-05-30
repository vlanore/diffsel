
#ifndef CODONSUFFSTAT_H
#define CODONSUFFSTAT_H

#include "PathSuffStat.hpp"
#include "CodonSubMatrixArray.hpp"
#include <typeinfo>

/*
class NucPathSuffStat : public PathSuffStat	{

	public:
	NucPathSuffStat() : PathSuffStat() {}
	~NucPathSuffStat() {}

	// assumes pathsuffstat is 61x61
	// collect the 4x4 path suff stat out of codonpathsuffstat
	void AddSuffStat(const MGCodonSubMatrix& codonstatespace, const PathSuffStat& codonpathsuffstat);
	void AddSuffStat(const MGCodonSubMatrix& codonstatespace, const PathSuffStatArray& codonpathsuffstat);
};
*/

class OmegaSuffStat : public PoissonSuffStat	{

	public:

	OmegaSuffStat() {}
	~OmegaSuffStat() {}

	// assumes pathsuffstat is 61x61
	// tease out syn and non-syn substitutions and sum up count and beta stats  
	void AddSuffStat(const MGOmegaCodonSubMatrix& codonsubmatrix, const PathSuffStat& pathsuffstat)	{
		pathsuffstat.AddOmegaSuffStat(*this,codonsubmatrix);
	}

	// summing over all entries of an array
	void AddSuffStat(const MGOmegaCodonSubMatrix& codonsubmatrix, const PathSuffStatArray& pathsuffstatarray)	{
		for (int i=0; i<pathsuffstatarray.GetSize(); i++)	{
			AddSuffStat(codonsubmatrix,pathsuffstatarray.GetVal(i));
		}
	}

};

class OmegaSuffStatArray : public PoissonSuffStatArray	{

	public:

	OmegaSuffStatArray(int insize) : PoissonSuffStatArray(insize) {}
	~OmegaSuffStatArray() {}

    OmegaSuffStat& GetOmegaSuffStat(int i)  {
        try {
            return dynamic_cast<OmegaSuffStat&> ((*this)[i]);
        }
        catch(std::bad_cast exp)    {
            std::cerr << "in OmegaSuffStatArray: bad cast exception\n";
            exit(1);
        }
    }

	void AddSuffStat(const MGOmegaHeterogeneousCodonSubMatrixArray& codonsubmatrixarray, const PathSuffStatArray& pathsuffstatarray)	{
		for (int i=0; i<pathsuffstatarray.GetSize(); i++)	{
            pathsuffstatarray.GetVal(i).AddOmegaSuffStat((*this)[i],*codonsubmatrixarray.GetMGOmegaCodonSubMatrix(i));
		}
	}
};

/*
class BranchOmegaSuffStatArray : public BranchPoissonSuffStatArray	{

	public:
	void AddSuffStat(CodonStateSpace* codonstatespace, BranchPathSuffStatArray* branchpathsuffstatarray)	{
		for (int i=0; i<GetNbranch(); i++)	{
			GetOmegaSuffStat(i).AddSuffStat(codonstatespace,branchpathsuffstatarray->GetVal(i));
		}
	}
};
*/

#endif
