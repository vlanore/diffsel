
#ifndef CODONSUFFSTAT_H
#define CODONSUFFSTAT_H

#include "PoissonSuffStat.hpp"
#include "PathSuffStat.hpp"
#include "CodonSubMatrix.hpp"

class NucPathSuffStat : public PathSuffStat	{

	public:
	NucPathSuffStat() : PathSuffStat() {}
	~NucPathSuffStat() {}

	// assumes pathsuffstat is 61x61
	// collect the 4x4 path suff stat out of codonpathsuffstat
	void AddSuffStat(const MGCodonSubMatrix& codonstatespace, const PathSuffStat& codonpathsuffstat);
	void AddSuffStat(const MGCodonSubMatrix& codonstatespace, const PathSuffStatArray& codonpathsuffstat);
};

class OmegaSuffStat : public PoissonSuffStat	{

	public:

	OmegaSuffStat() {}
	~OmegaSuffStat() {}

	// assumes pathsuffstat is 61x61
	// tease out syn and non-syn substitutions and sum up count and beta stats  
	void AddSuffStat(const MGOmegaCodonSubMatrix& codonsubmatrix, const PathSuffStat& pathsuffstat);

	// summing over all entries of an array
	void AddSuffStat(const MGOmegaCodonSubMatrix& codonsubmatrix, const PathSuffStatArray& pathsuffstatarray)	{
		for (int i=0; i<GetSize(); i++)	{
			AddSuffStat(codonstatespace,pathsuffstatarray->GetVal(i));
		}
	}
};

class OmegaSuffStatArray : public PoissonSuffStatArray	{

	public:

	void AddSuffStat(CodonStateSpace* codonstatespace, PathSuffStatArray* pathsuffstatarray)	{
		for (int i=0; i<GetSize(); i++)	{
			GetVal(i).AddSuffStat(codonstatespace,pathsuffstatarray->GetVal(i));
		}
	}
};

class BranchOmegaSuffStatArray : public BranchPoissonSuffStatArray	{

	public:

	void AddSuffStat(CodonStateSpace* codonstatespace, BranchPathSuffStatArray* branchpathsuffstatarray)	{
		for (int i=0; i<GetNbranch(); i++)	{
			GetVal(i).AddSuffStat(codonstatespace,branchpathsuffstatarray->GetVal(i));
		}
	}
};

#endif
