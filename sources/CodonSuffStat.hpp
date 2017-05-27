
#ifndef CODONSUFFSTAT_H
#define CODONSUFFSTAT_H

#include "SuffStat.hpp"
#include "CodonStateSpace.hpp"

class NucPathSuffStat : public PathSuffStat	{

	public:
	NucPathSuffStat() : PathSuffStat() {}
	~NucPathSuffStat() {}

	// assumes pathsuffstat is 61x61
	// collect the 4x4 path suff stat out of codonpathsuffstat
	void AddSuffStat(CodonStateSpace* codonstatespace, PathSuffStat* codonpathsuffstat);
};

class OmegaSuffStat : public PoissonSuffStat	{

	public:

	OmegaSuffStat() {}
	~OmegaSuffStat() {}

	// assumes pathsuffstat is 61x61
	// tease out syn and non-syn substitutions and sum up count and beta stats  
	void AddSuffStat(CodonStateSpace* codonstatespace, PathSuffStat* pathsuffstat);

	// summing over all entries of an array
	void AddSuffStat(CodonStateSpace* codonstatespace, PathSuffStatArray* pathsuffstatarray)	{
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
