
#ifndef POISSONSUFFSTAT_H
#define POISSONSUFFSTAT_H

#include "SuffStat.hpp"
#include "Array.hpp"

class PoissonSuffStat : public SuffStat	{

	public:

	PoissonSuffStat() {}
	~PoissonSuffStat() {}

	void Clear()	{
		count = beta = 0;
	}

	void IncrementCount()	{
		count++;
	}

	void AddCount(int in)	{
		count += in;
	}

	void AddBeta(double in)	{
		beta += in;
	}

	int GetCount() const {
		return count;
	}

	double GetBeta() const {
		return beta;
	}

	double GetLogProb(double rate) const {
		return count*log(rate) - beta*rate;
	}

	double GetMarginalLogProb(double shape, double scale) const {
		return shape*log(scale) - Random::logGamma(shape) - (shape + count)*log(scale + beta) + Random::logGamma(shape + count);
	}

	private:

	int count;
	double beta;
};

class PoissonSuffStatArray : public SimpleArray<PoissonSuffStat>	{

	public:

	PoissonSuffStatArray(int insize) : Array<PoissonSuffStat>(insize), SimpleArray<PoissonSuffStat>(insize) {}
	~PoissonSuffStatArray() {}

	void Clear()	{
		for (int i=0; i<GetSize(); i++)	{
			GetVal(i).Clear();
		}
	}

	double GetLogProb(const Array<double>* ratearray) const{
		double total = 0;
		for (int i=0; i<GetSize(); i++)	{
			total += GetVal(i).GetLogProb(ratearray->GetVal(i));
		}
		return total;
	}

	double GetMarginalLogProb(double shape, double scale)	const {
		double total = 0;
		/*
		for (int i=0; i<GetSize(); i++)	{
			total += GetVal(i).GetMarginalLogProb(shape,scale);
		}
		*/
		// factoring out prior factor
		for (int i=0; i<GetSize(); i++)	{
			int count = GetVal(i).GetCount();
			double beta = GetVal(i).GetBeta();
			total += -(shape+count)*log(scale+beta) + Random::logGamma(shape+count);
		}
		total += GetSize() * (shape*log(scale) - Random::logGamma(shape));
		return total;
	}
};

class BranchPoissonSuffStatArray : public SimpleBranchArray<PoissonSuffStat>	{

	public:

	BranchPoissonSuffStatArray(const Tree* intree) : Array<PoissonSuffStat>(intree->GetNbranch()), SimpleArray<PoissonSuffStat>(intree->GetNbranch()), BranchArray<PoissonSuffStat>(intree), SimpleBranchArray<PoissonSuffStat>(intree) {}
	~BranchPoissonSuffStatArray() {}
	
	void Clear()	{
		for (int i=0; i<GetNbranch(); i++)	{
			GetBranchVal(i).Clear();
		}
	}

	double GetLogProb(const BranchArray<double>* branchratearray) const {
		double total = 0;
		for (int i=0; i<GetNbranch(); i++)	{
			total += GetBranchVal(i).GetLogProb(branchratearray->GetBranchVal(i));
		}
		return total;
	}

	double GetMarginalLogProb(double shape, double scale) const {
		double total = 0;
		/*
		for (int i=0; i<GetNbranch(); i++)	{
			total += GetBranchVal(i).GetMarginalLogProb(shape,scale);
		}
		*/
		// factoring out prior factor
		for (int i=0; i<GetNbranch(); i++)	{
			int count = GetBranchVal(i).GetCount();
			double beta = GetBranchVal(i).GetBeta();
			total += -(shape+count)*log(scale+beta) + Random::logGamma(shape+count);
		}
		total += GetNbranch() * (shape*log(scale) - Random::logGamma(shape));
		return total;
	}
};

#endif
