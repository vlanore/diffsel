
#ifndef POISSONSUFFSTAT_H
#define POISSONSUFFSTAT_H

#include "SuffStat.hpp"
#include "Array.hpp"
#include <cmath>

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

class BranchPoissonSuffStatArray : public virtual PoissonSuffStatArray, public virtual BranchArray<PoissonSuffStat>  {

	public:

	BranchPoissonSuffStatArray(const Tree* intree) : Array<PoissonSuffStat>(intree->GetNbranch()), PoissonSuffStatArray(intree->GetNbranch()), BranchArray<PoissonSuffStat>(intree) {}
	~BranchPoissonSuffStatArray() {}
};

#endif
