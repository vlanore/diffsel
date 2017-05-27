
#ifndef POISSONSUFFSTAT_H
#define POISSONSUFFSTAT_H

#include "SuffStat.hpp"
#include "Array.hpp"
#include "BranchArray.hpp"
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

	void AddSuffStat(int incount, double inbeta)	{
		count += incount;
		beta += inbeta;
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

	PoissonSuffStatArray(int insize) : SimpleArray<PoissonSuffStat>(insize) {}
	~PoissonSuffStatArray() {}

	void Clear()	{
		for (int i=0; i<GetSize(); i++)	{
			(*this)[i].Clear();
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

class PoissonSuffStatBranchArray : public SimpleBranchArray<PoissonSuffStat>	{

	public:

	PoissonSuffStatBranchArray(const Tree* intree) : SimpleBranchArray<PoissonSuffStat>(intree) {}
	~PoissonSuffStatBranchArray() {}

	void Clear()	{
		for (int i=0; i<GetNbranch(); i++)	{
			(*this)[i].Clear();
		}
	}

	double GetLogProb(const ConstBranchArray<double>* ratearray) const{
		double total = 0;
		for (int i=0; i<GetNbranch(); i++)	{
			total += GetVal(i).GetLogProb(ratearray->GetVal(i));
		}
		return total;
	}

	double GetMarginalLogProb(double shape, double scale)	const {
		double total = 0;
		for (int i=0; i<GetNbranch(); i++)	{
			int count = GetVal(i).GetCount();
			double beta = GetVal(i).GetBeta();
			total += -(shape+count)*log(scale+beta) + Random::logGamma(shape+count);
		}
		total += GetNbranch() * (shape*log(scale) - Random::logGamma(shape));
		return total;
	}
};

#endif
