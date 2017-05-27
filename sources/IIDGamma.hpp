
#ifndef IIDGAMMA_H
#define IIDGAMMA_H

#include "Array.hpp"
#include "BranchArray.hpp"
#include "Random.hpp"
#include "PoissonSuffStat.hpp"

class GammaSuffStat : public SuffStat	{

	public:
	GammaSuffStat() {}
	~GammaSuffStat() {}

	void Clear()	{
		sum = 0;
		sumlog = 0;
		n = 0;
	}

	void AddSuffStat(double x, double logx, int c = 1)	{
		sum += x;
		sumlog += logx;
		n += c;
	}

	double GetLogProb(double shape, double scale)	{
		return n*(shape*log(scale) - Random::logGamma(shape)) + (shape-1)*sumlog - scale*sum;
	}
	
	private:

	double sum;
	double sumlog;
	int n;
};

class IIDGamma: public virtual SimpleArray<double>	{

	public: 

	IIDGamma(int insize, double inshape, double inscale) : SimpleArray<double>(insize), shape(inshape), scale(inscale)	{
		Sample();
	}

	~IIDGamma() {}

	double GetShape() const {return shape;}
	double GetScale() const {return scale;}

	void SetShape(double inshape)	{
		shape = inshape;
	}

	void SetScale(double inscale)	{
		scale = inscale;
	}

	void Sample()	{
		for (int i=0; i<GetSize(); i++)	{
			(*this)[i] = Random::Gamma(shape,scale);
		}
	}

	void GibbsResample(const PoissonSuffStatArray* suffstatarray)	{
		for (int i=0; i<GetSize(); i++)	{
			const PoissonSuffStat& suffstat = suffstatarray->GetVal(i);
			(*this)[i] = Random::Gamma(shape + suffstat.GetCount(), scale + suffstat.GetBeta());
		}
	}

	double GetLogProb()	{
		double total = 0;
		for (int i=0; i<GetSize(); i++)	{
			total += GetLogProb(i);
		}
		return total;
	}

	double GetLogProb(int index)	{
		return shape * log(scale) - Random::logGamma(shape) + (shape-1)*log((*this)[index]) - scale*(*this)[index];
	}

	void AddSuffStat(GammaSuffStat& suffstat)	{
		for (int i=0; i<GetSize(); i++)	{
			suffstat.AddSuffStat((*this)[i],log((*this)[i]));
		}
	}

	protected:
	double shape;
	double scale;
};

class BranchIIDGamma: public virtual SimpleBranchArray<double>	{

	public: 

	BranchIIDGamma(const Tree* intree, double inshape, double inscale) : SimpleBranchArray<double>(intree), shape(inshape), scale(inscale)	{
		Sample();
	}

	~BranchIIDGamma() {}

	double GetShape() const {return shape;}
	double GetScale() const {return scale;}

	void SetShape(double inshape)	{
		shape = inshape;
	}

	void SetScale(double inscale)	{
		scale = inscale;
	}

	void Sample()	{
		for (int i=0; i<GetNbranch(); i++)	{
			(*this)[i] = Random::Gamma(shape,scale);
		}
	}

	void GibbsResample(const PoissonSuffStatBranchArray* suffstatarray)	{
		for (int i=0; i<GetNbranch(); i++)	{
			const PoissonSuffStat& suffstat = suffstatarray->GetVal(i);
			(*this)[i] = Random::Gamma(shape + suffstat.GetCount(), scale + suffstat.GetBeta());
		}
	}

	double GetLogProb()	{
		double total = 0;
		for (int i=0; i<GetNbranch(); i++)	{
			total += GetLogProb(i);
		}
		return total;
	}

	double GetLogProb(int index)	{
		return shape * log(scale) - Random::logGamma(shape) + (shape-1)*log((*this)[index]) - scale*(*this)[index];
	}

	void AddSuffStat(GammaSuffStat& suffstat)	{
		for (int i=0; i<GetNbranch(); i++)	{
			suffstat.AddSuffStat((*this)[i],log((*this)[i]));
		}
	}

	protected:
	double shape;
	double scale;
};

#endif
