
#ifndef IIDGAMMA_H
#define IIDGAMMA_H

#include "Array.hpp"
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


class IIDGamma: public virtual Array<double>	{

	public: 

	IIDGamma(int insize, double inshape, double inscale) : Array(insize), array(insize), shape(inshape), scale(inscale)	{
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
			array[i] = Random::Gamma(shape,scale);
		}
	}

	void GibbsResample(const PoissonSuffStatArray* suffstatarray)	{
		for (int i=0; i<GetSize(); i++)	{
			const PoissonSuffStat& suffstat = suffstatarray->GetVal(i);
			array[i] = Random::Gamma(shape + suffstat.GetCount(), scale + suffstat.GetBeta());
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
		return shape * log(scale) - Random::logGamma(shape) + (shape-1)*log(array[index]) - scale*array[index];
	}

	void AddSuffStat(GammaSuffStat& suffstat)	{
		for (int i=0; i<GetSize(); i++)	{
			suffstat.AddSuffStat(array[i],log(array[i]));
		}
	}

	const double& GetVal(int index) const {
		return array[index];
	}
	double& GetVal(int index) {
		return array[index];
	}

	protected:
	vector<double> array;
	double shape;
	double scale;
};

class BranchIIDGamma: public virtual BranchArray<double>, public IIDGamma	{

	public: 

	BranchIIDGamma(const Tree* intree, double inshape, double inscale) : BranchArray<double>(intree), Array<double>(intree->GetNbranch()), IIDGamma(intree->GetNbranch(),inshape,inscale) {}
	~BranchIIDGamma() {}
};

#endif
