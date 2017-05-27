
#include "SuffStat.hpp"
#include "BranchArray.hpp"

class BranchIIDExpo	{

	public:

	double* branchlength;
	double lambda;
};

class BranchIIDGamma : public BranchArray<double>	{

	BranchIIDGamma(const Tree* intree, double* inshape, double* inscale) : BranchArray(intree) {

		shape = inshape;
		scale = inscale;
		array[0] = 0;
		for (int j=1; j<GetNbranch(); j++)	{
			double tmp = Random::Gamma(*shape,*scale);
			array[j] = new double(tmp);
		}
	}
 
	double GetLogProb()	{

		double total = 0;
		for (int j=1; j<GetNbranch(); j++)	{
			total += GetFastLogProb(j);
		}
		total += (GetNbranch() - 1) * ((*shape) * log(*scale) - Random::logGamma(*shape))
		return total;
	}

	double GetFastLogProb(int branchindex)	{
		return (*shape - 1) * log(*array[branchindex]) - (*scale) * *array[branchindex];
	}

	double GetLogProb(int branchindex)	{
		return (*shape) * log(*scale) - Random::logGamma(*shape) + (*shape - 1) * log(*array[branchindex]) - (*scale) * *array[branchindex];
	}

	double GetTotal()	{
		double total = 0;
		for (int j=1; j<GetNbranch(); j++)	{
			total += *array[branchindex];
		}
	}
};
