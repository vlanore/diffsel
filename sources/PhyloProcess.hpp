#ifndef PHYLOPROCESS_H
#define PHYLOPROCESS_H

#include <map>
#include "BranchSitePath.hpp"
#include "SequenceAlignment.hpp"
#include "Chrono.hpp"
#include "SubMatrix.hpp"
#include "Tree.hpp"

// PhyloProcess is a dispatcher:
// its responsibility is to create a random branch/site path
// for each branch/site pair

class PhyloProcess	{

public:

	// generic constructor
	PhyloProcess(const Tree* intree, const SequenceAlignment* indata, const BranchArray<double>* inbranchlength, const Array<double>* insiterate, const BranchSiteArray<SubMatrix>* insubmatrixarray, const Array<SubMatrix>* inrootsubmatrixarray);

	// branch and site homogeneous
	PhyloProcess(const Tree* intree, const SequenceAlignment* indata, const BranchArray<double>* inbranchlength, const Array<double>* insiterate, SubMatrix* insubmatrix);

	// branch homogeneous, site heterogeneous
	PhyloProcess(const Tree* intree, const SequenceAlignment* indata, const BranchArray<double>* inbranchlength, const Array<double>* insiterate, Array<SubMatrix>* insubmatrixarray);

	~PhyloProcess();

	// accessors

	double GetBranchLength(int branch) const {
		return branchlength->GetBranchVal(branch);
	}

	double GetSiteRate(int site) const {
		if (! siterate)	{
			return 1.0;
		}
		return siterate->GetVal(site);
	}

	const SubMatrix& GetSubMatrix(int branch, int site) const {
		return submatrixarray->GetVal(branch,site);
	}

	const double* GetRootFreq(int site) const {
		return rootsubmatrixarray->GetVal(site).GetStationary();
	}

	double SiteLogLikelihood(int site);
	double FastSiteLogLikelihood(int site);

	const StateSpace* GetStateSpace() const { return data->GetStateSpace(); }
	const TaxonSet* GetTaxonSet() const { return data->GetTaxonSet(); }

	int GetNsite() const { return data->GetNsite(); }
	int GetNtaxa() const { return data->GetNtaxa(); }

	int GetNstate() const { return Nstate;}

	const SequenceAlignment *GetData() const { return data; }
	int GetData(int taxon, int site) const { return data->GetState(taxon, site); }

	const Tree* GetTree() const { return tree;}
	const Link* GetRoot() const { return GetTree()->GetRoot(); }

	int GetMaxTrial() const { return maxtrial; }
	void SetMaxTrial(int i) { maxtrial = i; }

	void SetData(const SequenceAlignment *indata);
	void ClampData() { clampdata = true; }
	void UnclampData() { clampdata = false; }

	int& GetState(const Node *node, int site) { return statemap[node][site]; }

	void GetLeafData(SequenceAlignment *data);
	void RecursiveGetLeafData(const Link *from, SequenceAlignment *data);

	bool isDataCompatible(int taxon, int site, int state) const {
		return GetStateSpace()->isCompatible(GetData(taxon, site), state);
	}

	// probability, pruning, sampling
	double GetLogProb();                        // likelihood Felsenstein 1981
	double GetFastLogProb();                            // likelihood Felsenstein 1981
	double GetPathLogProb();                            // probability of the entire mapping
	double GetPathLogProb(const Link *from, int site);  // probability of the entire mapping

	double Move(double fraction);
	void DrawSites(double fraction); // draw a fraction of sites which will be resampled
	void ResampleSub();  // clamped Nielsen
	void ResampleSub(int site);

	// basic building blocks for compiling suffstats...

	void AddRootSuffStat(int site, PathSuffStat& suffstat);
	void AddPathSuffStat(const Link* link, int site, PathSuffStat& suffstat);
	void AddLengthSuffStat(const Link* link, int site, PoissonSuffStat& suffstat);
	// void AddPoissonSuffStat(const Link* link, int site, PoissonSuffStat& suffstat);

	// ... which are then used for looping over branches and sites

	// homogeneous across sites and branches
	void AddPathSuffStat(PathSuffStat& suffstat);
	void RecursiveAddPathSuffStat(const Link* from, PathSuffStat& suffstat);

	// heterogeneeous across sites, homogeneous across branches
	void AddPathSuffStat(PathSuffStatArray* suffstatarray);
	void RecursiveAddPathSuffStat(const Link* from, PathSuffStatArray* suffstatarray);

	// homogeneous across sites, heterogeneous across branches
	// void AddSuffStat(BranchPathSuffStatArray* branchsuffstatarray, PathSuffStat& rootsuffstat);

	// heterogeneous across sites and branches
	// void AddSuffStat(BranchSitePathSuffStatArray* branchsitesuffstatarray);

	// homogeneous across sites
	void AddLengthSuffStat(BranchPoissonSuffStatArray* branchlengthsuffstatarray);
	void RecursiveAddLengthSuffStat(const Link* from, BranchPoissonSuffStatArray* branchlengthsuffstatarray);

	// homogeneous across branches
	// void AddPoissonSuffStat(PoissonSuffStatArray* poissonsuffstatarray);

	void PostPredSample(bool rootprior = false);  // unclamped Nielsen
	void PostPredSample(int site, bool rootprior = false);
	// rootprior == true : root state drawn from stationary probability of the
	// process
	// rootprior == false: root state drawn from posterior distribution

	// various accessors

	bool isMissing(const Node *node, int site) const { return false; return missingmap[node][site]; }

	bool isMissing(const Link *link, int site) const {
		return false;
		return (missingmap[link->GetNode()][site] || missingmap[link->Out()->GetNode()][site]);
	}

	void CreateMissingMap();
	void RecursiveCreateMissingMap(const Link *from);
	bool FillMissingMap(const Link *from, int i);

	double *GetCondLikelihood(const Link *from) const { return condlmap[from]; }

	double GetPruningTime() const { return pruningchrono.GetTime(); }
	double GetResampleTime() const { return resamplechrono.GetTime(); }

	void Unfold();
	void Cleanup();

	protected:

	void RecursiveCreate(const Link *from);
	void RecursiveDelete(const Link *from);

	void RecursiveCreateTBL(const Link *from);
	void RecursiveDeleteTBL(const Link *from);

	void Pruning(const Link *from, int site);
	void ResampleSub(const Link *from, int site);
	void ResampleState();
	void ResampleState(int site);
	void PruningAncestral(const Link *from, int site);
	void PriorSample(const Link *from, int site, bool rootprior);
	void PriorSample();
	void RootPosteriorDraw(int site);

	// borrowed from phylobayes
	// where should that be?
	BranchSitePath* SamplePath(int stateup, int statedown, double time, double rate, const SubMatrix& matrix);
	BranchSitePath* SampleRootPath(int rootstate);
	BranchSitePath* ResampleAcceptReject(int maxtrial, int stateup, int statedown, double rate, double totaltime, const SubMatrix& matrix);
	BranchSitePath* ResampleUniformized(int stateup, int statedown, double rate, double totaltime, const SubMatrix& matrix);

	const Tree *tree;
	const SequenceAlignment *data;
	const BranchArray<double>* branchlength;
	const Array<double>* siterate;
	const BranchSiteArray<SubMatrix>* submatrixarray; 
	const Array<SubMatrix>* rootsubmatrixarray;
	bool allocsubmatrixarray;
	bool allocrootsubmatrixarray;

	int *sitearray;
	double *sitelnL;
	mutable std::map<const Link *, double *> condlmap;

	int Nstate;

	bool clampdata;

	BranchSitePath* GetPath(const Node* node, int site)	{
	    if (pathmap[node][site] == nullptr) {
		std::cerr << "error in phyloprocess::getpath: null path\n";
		exit(1);
	    }
	    return pathmap[node][site];
	}

	private:

	mutable std::map<const Node*, BranchSitePath **> pathmap;
	mutable std::map<const Node *, int *> statemap;
	mutable std::map<const Node *, bool *> missingmap;
	mutable std::map<const Node *, int> totmissingmap;

	int maxtrial;
	static const int unknown = -1;

	static const int DEFAULTMAXTRIAL = 100;

	mutable Chrono pruningchrono;
	mutable Chrono resamplechrono;
};

#endif  // PHYLOPROCESS_H
