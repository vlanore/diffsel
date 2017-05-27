
template<class T> Array	{

	public:
	Array(int size);
	virtual ~Array() {}

	int GetSize() {return size;}
	virtual const T& GetVal(int index) = 0;

	protected:
	int size;
};

template<class T> HomogeneousArray : public virtual Array<T>	{

	public:
	HomogeneousArray(int insize, const T* invalue) : Array<T>(size), value(invalue) {}
	~HomogeneousArray() {}

	virtual const T& GetVal(int index) const {return *value;}

	private:
	const T* value;
};

template<class T> BranchArray	{

	public:
	BranchArray<T>(const Tree* intree);
	virtual ~BranchArray<T> {}

	const Tree* GetTree() const {return tree;}
	int GetNbranch() const {return tree->GetNbranch();}

	virtual T& GetBranchVal(int branch) const = 0;

	protected:
	const Tree* tree;
};

template<class T> BranchHomogeneousArray : public virtual BranchArray<T> {

	public:

	BranchHomogeneousArray(const Tree* intree, const T* invalue) : BranchArray<T>(intree), value(invalue)	{}
	~BranchHomogeneousArray() {}

	const T& GetBranchVal(int branch) const {return *value;}

	private:
	const T* value;

	
};

// dangerous (depends on whether we opt for automatic or dynamic creation
template<class T> SimpleBranchArray : public virtual BranchArray<T> {

	public:

	SimpleBranchArray(const Tree* intree, const T& init) : BranchArray<T>(intree), array(intree->GetNbranch(),init) {}
	~SimpleBranchArray() {}
	private:
	vector<T> array;
};

template<class T> BranchPtrArray : public virtual BranchArray<T*>	{

	public:
	BranchPtrArray(const Tree* intree) : BranchArray<T*>(intree), array(intree->GetNbranch()) {
		Create();
	}

	virtual ~BranchPtrArray() {
		Delete();
	}

	void Create()	{
		RecursiveCreate(GetTree()->GetRoot());
	}

	void Delete()	{
		RecursiveDelete(GetTree()->GetRoot());
	}

	virtual void CreateBranchVal(const Link* link) = 0;
};

template<class T> BranchSiteArray	{

	public:
	BranchSiteArray(const Tree* intree, int insize) : tree(intree), size(insize) {}
	virtual ~BranchSiteArray {}

	const Tree* GetTree() const {return tree;}
	int GetNbranch() const {return tree->GetNbranch();}
	int GetSize() const {return size;}

	virtual const T& GetVal(int branch, int site) const = 0;

	protected:
	const Tree* tree;
	int size;
};


template<class T> HomogeneousBranchSiteArray : public virtual BranchSiteArray<T> {

	public:
	HomogeneousBranchSiteArray(const Tree* intree, int insize, const T* invalue) : BranchSiteArray<T>(intree,insize), value(invalue) {}
	~HomogeneousBranchSiteArray() {}

	const T& GetVal(int branch, int site) const {return *value;}

	private:
	const T* value;
};

template<class T> BranchHomogeneousSiteHeterogeneousArray : public virtual BranchSiteArray<T>	{

	public:
	BranchHomogeneousSiteHeterogeneousArray(const Tree* intree, const Array<T>* inarray) : BranchSiteArray<T>(intree,inarray->GetSize()), array(inarray) {}
	~BranchHomogeneousSiteHeterogeneousArray() {}

	const T& GetVal(int branch, int site) const {return array->GetVal(site);}

	private:
	const Array<T>* array;
};

template<class T> BranchHeterogeneousSiteHomogeneousArray : public virtual BranchSiteArray<T> {

	public
	BranchHeterogeneousSiteHomogeneousArray(const BranchArray<T>* inbrancharray, int insize) : BranchSiteArray<T>(inbrancharray->GetTree(),insize), brancharray(inbrancharray) {}
	~BranchHeterogeneousSiteHomogeneousArray() {}

	const T& GetVal(int branch, int site) const {return brancharray->GetBranchVal(branch);}

	private:
	const BranchArray<T>* brancharray;
};

class PhyloProcess	{

	double GetBranchLength(int branch);
	double GetSiteRate(int site);
	const SubMatrix& GetSubMatrix(int branch, int site);
	const Profile& GetRootFreq(int site);

	// basic building blocks for compiling suffstats...

	void AddRootSuffStat(int site, PathSuffStat& suffstat);
	void AddPathSuffStat(int branch, int site, PathSuffStat& suffstat);
	void AddLengthSuffStat(int branch, int site, RateSuffStat& suffstat);
	void AddRateSuffStat(int branch, int site, RateSuffStat& suffstat);

	// ... which are then used for looping over branches and sites

	// homogeneous across sites and branches
	void AddSuffStat(PathSuffStat& suffstat);

	// heterogeneeous across sites, homogeneous across branches
	void AddSuffStat(PathSuffStatArray* suffstatarray);

	// homogeneous across sites, heterogeneous across branches
	void AddSuffStat(BranchPathSuffStatArray* branchsuffstatarray);

	// heterogeneous across sites and branches
	void AddSuffStat(BranchSitePathSuffStatArray* branchsitesuffstatarray);

	// homogeneous across sites
	void AddLengthSuffStat(BranchRateSuffStatArray* branchlengthsuffstat);

	// homogeneous across branches
	void AddRateSuffStat(RateSuffStatArray* ratesuffstatarray);

	private:

	const BranchArray<double>* branchlength;
	const Array<double>* siterate;
	const BranchSiteArray<SubMatrix>* submatrix; 
	const Array<SubMatrix>* rootsubmatrix;
};

class PathSuffStat : public SuffStat	{

	friend class PhyloProcess;
	friend class BranchSitePath;

	public:

	PathSuffStat() {}
	~PathSuffStat() {}

	void Clear()	{
		rootcount.clear();
		paircount.clear();
		waitingtime.clear();
	}

	void AddSuffStat(PhyloProcess* phyloprocess)	{

		// loop over branches and sites
		// call phyloprocess->AddSuffStat(*this)
	}

	private:
	map<int,int> rootcount;
	map<pair<int,int>,int> paircount;
	map<int,double> waitingtime;
};

class PathSuffStatArray : public Array<PathSuffStat>	{

	public:

	PathSuffStatArray(int insize) : Array<PathSuffStat>(insize) {}
	~PathSuffStatArray();

	void Clear()	{
		for (int i=0; i<GetSize(); i++)	{
			GetVal(i).Clear();
		}
	}

	void AddSuffStat(PhyloProcess* phyloprocess)	{
};

class NucPathSuffStat : public PathSuffStat	{

	// assumes pathsuffstat is 61x61
	// collect the 4x4 path suff stat out of codonpathsuffstat
	void AddSuffStat(CodonStateSpace* codonstatespace, PathSuffStat* codonpathsuffstat);
};

class RateSuffStat : public SuffStat	{

	// count and beta
};

class OmegaSuffStat : public RateSuffStat	{

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

enum SuffStatType = {RATE,LENGTH};

class BranchRateSuffStatArray : public BranchArray<RateSuffStat>	{

	public:

	BranchRateSuffStatArray(const Tree* intree) : BranchArray<RateSuffStat>(intree) {}
	// BranchRateSuffStatArray(const BranchArray<double>* inbranchlength) : BranchArray<RateSuffStat>(inbranchlength->GetTree()), branchlength(inbranchlength) {}

	~BranchRateSuffStatArray() {}
	
	void ClearSuffStat();
	void AddSuffStat(PhyloProcess* inphyloprocess, SuffStatType type);

	double GetSuffStatLogProb(int branch, double length);

	double GetSuffStatCount(int branch);
	double GetSuffStatBeta(int branch);
};

class RateSuffStatArray : public Array<RateSuffStat>	{

};

class OmegaSuffStatArray : public RateSuffStatArray	{

	public:
	void AddSuffStat(CodonStateSpace* codonstatespace, PathSuffStatArray* pathsuffstatarray)	{
		for (int i=0; i<GetSize(); i++)	{
			GetVal(i).AddSuffStat(codonstatespace,pathsuffstatarray->GetVal(i));
		}
	}
};

class BranchOmegaSuffStatArray : public BranchRateSuffStatArray	{

	public:

	void AddSuffStat(CodonStateSpace* codonstatespace, BranchPathSuffStatArray* branchpathsuffstatarray)	{
		for (int i=0; i<GetNbranch(); i++)	{
			GetVal(i).AddSuffStat(codonstatespace,branchpathsuffstatarray->GetVal(i));
		}
	}
};

class IIDGamma: public Array<double>	{

	public: 

	IIDGamma(int insize, double inshape, double inscale) : array(size), shape(inshape), scale(inscale)	{
		Sample();
	}

	~IIDGamma() {}

	void SetShape(double inshape);
	void SetScale(double inscale);

	void Sample();

	void GibbsResample(const RateSuffStatArray* ratesuffstatarray)	{
		for (int i=0; i<GetSize(); i++)	{
			const RateSuffStat& suffstat = ratesuffstatarray->GetVal(i);
			array[i] = Random::Gamma(shape + suffstat.GetCount(), scale + suffstat.GetBeta());
		}
	}

	double GetLogProb();
	double GetLogProb(int index);
	double GetSuffStatLogProb(SuffStat& suffstat);

	void AddSuffStat(SomeSuffStatHere& suffstat);

	const double& GetVal(int index) const {
		return array[index];
	}

	private:
	vector<double> array;
	double shape;
	double scale;
};

class BranchIIDGamma: public virtual BranchArray<double>, public IIDGamma	{

	public: 

	BranchIIDGamma(const Tree* intree, double inshape, double inscale) : BranchArray<double>(intree), IIDGamma(intree->GetNbranch(),inshape,inscale) {}
	~BranchIIDGamma() {}

	GibbsResample(BranchRateSuffStatArray* branchratesuffstatarray)	{
		for (int i=0; i<GetNbranch(); i++)	{
			const RateSuffStat& suffstat = branchratesuffstatarray->GetVal(i);
			array[i] = Random::Gamma(shape + suffstat.GetCount(), scale + suffstat.GetBeta());
		}
	}

};

// pure interface
class SuffStat	{

};

class MGOmegaHeterogeneousCodonSubMatrixArray : public Array<SubMatrix>	{

	public:
	MGOmegaHeterogeneousCodonSubMatrixArray(CodonStateSpace* incodonstatespace, GTRSubMatrix* innucmatrix, const Array<double>* inomegaarray) : codonstatespace(incodonstatespace), nucmatrix(innucmatrix), omegaarray(inomegaarray), matrixarray(inomegaarray->GetSize()) {
		for (int i=0; i<matrixarray.size(); i++)	{
			matrixarray[i] = new MGOmegaCodonSubMatrix(nucmatrix,omegarray->GetVal(i));
		}
	}

	~MGOmegaHeterogeneousCodonSubMatrixArray()	{
		for (int i=0; i<matrixarray.size(); i++)	{
			delete matrixarray[i];
		}
	}
		
	int GetSize() {return matrixarray.size();}

	MGOmegaCodonSubMatrix& GetMGOmegaCodonSubMatrix(int site) const {
		return *matrixarray[site];
	}

	const GTRSubMatrix& GetNucMatrix() const {return *nucmatrix;}

	const SubMatrix& GetVal(int site) const {return *matrixarray[site];}

	private:

	CodonStateSpace* codonstatespace;
	const GTRSubMatrix* nucmatrix;
	const Array<double>* omegaarray;
	vector<MGOmegaCodonSubMatrix*> matrixarray;
};

class SiteOmegaModel	{

	double lambda;
	BranchIIDGamma* branchlength;
	
	double alpha;
	double beta;
	IIDGamma* omegaarray;

	vector<double> nucrelrate;
	vector<double> nucstat;
	GTRSubMatrix* nucmatrix;

	MGOmegaCodonSubMatrixArray* codonmatrixarray;
	BranchHomogeneousSiteHeterogeneousArray<SubMatrix>* matrixsystem;

	// suffstats

	BranchRateSuffStatArray* lengthsuffstat;
	PathSuffStatArray* pathsuffstat;
	RateSuffStatArray* omegasuffstat;
	NucPathSuffStat nucsuffstat;	

	SiteOmegaModel()	{

	}

	void MakeModel()	{

		lambda = 10.0;
		branchlength = new BranchIIDGamma(tree,1.0,lambda);

		alpha = beta = 1.0;
		omegaarray = new IIDGamma(GetNsite(),alpha,beta);

		nucrelrate ...
		nucstat ...
		nucmatrix = neew GTRSubMatrix(Nnuc,nucrelrate,nucstat,true);

		codonmatrixarray = new MGOmegaHeterogeneousCodonSubMatrixArray(codonstatespace,nucmatrix,omegaarray);

		matrixsystem = new BranchHomogeneousSiteHeterogeneousArray<SubMatrix>(codonmatrixarray);

		phyloprocess = new PhyloProcess(tree,data,branchlength,0,matrixsystem,0,codonmatrixarray);

		lengthsuffstatarray = new BranchRateSuffStatArray(tree);
		pathsuffstatarray = new PathSuffStatArray(GetNsite());
		omegasuffstatarray = new RateSuffStatArray(GetNsite());
	}

	void Move()	{

		phyloprocess->ResampleSub();

		lengthsuffstatarray->Clear();
		phyloprocess->AddLengthSuffStat(lengthsuffstatarray);
		branchlength->GibbsResample(lengthsuffstatarray);
		MoveLambda();

		pathsuffstatarray->Clear();
		phyloprocess->AddPathSuffStat(pathsuffstatarray);
		omegasuffstatarray->Clear();
		omegasuffstatarray->AddSuffStat(codonstatespace,pathsuffstatarray);
		omegaarray->GibbsResample(omegasuffstatarray);

		MoveAlpha();
		MoveBeta();

		nucsuffstat->Clear();
		nucsuffstat->Add(codonstatespace,pathsuffstatarray);
		MoveNucRates();
	}

};

class MultinomialAllocationVector {

	MultinomialAllocationVector(int insize, const vector<double>* inweight);

	void SampleAlloc();
	
	int GetSize() {return allocvector.size();}
	int GetNcomponent() {return weight.size();}

	const vector<int>& GetAllocationVector() const {return allocvector;}

	int GetAlloc(int i)	{
		return allocvector[i];
	}

	void AddSuffStat(vector<int>& occupancy);
	const vector<int>& GetOccupancyVector();
	void UpdateOccupancy();

	private:
	int size;
	int ncomp;
	vector<double>* weight;
	vector<int> allocvector;
	vector<int> occupancy;

};

template<class T> FiniteMixture : public virtual Array<T>	{

	public:

	FiniteMixture(const Array<T>* incomponents, const MultinomialAllocationVector* inalloc) : components(incompoments), alloc(inalloc)	{
	}

	~FiniteMixture() {}

	int GetSize() {return alloc->GetSize();}

	const T& GetVal(int index) const	{
		return components->GetVal(alloc->GetAlloc(i));
	}
}

class SiteOmegaMixtureModel	{

	double lambda;
	BranchIIDGamma* branchlength;
	
	int K;
	double alpha;
	double beta;
	IIDGamma* omegaarray;

	vector<double> nucrelrate;
	vector<double> nucstat;
	GTRSubMatrix* nucmatrix;

	MGOmegaCodonSubMatrixArray* codonmatrixarray;

	vector<double> weight;
	MultinomialAllocationVector alloc;

	BranchHomogeneousSiteHeterogeneousArray<SubMatrix>* matrixsystem;

	SiteOmegaMixtureModel()	{

	}

	void MakeModel()	{

		lambda = 10.0;
		branchlength = new BranchIIDGamma(tree,1.0,lambda);

		K = 5;
		alpha = beta = 1.0;
		omegaarray = new IIDGamma(K,alpha,beta);

		nucrelrate ...
		nucstat ...
		nucmatrix = neew GTRSubMatrix(Nnuc,nucrelrate,nucstat,true);

		codonmatrixarray = new MGOmegaHeterogeneousCodonSubMatrixArray(codonstatespace,nucmatrix,omegaarray);

		weight = vector<double>(K);
		alloc = new MultinomialAllocationVector(GetNsite(),weight);

		// potentially useful for auxiliary MCMC work
		// but otherwise everything should go through the mixture of matrices defined below
		// omegamixture = new FiniteMixture(omegarray,alloc);

		matrixmixture = new FiniteMixture(codonmatrixarray,alloc);

		matrixsystem = new BranchHomogeneousSiteHeterogeneousArray<SubMatrix>(matrixmixture);

		phyloprocess = new PhyloProcess(tree,data,branchlength,0,matrixsystem,0,matrixmixture);
	}

};
