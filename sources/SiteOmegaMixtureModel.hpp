
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

