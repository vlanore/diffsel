

template<class T> Array	{

	public:

	Array(int insize) : array(insize,0) {}
	virtual ~Array() {}

	int GetSize() {return array->size();}
	T* GetVal(int i) {return array[i];}

	protected:

	vector<T*> array;
};

template<class T> BranchArray : public Array<T> {

	public:

	BranchArray(const Tree* intree) : Array(intree->GetNbranch), tree(intree) {}

	const Tree* GetTree() {return tree;}
	const T* GetBranchVal(const Branch* branch) const {return array[branch->GetIndex()];}
	int GetNbranch() {return tree->GetNbranch();}

	T& GetBranchVal(const Branch* branch) const {return *array[branch->GetIndex()];}

	// const T& GetBranchVal(const Branch* branch) const {return *array[branch->GetIndex()];}

	virtual const T& GetBranchVal(const Branch* branch) = 0;

	private:

	const Tree* tree;
};

class BranchIIDExponential : public BranchArray<double>	{

	public:

	BranchIIDExponential(const Tree* intree, double inlambda) : BranchArray<double>(intree), lambda(inlambda)	{
		branchlength = vector<double>(tree->GetNbranch()];
		suffstatcount = vector<int>[tree->GetNbranch()];
		suffstatbeta = vector<double>[tree->GetNbranch()];
		Sample();
		for (int i=0; i<GetNbranch(); i++)	{
			array[i] = &branchlength[i];
		}
	}

	~BranchIIDExponential()	{}

	void Sample()	{
		for (int i=0; i<GetNbranch(); i++)	{
			branchlength[i] = Random::sExpo() / lambda;
		}
	}

	double GetLambda()	{
		return lambda;
	}

	void SetLambda(double inlambda)	{
		lambda = inlambda;
	}

	// log probability
	// sum over branches
	double GetLogProb();
	// for a given branch
	double GetLogProb(int branchindex);

	// suff stats for this prior probability
	void AddSuffStat(double& sum, double& sumoflogs, int& nbranch)

	// suff stat moves

	// based on suffstat
	// void GibbsSample();

	void SetPhyloProcess(PhyloProcess* in phyloprocess);

	void CollectSuffStat();
	void ClearSuffStat();
	void RecursiveCollectSuffStat();
	void ResampleBranchLengths();

	private:
	vector<double> branchlength;
	double lambda;
	vector<int> suffstatcount;
	vector<double> suffstatbeta;

	PhyloProcess* inphyloprocess;
};

template<class T> IIDArray	{

	virtual void CreateComponent(int k) = 0;
	virtual void DeleteComponent(int k) = 0;
	virtual void SampleComponent(int k) = 0;
	
	const T& GetComponent(int k) const;
	T& GetComponent(int k);
};

class Dirichlet	{

	void Resample(vector<int>& suffstat);
};


template<class T> class FiniteMixture : public virtual Array<T>	{

	FiniteMixture(int insize, const Array<T>* incomponents, vector<double>* inweights) : 
				Array<T>(size), components(incomponents), weights(inweights), alloc(insize), occupancy(incomponents->GetSize()) {
		
		SampleAlloc();
		UpdateAlloc();
	}

	// FiniteMixture(const vector<T>& components, vector<double>& weights);
	virtual ~FiniteMixture();

	T* GetVal(int i)	{
		return component[alloc[i]];
	}

	void UpdateAlloc()	{
		for (int i=0; i<this->GetSize(); i++)	{
			array[i] = &component[alloc[i]];
		}
	}

	void UpdateOccupancy()	{
		for (int k=0; k<occupancy.size(); k++)	{
			occupancy[k] = 0;
		}
		for (int i=0; i<GetSize(); i++)	{
			occupancy[alloc[i]]++;
		}
	}

	void AddWeightSuffStat(vector<int>& occupancysuffstat)	{
		UpdateOccupancy();
		for (int k=0; k<occupancy.size(); k++)	{
			occupancysuffstat[k] += occupancy[k];
		}
	}

	void RemoveItem(int i)	{
		occupancy[alloc[i]]--;
		alloc[i] = -1;
	}
	void AddItem(int i, int k)	{
		alloc[i] = k;
		occupancy[alloc[i]]++;
	}

	virtual double AllocationLogProb(int i, int k) = 0;
	// uses AllocationLogProb
	void GibbsResample();

	private:

	const Array<T>* components;
	vector<double>* weights;
	vector<int> alloc;
	vector<int> occupancy;
};

