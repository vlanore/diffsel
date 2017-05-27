

template<class T> class BranchArray	{

	public:

	BranchArray(const Tree* intree) : tree(intree), Nbranch(tree->GetNbranch()), array(tree->GetNbranch())	{

	}

	int GetNbranch() {return Nbranch;}

	T& GetBranchVal(int branchindex)	{return array[index];}

	private:

	const Tree* tree;
	int Nbranch;
	vector<T*> array;
};


BranchIIDExpo : public BranchArray<double> {


	
	double* branchlength;

};

