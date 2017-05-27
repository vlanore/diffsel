
#ifndef BRANCHARRAY_H
#define BRANCHARRAY_H

#include "Tree.hpp"
#include <vector>

template<class T> class ConstBranchArray	{

	public:
	ConstBranchArray(const Tree* intree) : tree(intree) {}
	virtual ~ConstBranchArray() {}

	const Tree* GetTree() const {return tree;}
	int GetNbranch() const {return tree->GetNbranch();}

	virtual const T& GetVal(int index) const = 0;

	protected:
	const Tree* tree;
};

template<class T> class BranchArray : public ConstBranchArray<T>	{

	public:
	BranchArray(const Tree* intree) : ConstBranchArray<T>(intree) {}
	~BranchArray() {}

	virtual T& operator[](int index) = 0;
};

template<class T> class HomogeneousBranchArray : public ConstBranchArray<T>	{

	public:
	HomogeneousBranchArray(const Tree* intree, const T& invalue) : ConstBranchArray<T>(intree), value(invalue) {}
	~HomogeneousBranchArray() {}

	const T& GetVal(int index) const override {return value;}

	private:
	const T& value;
};

template<class T> class SimpleBranchArray : public BranchArray<T>	{

	public:
	SimpleBranchArray(const Tree* intree) : BranchArray<T>(intree), array(intree->GetNbranch()) {}
	virtual ~SimpleBranchArray() {}

	T& operator[](int index) override {return array[index];}
	const T& GetVal(int index) const override {return array[index];}

	private:
	vector<T> array;
};

#endif
