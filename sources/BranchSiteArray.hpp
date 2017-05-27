
#ifndef BRANCHSITEARRAY_H
#define BRANCHSITEARRAY_H

#include "Array.hpp"
#include "BranchArray.hpp"

template<class T> class ConstBranchSiteArray	{

	public:
	ConstBranchSiteArray(const Tree* intree, int insize) : tree(intree), size(insize) {}
	virtual ~ConstBranchSiteArray() {}

	const Tree* GetTree() const {return tree;}
	int GetNbranch() const {return tree->GetNbranch();}
	int GetSize() const {return size;}

	virtual const T& GetVal(int branch, int site) const = 0;

	protected:
	const Tree* tree;
	int size;
};

template<class T> class BranchSiteArray : public ConstBranchSiteArray<T> {

	public:
	BranchSiteArray(const Tree* intree, int insize) : ConstBranchSiteArray<T>(intree,insize) {}
	virtual ~BranchSiteArray() {}

	virtual T& operator()(int branch, int site) = 0;
};

template<class T> class HomogeneousBranchSiteArray : public virtual ConstBranchSiteArray<T> {

	public:
	HomogeneousBranchSiteArray(const Tree* intree, int insize, const T& invalue) : ConstBranchSiteArray<T>(intree,insize), value(invalue) {}
	~HomogeneousBranchSiteArray() {}

	const T& GetVal(int branch, int site) const override {return value;}

	private:
	const T& value;
};

template<class T> class BranchHomogeneousSiteHeterogeneousArray : public virtual ConstBranchSiteArray<T>	{

	public:
	BranchHomogeneousSiteHeterogeneousArray(const Tree* intree, const ConstArray<T>* inarray) : ConstBranchSiteArray<T>(intree,inarray->GetSize()), array(inarray) {}
	~BranchHomogeneousSiteHeterogeneousArray() {}

	const T& GetVal(int branch, int site) const override {return array->GetVal(site);}

	private:
	const ConstArray<T>* array;
};

template<class T> class BranchHeterogeneousSiteHomogeneousArray : public virtual ConstBranchSiteArray<T> {

	public:
	BranchHeterogeneousSiteHomogeneousArray(const ConstBranchArray<T>* inbrancharray, int insize) : ConstBranchSiteArray<T>(inbrancharray->GetTree(),insize), brancharray(inbrancharray) {}
	~BranchHeterogeneousSiteHomogeneousArray() {}

	const T& GetVal(int branch, int site) const override {return brancharray->GetVal(branch);}

	private:
	const ConstBranchArray<T>* brancharray;
};

#endif
