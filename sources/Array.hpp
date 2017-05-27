
#ifndef ARRAY_H
#define ARRAY_H

#include "Tree.hpp"
#include <vector>

template<class T> class Array	{

	public:
	Array(int insize) : size(insize) {}
	virtual ~Array() {}

	int GetSize() const {return size;}
	virtual T& GetVal(int index) = 0;
	virtual const T& GetVal(int index) const = 0;

	protected:
	int size;
};

template<class T> class HomogeneousArray : public virtual Array<T>	{

	public:
	HomogeneousArray(int insize, T* invalue) : Array<T>(insize), value(invalue) {}
	~HomogeneousArray() {}

	T& GetVal(int index) override {return *value;}
	const T& GetVal(int index) const override {return *value;}

	private:
	T* value;
};

template<class T> class BranchArray	{

	public:
	BranchArray(const Tree* intree) : tree(intree) {}
	virtual ~BranchArray() {}

	const Tree* GetTree() const {return tree;}
	int GetNbranch() const {return tree->GetNbranch();}

	virtual T& GetBranchVal(int branch) = 0;
	virtual const T& GetBranchVal(int branch) const = 0;

	protected:
	const Tree* tree;
};

template<class T> class BranchHomogeneousArray : public virtual BranchArray<T> {

	public:

	BranchHomogeneousArray(const Tree* intree, T* invalue) : BranchArray<T>(intree), value(invalue)	{}
	~BranchHomogeneousArray() {}

	T& GetBranchVal(int branch) override {return *value;}
	const T& GetBranchVal(int branch) const override {return *value;}

	private:
	T* value;

	
};

template<class T> class SimpleArray : public virtual Array<T>	{

	public:
	SimpleArray(int insize) : Array<T>(insize), array(insize) {}
	virtual ~SimpleArray() {}

	T& GetVal(int index) override {return array[index];}
	const T& GetVal(int index) const override {return array[index];}

	private:
	vector<T> array;
};

template<class T> class SimpleBranchArray : public virtual SimpleArray<T>, public virtual BranchArray<T> {

	public:

	SimpleBranchArray(const Tree* intree) : Array<T>(intree->GetNbranch()), SimpleArray<T>(intree->GetNbranch()), BranchArray<T>(intree) {}
	~SimpleBranchArray() {}

	T& GetBranchVal(int index) override {return this->GetVal(index);}
	const T& GetBranchVal(int index) const override {return this->GetVal(index);}

};

/*
template<class T> class BranchPtrArray : public virtual BranchArray<T*>	{

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
*/

template<class T> class BranchSiteArray	{

	public:
	BranchSiteArray(const Tree* intree, int insize) : tree(intree), size(insize) {}
	virtual ~BranchSiteArray() {}

	const Tree* GetTree() const {return tree;}
	int GetNbranch() const {return tree->GetNbranch();}
	int GetSize() const {return size;}

	virtual T& GetVal(int branch, int site) = 0;
	virtual const T& GetVal(int branch, int site) const = 0;

	protected:
	const Tree* tree;
	int size;
};


template<class T> class HomogeneousBranchSiteArray : public virtual BranchSiteArray<T> {

	public:
	HomogeneousBranchSiteArray(const Tree* intree, int insize, T* invalue) : BranchSiteArray<T>(intree,insize), value(invalue) {}
	~HomogeneousBranchSiteArray() {}

	T& GetVal(int branch, int site) override {return *value;}
	const T& GetVal(int branch, int site) const override {return *value;}

	private:
	T* value;
};

template<class T> class BranchHomogeneousSiteHeterogeneousArray : public virtual BranchSiteArray<T>	{

	public:
	BranchHomogeneousSiteHeterogeneousArray(const Tree* intree, Array<T>* inarray) : BranchSiteArray<T>(intree,inarray->GetSize()), array(inarray) {}
	~BranchHomogeneousSiteHeterogeneousArray() {}

	T& GetVal(int branch, int site) override {return array->GetVal(site);}
	const T& GetVal(int branch, int site) const override {return array->GetVal(site);}

	private:
	Array<T>* array;
};

template<class T> class BranchHeterogeneousSiteHomogeneousArray : public virtual BranchSiteArray<T> {

	public:
	BranchHeterogeneousSiteHomogeneousArray(BranchArray<T>* inbrancharray, int insize) : BranchSiteArray<T>(inbrancharray->GetTree(),insize), brancharray(inbrancharray) {}
	~BranchHeterogeneousSiteHomogeneousArray() {}

	T& GetVal(int branch, int site) override {return brancharray->GetBranchVal(branch);}
	const T& GetVal(int branch, int site) const override {return brancharray->GetBranchVal(branch);}

	private:
	BranchArray<T>* brancharray;
};

#endif

