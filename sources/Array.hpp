
#ifndef ARRAY_H
#define ARRAY_H

#include <vector>

template<class T> class ConstArray	{

	public:
	ConstArray(int insize) : size(insize) {}
	virtual ~ConstArray() {}

	int GetSize() const {return size;}

	virtual const T& GetVal(int index) const = 0;

	protected:
	int size;
};

template<class T> class Array : public ConstArray<T>	{

	public:
	Array(int insize) : ConstArray<T>(insize) {}
	~Array() {}

	virtual T& operator[](int index) = 0;
};

template<class T> class HomogeneousArray : public ConstArray<T>	{

	public:
	HomogeneousArray(int insize, const T& invalue) : ConstArray<T>(insize), value(invalue) {}
	~HomogeneousArray() {}

	const T& GetVal(int index) const override {return value;}

	private:
	const T& value;
};

template<class T> class SimpleArray : public Array<T>	{

	public:
	SimpleArray(int insize) : Array<T>(insize), array(insize) {}
	virtual ~SimpleArray() {}

	T& operator[](int index) override {return array[index];}
	const T& GetVal(int index) const override {return array[index];}

	private:
	vector<T> array;
};


#endif

