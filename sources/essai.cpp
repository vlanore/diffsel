
#include <vector>
#include <iostream>

class real	{

	public:
	real(double inval) : val(inval) {}

	// this does not work
	double& getVal() {return val;}
	const double& getVal() const {return val;}

	private:
	double val;
};

int main()	{

	real a(2.0);
	std::cerr << a.getVal() << '\n';
}

