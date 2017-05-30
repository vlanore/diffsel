#include <cmath>
#include <fstream>
#include "SiteOmegaModel.hpp"
using namespace std;

int main(int argc, char* argv[])	{

	string datafile = argv[1];
	string treefile = argv[2];
	string name = argv[3];

	SiteOmegaModel* model = new SiteOmegaModel(datafile,treefile);
	ofstream os((name + ".trace").c_str());
	model->TraceHeader(os);
	os.flush();
	while(1)	{
		model->Move();
		model->Trace(os);
		os.flush();
	}
}



