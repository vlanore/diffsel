#include <cmath>
#include <fstream>
#include "MultiGeneSingleOmegaModel.hpp"
using namespace std;

int main(int argc, char* argv[])	{

	string genelistfile = argv[1];
	string treefile = argv[2];
	string name = argv[3];

	MultiGeneSingleOmegaModel* model = new MultiGeneSingleOmegaModel(genelistfile,treefile);
	ofstream os((name + ".trace").c_str());
	model->TraceHeader(os);
	os.flush();
	while(1)	{
		model->Move();
		model->Trace(os);
		os.flush();
	}
}


