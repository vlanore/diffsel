#include <cmath>
#include <fstream>
#include "SingleOmegaModel.hpp"
using namespace std;

int main(int, char* argv[]) {
    string datafile = argv[1];
    string treefile = argv[2];
    string name = argv[3];

    SingleOmegaModel* model = new SingleOmegaModel(datafile, treefile);
    ofstream os((name + ".trace").c_str());
    model->TraceHeader(os);
    os.flush();
    while (1) {
        model->Move();
        model->Trace(os);
        os.flush();
    }
}
