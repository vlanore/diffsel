#define DOCTEST_CONFIG_IMPLEMENT
#include <cmath>
#include <fstream>
#include "SingleOmegaModel.hpp"
using namespace std;

int main(int argc, char* argv[]) {
    cerr << "Argc " << argc << endl;
    string datafile = argv[1];
    string treefile = argv[2];
    string name = argv[3];
    int it{-1};
    if (argc > 4) {
        cerr << "Dun dun dun\n";
        it = atoi(argv[4]);
    }

    SingleOmegaModel* model = new SingleOmegaModel(datafile, treefile);
    ofstream os((name + ".trace").c_str());
    model->TraceHeader(os);
    os.flush();
    auto i = 0;
    while (i < it or it == -1) {
        model->Move();
        model->Trace(os);
        os.flush();
        i++;
    }
}
