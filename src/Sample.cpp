#include "Sample.hpp"
#include <utility>
#include "ProbModel.hpp"

Sample::Sample(string filename, int in_burnin, int in_every, int in_until) {
    burnin = in_burnin;
    every = in_every;
    until = in_until;
    name = filename;
}

void Sample::OpenChainFile() {
    if (until == -1) {
        until = chainsize;
    }
    if (until > chainsize) {
        cerr << "number of points saved is less than " << until << '\n';
        until = chainsize;
    }
    if (burnin == -1) {
        burnin = chainsize / 10;
    }
    size = (until - burnin) / every;
    if (size <= 0) {
        cerr << "error : chain not long enough\n";
        exit(1);
    }
    currentpoint = 0;

    ifstream newstream{(name + ".chain").c_str()};
    chain_is.swap(newstream);

    if (!chain_is) {
        cerr << "error: cannot find file " << name << ".chain\n";
        exit(1);
    }
    for (int i = 0; i < burnin; i++) {
        model->FromStream(chain_is);
    }
}

void Sample::GetNextPoint() {
    if (currentpoint == size) {
        cerr << "error in Sample::GetNextPoint: going past last points\n";
        exit(1);
    }
    if (currentpoint) {
        for (int i = 0; i < every - 1; i++) {
            model->FromStream(chain_is);
        }
    }
    model->FromStream(chain_is);
    currentpoint++;
}
