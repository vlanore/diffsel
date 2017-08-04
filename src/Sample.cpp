/*Copyright or © or Copr. Centre National de la Recherche Scientifique (CNRS) (2017-06-14).
Contributors:
* Nicolas LARTILLOT - nicolas.lartillot@univ-lyon1.fr

This software is a computer program whose purpose is to detect convergent evolution using Bayesian
phylogenetic codon models.

This software is governed by the CeCILL-C license under French law and abiding by the rules of
distribution of free software. You can use, modify and/ or redistribute the software under the terms
of the CeCILL-C license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and rights to copy, modify and redistribute
granted by the license, users are provided only with a limited warranty and the software's author,
the holder of the economic rights, and the successive licensors have only limited liability.

In this respect, the user's attention is drawn to the risks associated with loading, using,
modifying and/or developing or reproducing the software by the user in light of its specific status
of free software, that may mean that it is complicated to manipulate, and that also therefore means
that it is reserved for developers and experienced professionals having in-depth computer knowledge.
Users are therefore encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or data to be ensured and,
more generally, to use and operate it in the same conditions as regards security.

The fact that you are presently reading this means that you have had knowledge of the CeCILL-C
license and that you accept its terms.*/

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

    chain_is = ifstream{(name + ".chain").c_str()};  // move construction (c++11)

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
