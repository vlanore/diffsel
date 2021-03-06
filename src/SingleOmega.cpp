/*Copyright or © or Copr. Centre National de la Recherche Scientifique (CNRS) (2017-06-14).
Contributors:
* Nicolas LARTILLOT - nicolas.lartillot@univ-lyon1.fr
* Vincent LANORE - vincent.lanore@univ-lyon1.fr

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

#define DOCTEST_CONFIG_IMPLEMENT
#include <cmath>
#include <fstream>
#include "SingleOmegaModel.hpp"
using namespace std;

int main(int argc, char* argv[]) {
    string datafile = argv[1];
    string treefile = argv[2];
    string name = argv[3];
    int it{-1};
    if (argc > 4) {
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
