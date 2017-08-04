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

#include "Chain.hpp"
#include <fstream>
#include <iostream>
#include "Chrono.hpp"
#include "ProbModel.hpp"
using namespace std;

Chain::Chain() {
    every = 1;
    until = -1;
    size = 0;
    model = nullptr;
    name = "";
}

void Chain::MakeFiles(int force) {
    if (ifstream((name + ".param").c_str()) && (force == 0)) {
        cerr << "already existing chain, cannot override (unless in forcing mode)\n";
        exit(1);
    }
    // ofstream param_os((name + ".param").c_str());
    // ofstream chain_os((name + ".chain").c_str());
    // ofstream mon_os((name + ".monitor").c_str());
    ofstream trace_os((name + ".trace").c_str());
    model->TraceHeader(trace_os);
}

void Chain::Monitor() {
    ofstream trace_os((name + ".trace").c_str(), ios_base::app);
    model->Trace(trace_os);
    ofstream mon_os((name + ".monitor").c_str());
    // ofstream mon_det_os((name + ".details").c_str());
    model->Monitor(mon_os);
}

void Chain::SavePoint() {
    ofstream chain_os((name + ".chain").c_str(), ios_base::app);
    model->ToStream(chain_os);
    size++;
}

void Chain::Reset(int force) {
    size = 0;
    MakeFiles(force);
    Save();
}

void Chain::Move() {
    for (int i = 0; i < every; i++) {
        model->Move();
    }
    SavePoint();
    Save();
    Monitor();
}

void Chain::Start() {
    ofstream run_os((name + ".run").c_str());
    run_os << 1 << '\n';
    run_os.close();
    Run();
}

int Chain::GetRunningStatus() {
    ifstream run_is((name + ".run").c_str());
    int run;
    run_is >> run;
    return run;
}

void Chain::Run() {
#define DEBUG 2
#if DEBUG > 0
    int i = 0;
    MeasureTime timer;
#endif
    while ((GetRunningStatus() != 0) && ((until == -1) || (size <= until))) {
        Chrono chrono;
        chrono.Reset();
        chrono.Start();
        Move();
        chrono.Stop();
#if DEBUG > 0
        timer << "Iteration " << i * every << ". ";
        timer.print<0>();
        i++;
#endif
        /*
        ofstream check_os((name + ".time").c_str());
        check_os << chrono.GetTime() / 1000 << '\n';
        */
    }
    ofstream run_os((name + ".run").c_str());
    run_os << 0 << '\n';
}
