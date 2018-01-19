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

#define DOCTEST_CONFIG_IMPLEMENT

#include <fstream>
#include "Chain.hpp"
#include "DiffSelSparseModel.hpp"

using namespace std;

class DiffSelSparseChain : public Chain {
  private:
    // Chain parameters
    string modeltype, datafile, treefile;
    int codonmodel, category, level, fixglob, fixvar;

  public:
    DiffSelSparseModel::AcceptStats stats;

    DiffSelSparseModel* GetModel() { return static_cast<DiffSelSparseModel*>(model.get()); }

    string GetModelType() override { return modeltype; }

    DiffSelSparseChain(const string& indata, const string& intree, int incategory, int inlevel,
                       int inevery, int inuntil, int infixglob, int infixvar, int incodonmodel,
                       string inname, int force)
        : modeltype("DIFFSEL"),
          datafile(indata),
          treefile(intree),
          codonmodel(incodonmodel),
          category(incategory),
          level(inlevel),
          fixglob(infixglob),
          fixvar(infixvar) {
        SetEvery(inevery);
        SetUntil(inuntil);
        SetName(inname);
        New(force);
    }

    explicit DiffSelSparseChain(string filename) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        model = std::unique_ptr<DiffSelSparseModel>(new DiffSelSparseModel(
            datafile, treefile, category, level, fixglob, fixvar, codonmodel, true, &stats));
        cerr << "-- Reset" << endl;
        Reset(force);
        cerr << "-- New ok\n";
        model->Trace(cerr);
    }

    void Open() override {
        ifstream is((name + ".param").c_str());
        if (!is) {
            cerr << "-- Error : cannot find file : " << name << ".param\n";
            exit(1);
        }
        is >> modeltype;
        is >> datafile >> treefile >> category >> level;
        is >> fixglob >> fixvar >> codonmodel;
        int tmp;
        is >> tmp;
        if (tmp) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> every >> until >> size;

        if (modeltype == "DIFFSEL") {
            model = std::unique_ptr<DiffSelSparseModel>(new DiffSelSparseModel(
                datafile, treefile, category, level, fixglob, fixvar, codonmodel, true));
        } else {
            cerr << "-- Error when opening file " << name
                 << " : does not recognise model type : " << modeltype << '\n';
            exit(1);
        }
        model->FromStream(is);

        model->Update();
        cerr << size << "-- Points saved\n";
        model->Trace(cerr);
    }

    void Save() override {
        ofstream param_os((name + ".param").c_str());
        param_os << GetModelType() << '\n';
        param_os << datafile << '\t' << treefile << '\t' << category << '\t' << level << '\n';
        param_os << fixglob << '\t' << fixvar << '\t' << codonmodel << '\n';
        param_os << 0 << '\n';
        param_os << every << '\t' << until << '\t' << size << '\n';

        model->ToStream(param_os);
    }


    void MakeFiles(int force = 0) override {
        Chain::MakeFiles(force);
        std::ofstream os((name + ".fitness").c_str());
    }

    void SavePoint() override {
        Chain::SavePoint();
        std::ofstream os((name + ".fitness").c_str(), ios_base::app);
        GetModel()->PrintFitnessVector(os);
    }

    void Move() override {
        for (int i = 0; i < every; i++) {
            model->Move();
        }

        SavePoint();
        Save();
        Monitor();
    }
};


int main(int argc, char* argv[]) {
    cerr << "-- Parsing command line arguments\n";

    // this is an already existing chain on the disk; reopen and restart
    if (argc == 2 && argv[1][0] != '-') {
        string name = argv[1];
        cerr << "-- Trying to reopen existing chain named " << name << " on disk\n";
        DiffSelSparseChain* chain = new DiffSelSparseChain(name);
        cerr << "start\n";
        chain->Start();
        cerr << "chain stopped\n";
    }

    // this is a new chain
    else {
        string datafile = "";
        string treefile = "";
        int ncond = 2;
        int nlevel = 2;
        int fixglob = 1;
        int fixvar = 1;
        int codonmodel = 1;

        string name = "";
        int every = 1;
        int until = -1;

        try {
            if (argc == 1) {
                throw(0);
            }

            int i = 1;
            while (i < argc) {
                string s = argv[i];

                if (s == "-d") {
                    i++;
                    datafile = argv[i];
                } else if ((s == "-t") || (s == "-T")) {
                    i++;
                    treefile = argv[i];
                } else if (s == "-ncond") {
                    i++;
                    ncond = atoi(argv[i]);
                } else if (s == "-nlevel") {
                    i++;
                    nlevel = atoi(argv[2]);
                    if ((nlevel != 1) && (nlevel != 2)) {
                        cerr << "error: nlevel should be 1 or 2\n";
                        exit(1);
                    }
                } else if (s == "-fixvar") {
                    fixvar = 1;
                } else if (s == "-freevar") {
                    fixvar = 0;
                } else if (s == "-ms") {
                    codonmodel = 1;
                } else if (s == "-sr") {
                    codonmodel = 0;
                } else if ((s == "-x") || (s == "-extract")) {
                    i++;
                    if (i == argc) throw(0);
                    every = atoi(argv[i]);
                    i++;
                    if (i == argc) throw(0);
                    until = atoi(argv[i]);
                } else {
                    if (i != (argc - 1)) {
                        throw(0);
                    }
                    name = argv[i];
                }
                i++;
            }
        } catch (...) {
            cerr << "error in command\n"
                 << "usage is:\n"
                 << "\tflatdiffsel_bin -d datafile -t treefile [-ncond i] [-nlevel 1|2] "
                    "[-fixvar|-freevar] [-ms|-sr] [-x every until] name\n";
            exit(1);
        }
        DiffSelSparseChain* chain =
            new DiffSelSparseChain(datafile, treefile, ncond, nlevel, every, until, fixglob, fixvar,
                                   codonmodel, name, true);
        cerr << "start\n";
        chain->Start();
        cerr << "chain stopped\n";

        for (auto p : chain->stats.d) {
            cout << setw(35) << p.first << " -> "
                 << 100. * accumulate(p.second.begin(), p.second.end(), 0.) / p.second.size()
                 << "%\n";
        }
        delete chain;
    }
}
