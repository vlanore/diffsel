#include <cmath>
#include <fstream>
#include "Chain.hpp"
#include "DiffSelModel.hpp"
using namespace std;


class DiffSelChain : public Chain {
  private:
    // Chain parameters
    string modeltype, datafile, treefile;
    int codonmodel, category, level, fixglob, fixvar;

  public:
    DiffSelModel* GetModel() { return static_cast<DiffSelModel*>(model); }

    string GetModelType() override { return modeltype; }

    DiffSelChain(string indata, string intree, int incategory, int inlevel, int inevery,
                 int inuntil, int infixglob, int infixvar, int incodonmodel, string inname,
                 int force)
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

    DiffSelChain(string filename) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        model = new DiffSelModel(datafile, treefile, category, level, fixglob, fixvar, codonmodel,
                                 true);
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
            model = new DiffSelModel(datafile, treefile, category, level, fixglob, fixvar,
                                     codonmodel, false);
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
        DiffSelChain* chain = new DiffSelChain(name);
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
            cerr << "error in command\n";
            exit(1);
        }
        DiffSelChain* chain = new DiffSelChain(datafile, treefile, ncond, nlevel, every, until,
                                               fixglob, fixvar, codonmodel, name, true);
        cerr << "start\n";
        chain->Start();
        cerr << "chain stopped\n";
    }
}
