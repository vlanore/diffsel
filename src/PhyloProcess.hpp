#ifndef PHYLOPROCESS_H
#define PHYLOPROCESS_H

#include <map>
#include "BranchSitePath.hpp"
#include "Chrono.hpp"
#include "SequenceAlignment.hpp"
#include "SubMatrix.hpp"
#include "Tree.hpp"

// PhyloProcess is a dispatcher:
// its responsibility is to create a random branch/site path
// for each branch/site pair

class PhyloProcess {
  public:
    PhyloProcess(Tree* intree, SequenceAlignment* indata, double* inbranchlength,
                 double* insiterate, SubMatrix*** insubmatrix, double** inrootfreq,
                 SubMatrix** inrootsubmatrix);

    ~PhyloProcess();

    // accessors

    double GetBranchLength(const Branch* branch) { return branchlength[branch->GetIndex()]; }

    double GetSiteRate(int site) {
        if (siterate == nullptr) {
            return 1.0;
        }
        return siterate[site];
    }

    SubMatrix* GetSubMatrix(const Branch* branch, int site) {
        if (!branch) {
            cerr << "error in PhyloProcess::GetSubMatrix: called on root\n";
            exit(1);
        }
        SubMatrix* tmp = submatrix[branch->GetIndex()][site];
        if (!tmp) {
            cerr << "error in PhyloProcess::GetSubMatrix: no sub matrix\n";
            exit(1);
        }
        return tmp;
        // return submatrix[branch->GetIndex()][site];
    }

    const double* GetRootFreq(int site) {
        if (rootfreq) {
            return rootfreq[site];
        }
        if (!rootsubmatrix) {
            std::cerr << "error in PhyloProcess::GetRootFreq\n";
            exit(1);
        }
        return rootsubmatrix[site]->OldGetStationary();
    }


    double SiteLogLikelihood(int site);
    double FastSiteLogLikelihood(int site);

    const StateSpace* GetStateSpace() { return data->GetStateSpace(); }
    const TaxonSet* GetTaxonSet() { return data->GetTaxonSet(); }

    int GetNsite() { return data->GetNsite(); }
    int GetNtaxa() { return data->GetNtaxa(); }

    int GetNstate() { return Nstate; }

    SequenceAlignment* GetData() { return data; }
    int GetData(int taxon, int site) { return data->GetState(taxon, site); }

    Tree* GetTree() { return tree; }
    Link* GetRoot() { return GetTree()->GetRoot(); }

    int GetMaxTrial() { return maxtrial; }
    void SetMaxTrial(int i) { maxtrial = i; }

    void SetData(SequenceAlignment* indata);
    void ClampData() { clampdata = true; }
    void UnclampData() { clampdata = false; }

    int& GetState(const Node* node, int site) { return statemap[node][site]; }

    void GetLeafData(SequenceAlignment* data);
    void RecursiveGetLeafData(const Link* from, SequenceAlignment* data);
    bool isDataCompatible(int taxon, int site, int state) {
        return GetStateSpace()->isCompatible(GetData(taxon, site), state);
    }

    // probability, pruning, sampling
    double GetLogProb();                                // likelihood Felsenstein 1981
    double GetFastLogProb();                            // likelihood Felsenstein 1981
    double GetPathLogProb();                            // probability of the entire mapping
    double GetPathLogProb(const Link* from, int site);  // probability of the entire mapping

    double Move(double fraction);
    void DrawSites(double fraction);  // draw a fraction of sites which will be resampled
    void ResampleSub();               // clamped Nielsen
    void ResampleSub(int site);

    void AddSuffStat(int site, const Link* link, PathSuffStat& suffstat);
    void AddRootSuffStat(int site, PathSuffStat& suffstat);
    void AddLengthSuffStat(int site, const Link* link, int& count, double& beta);

    void PostPredSample(bool rootprior = false);  // unclamped Nielsen
    void PostPredSample(int site, bool rootprior = false);
    // rootprior == true : root state drawn from stationary probability of the
    // process
    // rootprior == false: root state drawn from posterior distribution

    // various accessors

    bool isMissing(const Node* node, int site) {
        return false;
        return missingmap[node][site];
    }

    bool isMissing(const Link* link, int site) {
        return false;
        return (missingmap[link->GetNode()][site] || missingmap[link->Out()->GetNode()][site]);
    }

    void CreateMissingMap();
    void RecursiveCreateMissingMap(const Link* from);
    bool FillMissingMap(const Link* from, int i);

    double* GetCondLikelihood(const Link* from) { return condlmap[from]; }

    double GetPruningTime() { return pruningchrono.GetTime(); }
    double GetResampleTime() { return resamplechrono.GetTime(); }


    void Unfold();
    void Cleanup();

  protected:
    void RecursiveCreate(const Link* from);
    void RecursiveDelete(const Link* from);

    void RecursiveCreateTBL(const Link* from);
    void RecursiveDeleteTBL(const Link* from);

    /*
    void BackwardPropagate(const double *down, double *up) {
        GetSubMatrix()->BackwardPropagate(down, up, GetTime() * GetRate());
    }
    void ForwardPropagate(const double *up, double *down) {
        GetSubMatrix()->ForwardPropagate(up, down, GetTime() * GetRate());
    }
    */

    void Pruning(const Link* from, int site);
    void ResampleSub(const Link* from, int site);
    void ResampleState();
    void ResampleState(int site);
    void PruningAncestral(const Link* from, int site);
    void PriorSample(const Link* from, int site, bool rootprior);
    void PriorSample();
    void RootPosteriorDraw(int site);

    // borrowed from phylobayes
    // where should that be?
    BranchSitePath* SamplePath(int stateup, int statedown, double time, double rate,
                               SubMatrix* matrix);
    BranchSitePath* SampleRootPath(int rootstate);
    BranchSitePath* ResampleAcceptReject(int maxtrial, int stateup, int statedown, double rate,
                                         double totaltime, SubMatrix* matrix);
    BranchSitePath* ResampleUniformized(int stateup, int statedown, double rate, double totaltime,
                                        SubMatrix* matrix);

    Tree* tree;
    SequenceAlignment* data;
    double* branchlength;
    double* siterate;
    SubMatrix*** submatrix;
    double** rootfreq;
    SubMatrix** rootsubmatrix;

    int* sitearray;
    double* sitelnL;
    std::map<const Link*, double*> condlmap;

    int Nstate;

    bool clampdata;

    BranchSitePath* GetPath(const Node* node, int site) {
        if (pathmap[node][site] == nullptr) {
            std::cerr << "error in phyloprocess::getpath: null path\n";
            exit(1);
        }
        return pathmap[node][site];
    }

  private:
    std::map<const Node*, BranchSitePath**> pathmap;
    std::map<const Node*, int*> statemap;
    std::map<const Node*, bool*> missingmap;
    std::map<const Node*, int> totmissingmap;

    int maxtrial;
    static const int unknown = -1;

    static const int DEFAULTMAXTRIAL = 100;

    Chrono pruningchrono;
    Chrono resamplechrono;
};

#endif  // PHYLOPROCESS_H
