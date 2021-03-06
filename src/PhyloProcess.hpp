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
            std::cerr << "error in PhyloProcess::GetSubMatrix: called on root\n";
            exit(1);
        }
        SubMatrix* tmp = submatrix[branch->GetIndex()][site];
        if (!tmp) {
            std::cerr << "error in PhyloProcess::GetSubMatrix: no sub matrix\n";
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

    bool isMissing(const Node* /*node*/, int /*site*/) {
        return false;
        // return missingmap[node][site];
    }

    bool isMissing(const Link* /*link*/, int /*site*/) {
        return false;
        // return (missingmap[link->GetNode()][site] || missingmap[link->Out()->GetNode()][site]);
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
