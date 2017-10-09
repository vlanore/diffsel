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

#include <fstream>
#include "CodonSequenceAlignment.hpp"
#include "CodonSubMatrix.hpp"
#include "GTRSubMatrix.hpp"
#include "PhyloProcess.hpp"
#include "ProbModel.hpp"
#include "Tree.hpp"

const int Nrr = Nnuc * (Nnuc - 1) / 2;
const int Nstate = 61;

class SingleOmegaModel : public ProbModel {
    Tree* tree;
    FileSequenceAlignment* data;
    const TaxonSet* taxonset;
    CodonSequenceAlignment* codondata;

    int Nsite;
    int Ntaxa;
    int Nbranch;

    double lambda;
    double* branchlength;
    int* branchlengthcount;
    double* branchlengthbeta;

    double* nucstat;
    double* nucrelrate;
    GTRSubMatrix* nucmatrix;

    double omega;
    MGOmegaCodonSubMatrix* codonmatrix;

    SubMatrix*** phylosubmatrix;
    SubMatrix** rootsubmatrix;

    PhyloProcess* phyloprocess;

    PathSuffStat suffstat;
    double suffstatlogprob;
    double bksuffstatlogprob;

  public:
    SingleOmegaModel(std::string datafile, std::string treefile) {
        data = new FileSequenceAlignment(datafile);
        codondata = new CodonSequenceAlignment(data, true);

        Nsite = codondata->GetNsite();  // # columns
        Ntaxa = codondata->GetNtaxa();

        std::cerr << "-- Number of sites: " << Nsite << std::endl;

        taxonset = codondata->GetTaxonSet();

        // get tree from file (newick format)
        tree = new Tree(treefile);

        // check whether tree and data fits together
        tree->RegisterWith(taxonset);

        tree->SetIndices();
        Nbranch = tree->GetNbranch();

        std::cerr << "number of taxa : " << Ntaxa << '\n';
        std::cerr << "number of branches : " << Nbranch << '\n';
        std::cerr << "-- Tree and data fit together\n";

        Allocate();
        std::cerr << "-- unfold\n";
        phyloprocess->Unfold();
        std::cerr << phyloprocess->GetLogProb() << '\n';
        std::cerr << "-- mapping substitutions\n";
        phyloprocess->ResampleSub();
        std::cerr << "-- collect suffstat\n";
        CollectSuffStat();
        Trace(std::cerr);
    }

    void Allocate() {
        lambda = 10;
        branchlength = new double[Nbranch];
        for (int j = 0; j < Nbranch; j++) {
            branchlength[j] = Random::sExpo() / lambda;
        }

        branchlengthcount = new int[Nbranch];
        branchlengthbeta = new double[Nbranch];

        nucrelrate = new double[Nrr];
        double totrr = 0;
        for (int k = 0; k < Nrr; k++) {
            nucrelrate[k] = Random::sExpo();
            totrr += nucrelrate[k];
        }
        for (int k = 0; k < Nrr; k++) {
            nucrelrate[k] /= totrr;
        }

        nucstat = new double[Nnuc];
        double totstat = 0;
        for (int k = 0; k < Nnuc; k++) {
            nucstat[k] = Random::sGamma(1.0);
            totstat += nucstat[k];
        }
        for (int k = 0; k < Nnuc; k++) {
            nucstat[k] /= totstat;
        }

        nucmatrix = new GTRSubMatrix(Nnuc, nucrelrate, nucstat, true);
        omega = 1.0;
        codonmatrix = new MGOmegaCodonSubMatrix((CodonStateSpace*)codondata->GetStateSpace(),
                                                nucmatrix, omega);

        // codon matrices
        // per branch and per site
        // (array of ptrs based on condsubmatrixarray)
        phylosubmatrix = new SubMatrix**[Nbranch];
        for (int j = 0; j < Nbranch; j++) {
            phylosubmatrix[j] = new SubMatrix*[Nsite];
            for (int i = 0; i < Nsite; i++) {
                phylosubmatrix[j][i] = codonmatrix;
            }
        }

        rootsubmatrix = new SubMatrix*[Nsite];
        for (int i = 0; i < Nsite; i++) {
            rootsubmatrix[i] = codonmatrix;
        }

        // phyloprocess
        phyloprocess =
            new PhyloProcess(tree, codondata, branchlength, 0, phylosubmatrix, 0, rootsubmatrix);
    }

    void UpdateNucMatrix() {
        nucmatrix->CopyStationary(nucstat);
        nucmatrix->CorruptMatrix();
    }

    void UpdateCodonMatrix() {
        codonmatrix->SetOmega(omega);
        codonmatrix->CorruptMatrix();
    }

    void CollectSuffStat() {
        suffstat.Clear();
        RecursiveCollectSuffStat(tree->GetRoot());
    }

    void RecursiveCollectSuffStat(const Link* from) {
        if (from->isRoot()) {
            for (int i = 0; i < Nsite; i++) {
                phyloprocess->AddRootSuffStat(i, suffstat);
            }
        } else {
            for (int i = 0; i < Nsite; i++) {
                phyloprocess->AddSuffStat(i, from->Out(), suffstat);
            }
        }
        for (const Link* link = from->Next(); link != from; link = link->Next()) {
            RecursiveCollectSuffStat(link->Out());
        }
    }

    void CollectLengthSuffStat() {
        ClearLengthSuffStat();
        RecursiveCollectLengthSuffStat(tree->GetRoot());
    }

    void ClearLengthSuffStat() {
        for (int j = 0; j < Nbranch; j++) {
            branchlengthcount[j] = 0;
            branchlengthbeta[j] = 0;
        }
    }

    void RecursiveCollectLengthSuffStat(const Link* from) {
        if (!from->isRoot()) {
            for (int i = 0; i < Nsite; i++) {
                phyloprocess->AddLengthSuffStat(i, from->Out(),
                                                branchlengthcount[from->GetBranch()->GetIndex()],
                                                branchlengthbeta[from->GetBranch()->GetIndex()]);
            }
        }
        for (const Link* link = from->Next(); link != from; link = link->Next()) {
            RecursiveCollectLengthSuffStat(link->Out());
        }
    }

    void UpdateSuffStatLogProb() {
        UpdateNucMatrix();
        UpdateCodonMatrix();
        suffstatlogprob = codonmatrix->SuffStatLogProb(&suffstat);
    }

    double GetSuffStatLogProb() { return suffstatlogprob; }

    void BackupSuffStatLogProb() { bksuffstatlogprob = suffstatlogprob; }

    void RestoreSuffStatLogProb() { suffstatlogprob = bksuffstatlogprob; }

    double OmegaLogProb() { return -omega; }

    double LambdaLogProb() { return -lambda / 10; }

    double LengthLogProb() { return Nbranch * log(lambda) - lambda * GetTotalLength(); }

    double Move() {
        phyloprocess->ResampleSub();

        int nrep = 30;

        for (int rep = 0; rep < nrep; rep++) {
            CollectLengthSuffStat();
            MoveBranchLength();
            MoveLambda(1.0, 10);
            MoveLambda(0.3, 10);

            CollectSuffStat();
            UpdateSuffStatLogProb();

            MoveOmega(0.3, 3);
            MoveOmega(0.1, 3);
            MoveOmega(0.03, 3);

            MoveRR(0.1, 1, 3);
            MoveRR(0.03, 3, 3);
            MoveRR(0.01, 3, 3);

            MoveNucStat(0.1, 1, 3);
            MoveNucStat(0.01, 1, 3);
        }
        return 1.0;
    }

    double MoveRR(double tuning, int n, int nrep) {
        double nacc = 0;
        double ntot = 0;
        double bk[Nrr];
        for (int rep = 0; rep < nrep; rep++) {
            for (int l = 0; l < Nrr; l++) {
                bk[l] = nucrelrate[l];
            }
            BackupSuffStatLogProb();
            double deltalogprob = -GetSuffStatLogProb();
            double loghastings = Random::ProfileProposeMove(nucrelrate, Nrr, tuning, n);
            deltalogprob += loghastings;
            UpdateSuffStatLogProb();
            deltalogprob += GetSuffStatLogProb();
            int accepted = (log(Random::Uniform()) < deltalogprob);
            if (accepted) {
                nacc++;
            } else {
                for (int l = 0; l < Nrr; l++) {
                    nucrelrate[l] = bk[l];
                }
                RestoreSuffStatLogProb();
            }
            ntot++;
        }
        return nacc / ntot;
    }

    double MoveNucStat(double tuning, int n, int nrep) {
        double nacc = 0;
        double ntot = 0;
        double bk[Nnuc];
        for (int rep = 0; rep < nrep; rep++) {
            for (int l = 0; l < Nnuc; l++) {
                bk[l] = nucstat[l];
            }
            BackupSuffStatLogProb();
            double deltalogprob = -GetSuffStatLogProb();
            double loghastings = Random::ProfileProposeMove(nucstat, Nnuc, tuning, n);
            deltalogprob += loghastings;
            UpdateSuffStatLogProb();
            deltalogprob += GetSuffStatLogProb();
            int accepted = (log(Random::Uniform()) < deltalogprob);
            if (accepted) {
                nacc++;
            } else {
                for (int l = 0; l < Nnuc; l++) {
                    nucstat[l] = bk[l];
                }
                RestoreSuffStatLogProb();
            }
            ntot++;
        }
        return nacc / ntot;
    }

    double MoveBranchLength() {
        for (int j = 0; j < Nbranch; j++) {
            branchlength[j] =
                Random::Gamma(1.0 + branchlengthcount[j], lambda + branchlengthbeta[j]);
            if (!branchlength[j]) {
                std::cerr << "error: resampled branch length is 0\n";
                exit(1);
            }
        }
        return 1.0;
    }

    double MoveLambda(double tuning, int nrep) {
        double nacc = 0;
        double ntot = 0;
        for (int rep = 0; rep < nrep; rep++) {
            double deltalogprob = -LambdaLogProb() - LengthLogProb();
            double m = tuning * (Random::Uniform() - 0.5);
            double e = exp(m);
            lambda *= e;
            deltalogprob += LambdaLogProb() + LengthLogProb();
            deltalogprob += m;
            int accepted = (log(Random::Uniform()) < deltalogprob);
            if (accepted) {
                nacc++;
            } else {
                lambda /= e;
            }
            ntot++;
        }
        return nacc / ntot;
    }

    double MoveOmega(double tuning, int nrep) {
        double nacc = 0;
        double ntot = 0;
        for (int rep = 0; rep < nrep; rep++) {
            BackupSuffStatLogProb();
            double deltalogprob = -OmegaLogProb() - GetSuffStatLogProb();
            double m = tuning * (Random::Uniform() - 0.5);
            double e = exp(m);
            omega *= e;
            UpdateSuffStatLogProb();
            deltalogprob += OmegaLogProb() + GetSuffStatLogProb();
            deltalogprob += m;
            int accepted = (log(Random::Uniform()) < deltalogprob);
            if (accepted) {
                nacc++;
            } else {
                omega /= e;
                RestoreSuffStatLogProb();
            }
            ntot++;
        }
        return nacc / ntot;
    }

    // summary statistics

    double GetTotalLength() {
        double tot = 0;
        for (int j = 0; j < Nbranch; j++) {
            tot += branchlength[j];
        }
        return tot;
    }

    double GetLogPrior() {
        double total = 0;
        total += LambdaLogProb();
        total += LengthLogProb();
        total += OmegaLogProb();
        return total;
    }

    double GetLogLikelihood() { return phyloprocess->GetLogProb(); }

    double GetEntropy(double* profile, int dim) {
        double tot = 0;
        for (int i = 0; i < dim; i++) {
            tot -= (profile[i] < 1e-6) ? 0 : profile[i] * log(profile[i]);
        }
        return tot;
    }

    void TraceHeader(std::ostream& os) {
        os << "#logprior\tlnL\tlength\tlambda\t";
        os << "omega\t";
        os << "statent\t";
        os << "rrent\n";
    }

    void Trace(std::ostream& os) {
        os << GetLogPrior() << '\t';
        os << GetLogLikelihood() << '\t';
        os << GetTotalLength() << '\t';
        os << lambda << '\t';
        os << omega << '\t';
        os << GetEntropy(nucstat, Nnuc) << '\t';
        os << GetEntropy(nucrelrate, Nrr) << '\n';
    }

    void Monitor(std::ostream&) {}

    void FromStream(std::istream&) {}
    void ToStream(std::ostream&) {}
};
