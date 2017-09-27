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


#include "CodonSequenceAlignment.hpp"
#include "GTRSubMatrix.hpp"
#include "MSCodonSubMatrix.hpp"
#include "PhyloProcess.hpp"
#include "ProbModel.hpp"
#include "Tree.hpp"

const int Nrr = Nnuc * (Nnuc - 1) / 2;
const int Nstate = 61;

using AAProfile = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using BMatrix = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>;

void InitUniformDirichlet(Eigen::VectorXd& v) {
    double tot = 0.;
    for (int i = 0; i < v.size(); i++) {
        v[i] = Random::sExpo();
        tot += v[i];
    }
    for (int i = 0; i < v.size(); i++) {
        v[i] /= tot;
    }
}

class DiffSelSparseModel : public ProbModel {
    // -----
    // model selectors
    // -----

    int fixglob;
    int fixvar;
    int codonmodel;

    // -----
    // external parameters
    // -----

    Tree* tree;
    FileSequenceAlignment* data;
    const TaxonSet* taxonset;
    CodonSequenceAlignment* codondata;

    // number of sites
    int Nsite;
    int Ntaxa;
    int Nbranch;

    // number of diff sel categories
    int Ncond;

    // number of levels of the model
    int Nlevel;

    // -----
    //  model structure
    // -----

    // iid exponential branch lengths of rate lambda
    double lambda;
    double* branchlength;

    // suff stats associated with branch lengths
    int* branchlengthcount;
    double* branchlengthbeta;

    // which branch is under which condition
    int* branchalloc;

    // nucleotide exchange rates, eq frequencies and GTR substitution matrix
    double* nucrelrate;
    double* nucstat;
    GTRSubMatrix* nucmatrix;


    double fitness_shape;
    AAProfile fitness_inv_rates;
    // differential selection factors across conditions k=1..Ncond and across sites
    // Ncond * Nsite * Naa
    std::vector<Eigen::MatrixXd> fitness;

    double prob_conv_m { 0.1 };
    double prob_conv_v { 0.5 };
    Eigen::VectorXd prob_conv;  // indexed by condition
    vector<BMatrix> ind_conv;   // indexed by condition * sites * aa;


    // codon substitution matrices
    // across conditions and across sites
    // Ncond * Nsite
    CodonSubMatrix*** condsubmatrixarray;

    // branch- and site-pointers over substitution matrices (for phyloprocess)
    SubMatrix*** phylosubmatrix;

    // site-pointers over substitution matrices of condition 0
    // their equilibrium frequency vectors are used for the root
    // (for phyloprocess)
    SubMatrix** rootsubmatrix;

    // suff stats, for each site and under each condition
    // Ncond * Nsite
    PathSuffStat** suffstatarray;

    // storing cond/site suff stat log probs
    double** sitecondsuffstatlogprob;
    double** bksitecondsuffstatlogprob;

    PhyloProcess* phyloprocess;

  public:
    DiffSelSparseModel(const std::string& datafile, const std::string& treefile, int Ncond,
                       int inNlevel, int infixglob, int infixvar, int incodonmodel, bool sample)
        : Ncond(Ncond) {
        fixglob = infixglob;
        if (!fixglob) {
            cerr << "error: free hyperparameters for baseline (global profile) not yet "
                    "implemented\n";
            exit(1);
        }
        fixvar = infixvar;
        codonmodel = incodonmodel;

        Nlevel = inNlevel;
        if (Nlevel != 2) {
            std::cerr << "-- Error: Nlevel should be equal to 2\n";
            exit(1);
        }

        ReadFiles(datafile, treefile);

        // specifies which condition for which branch
        MakeBranchAllocations();
        std::cerr << "-- conditions over branches ok\n";

        // model allocation
        Allocate();
        std::cerr << "-- model allocation ok\n";

        // unfold phyloprocess (allocate conditional likelihood vectors, etc)
        std::cerr << "-- unfolding\n";
        phyloprocess->Unfold();

        if (sample) {
            // stochastic mapping of substitution histories
            std::cerr << "-- mapping substitutions\n";
            phyloprocess->ResampleSub();
        }
    }

    DiffSelSparseModel(const DiffSelSparseModel&) = delete;

    ~DiffSelSparseModel() {
        // delete tree;
        delete[] branchlength;
        delete[] branchlengthcount;
        delete[] branchlengthbeta;
        delete[] branchalloc;
        delete[] nucrelrate;
        delete[] nucstat;
        delete nucmatrix;
        delete[] condsubmatrixarray;
        delete[] phylosubmatrix;
        delete[] rootsubmatrix;
        delete[] suffstatarray;
        delete[] sitecondsuffstatlogprob;
        delete[] bksitecondsuffstatlogprob;
        delete phyloprocess;
    }

    void ReadFiles(string datafile, string treefile) {
        // nucleotide sequence alignment
        data = new FileSequenceAlignment(datafile);

        // translated into codon sequence alignment
        codondata = new CodonSequenceAlignment(data, true);

        Nsite = codondata->GetNsite();  // # columns
        Ntaxa = codondata->GetNtaxa();

        std::cerr << "-- Number of sites: " << Nsite << std::endl;

        taxonset = codondata->GetTaxonSet();

        // get tree from file (newick format)
        tree = new Tree(treefile);

        // check whether tree and data fits together
        tree->RegisterWith(taxonset);

        // traversal of the tree, so as to number links, branches and nodes
        // convention is: branches start at 1 (branch number 0 is the null branch behind the root)
        // nodes start at 0 (for the root), and nodes 1..Ntaxa are tip nodes (corresponding to taxa
        // in sequence alignment)
        tree->SetIndices();
        Nbranch = tree->GetNbranch();

        std::cerr << "-- Number of taxa : " << Ntaxa << '\n';
        std::cerr << "-- Number of branches : " << Nbranch << '\n';

        std::cerr << "-- Tree and data fit together\n";
    }

    void Allocate() {
        // ----------
        // construction of the model
        // ----------

        // allocating data structures and sampling initial configuration

        // branch lengths

        lambda = 10;
        branchlength = new double[Nbranch];
        for (int j = 0; j < Nbranch; j++) {
            branchlength[j] = Random::sExpo() / lambda;
        }

        branchlengthcount = new int[Nbranch];
        branchlengthbeta = new double[Nbranch];

        // nucleotide matrix

        nucrelrate = new double[Nrr];
        // sample relrate from uniform Dirichlet distribution
        double totrr = 0;
        for (int k = 0; k < Nrr; k++) {
            nucrelrate[k] = Random::sExpo();
            totrr += nucrelrate[k];
        }
        for (int k = 0; k < Nrr; k++) {
            nucrelrate[k] /= totrr;
        }

        nucstat = new double[Nnuc];
        // sample nucstat from uniform Dirichlet distribution
        double totstat = 0;
        for (int k = 0; k < Nnuc; k++) {
            nucstat[k] = Random::sGamma(1.0);
            totstat += nucstat[k];
        }
        for (int k = 0; k < Nnuc; k++) {
            nucstat[k] /= totstat;
        }

        // normalized (true) GTR nucleotide substitution matrix
        nucmatrix = new GTRSubMatrix(Nnuc, nucrelrate, nucstat, true);

        fitness_shape = Random::sExpo();
        InitUniformDirichlet(fitness_inv_rates);

        fitness = std::vector<Eigen::MatrixXd>(Ncond, Eigen::MatrixXd(Nsite, Naa));
        for (auto& G_k : fitness) {
            for (int i = 0; i < Nsite; i++)
                for (int aa = 0; aa < Naa; aa++) {
                    G_k(i, aa) =
                        Random::Gamma(fitness_shape, fitness_shape / fitness_inv_rates[aa]);
                }
        }

        for (int k = 0; k < Ncond; k++) {
            prob_conv[k] = Random::BetaMV(prob_conv_m, prob_conv_v);
        }

        ind_conv = vector<BMatrix>(Ncond, BMatrix(Nsite, Naa));
        for (int k = 0; k < Ncond; k++)
            for (int i = 0; i < Nsite; i++)
                for (int aa = 0; aa < Naa; aa++) {
                    ind_conv[k](i, aa) = (Random::Uniform() < prob_conv[k]);
                }

        // codon matrices
        // per condition and per site
        condsubmatrixarray = new CodonSubMatrix**[Ncond];
        for (int k = 0; k < Ncond; k++) {
            condsubmatrixarray[k] = new CodonSubMatrix*[Nsite];
            for (int i = 0; i < Nsite; i++) {
                if (codonmodel == 0) {
                    condsubmatrixarray[k][i] = new MGSRFitnessSubMatrix(
                        (CodonStateSpace*)codondata->GetStateSpace(), nucmatrix, fitness[0].row(i),
                        fitness[k].row(i), ind_conv[k].row(i), false);
                } else {
                    condsubmatrixarray[k][i] = new MGMSFitnessSubMatrix(
                        (CodonStateSpace*)codondata->GetStateSpace(), nucmatrix, fitness[0].row(i),
                        fitness[k].row(i), ind_conv[k].row(i), false);
                }
            }
        }

        // arrays of pointers for phyloprocess
        // sub matrices per branch and per site
        phylosubmatrix = new SubMatrix**[Nbranch];
        for (int j = 0; j < Nbranch; j++) {
            int k = branchalloc[j];
            phylosubmatrix[j] = new SubMatrix*[Nsite];
            for (int i = 0; i < Nsite; i++) {
                phylosubmatrix[j][i] = condsubmatrixarray[k][i];
            }
        }
        // sub matrices per site (for root equilibrium frequencies)
        rootsubmatrix = new SubMatrix*[Nsite];
        for (int i = 0; i < Nsite; i++) {
            rootsubmatrix[i] = condsubmatrixarray[0][i];
        }

        // create phyloprocess
        phyloprocess =
            new PhyloProcess(tree, codondata, branchlength, 0, phylosubmatrix, 0, rootsubmatrix);

        // create suffstat arrays
        suffstatarray = new PathSuffStat*[Ncond];
        for (int k = 0; k < Ncond; k++) {
            suffstatarray[k] = new PathSuffStat[Nsite];
        }

        // -----
        // auxiliary arrays to speed up computation
        // ----

        // storing suff stat log probs
        sitecondsuffstatlogprob = new double*[Ncond];
        for (int k = 0; k < Ncond; k++) {
            sitecondsuffstatlogprob[k] = new double[Nsite];
        }
        bksitecondsuffstatlogprob = new double*[Ncond];
        for (int k = 0; k < Ncond; k++) {
            bksitecondsuffstatlogprob[k] = new double[Nsite];
        }
    }

    void MakeBranchAllocations() {
        delete[] branchalloc;
        branchalloc = new int[Nbranch];

        // default pre-initialization
        for (int j = 0; j < Nbranch; j++) {
            branchalloc[j] = -1;
        }

        // root allocation
        branchalloc[0] = 0;

        RecursiveMakeBranchAllocations(tree->GetRoot());

        // check that all branches have been correctly initialized
        for (int j = 0; j < Nbranch; j++) {
            if ((branchalloc[j] < 0) || (branchalloc[j] >= Ncond)) {
                std::cerr << "error in make branch allocation\n";
                cerr << j << '\t' << branchalloc[j] << '\n';
                exit(1);
            }
        }
    }

    void RecursiveMakeBranchAllocations(const Link* from) {
        if (!from->isRoot()) {
            int k = atoi(from->GetBranch()->GetName().c_str());
            if (k >= Ncond) {
                k = Ncond - 1;
            }
            if (k < 0) {
                std::cerr << "error : allocation out of bound\n";
                std::cerr << "k" << '\t' << "Ncond" << '\n';
                exit(1);
            }
            branchalloc[from->GetBranch()->GetIndex()] = k;
        }
        for (const Link* link = from->Next(); link != from; link = link->Next()) {
            RecursiveMakeBranchAllocations(link->Out());
        }
    }

    void CorruptNucMatrix() {
        nucmatrix->CopyStationary(nucstat);
        nucmatrix->CorruptMatrix();
    }

    void CorruptCodonMatrices() {
        for (int i = 0; i < Nsite; i++) {
            CorruptSiteCodonMatrices(i);
        }
    }

    void CorruptSiteCodonMatrices(int i) {
        for (int k = 0; k < Ncond; k++) {
            condsubmatrixarray[k][i]->CorruptMatrix();
        }
    }

    void Update() override {
        cerr << "in diffsel update\n";
        CorruptNucMatrix();
        CorruptCodonMatrices();
        phyloprocess->GetLogProb();
    }

    void UpdateAll() {
        CorruptNucMatrix();
        for (int i = 0; i < Nsite; i++) {
            UpdateSite(i);
        }
    }

    void BackupAll() {
        for (int i = 0; i < Nsite; i++) {
            BackupSite(i);
        }
    }

    void RestoreAll() {
        for (int i = 0; i < Nsite; i++) {
            RestoreSite(i);
        }
    }

    void UpdateSite(int i) {
        CorruptSiteCodonMatrices(i);
        for (int k = 0; k < Ncond; k++) {
            sitecondsuffstatlogprob[k][i] = SiteCondSuffStatLogProb(i, k);
        }
    }

    void CheckSite(int i) {
        CorruptSiteCodonMatrices(i);
        for (int k = 0; k < Ncond; k++) {
            if (fabs(sitecondsuffstatlogprob[k][i] - SiteCondSuffStatLogProb(i, k)) > 1e-8) {
                cerr << "error for site " << i << '\n';
                cerr << sitecondsuffstatlogprob[k][i] << SiteCondSuffStatLogProb(i, k) << '\n';
                exit(1);
            }
        }
    }

    void CheckAll() {
        CorruptNucMatrix();
        for (int i = 0; i < Nsite; i++) {
            CheckSite(i);
        }
    }

    void BackupSite(int i) {
        for (int k = 0; k < Ncond; k++) {
            bksitecondsuffstatlogprob[k][i] = sitecondsuffstatlogprob[k][i];
        }
    }

    void RestoreSite(int i) {
        for (int k = 0; k < Ncond; k++) {
            sitecondsuffstatlogprob[k][i] = bksitecondsuffstatlogprob[k][i];
        }
    }


    // ---------------
    // suff stat log probs
    // ---------------

    double GetSuffStatLogProb() {
        double total = 0;
        for (int i = 0; i < Nsite; i++) {
            total += GetSiteSuffStatLogProb(i);
        }
        return total;
    }

    double GetSiteSuffStatLogProb(int i) {
        double total = 0;
        for (int k = 0; k < Ncond; k++) {
            total += sitecondsuffstatlogprob[k][i];
        }
        return total;
    }

    double SiteCondSuffStatLogProb(int i, int k) {
        return condsubmatrixarray[k][i]->SuffStatLogProb(&suffstatarray[k][i]);
    }

    // ---------------
    // log priors
    // ---------------

    double GlobalProfileLogProb() {
        // uniform dirichlet: log prob is constant
        return Nsite * Random::logGamma((double)Naa);
        /*
        double total = 0;
        for (int i=0; i<Nsite; i++)	{
                total += SiteGlobalProfileLogProb(i);
        }
        return total;
        */
    }


    double LambdaLogProb() { return -log(10.0) - lambda / 10; }

    double LengthLogProb() {
        // iid exp of rate lambda: total length is suff stat
        return Nbranch * log(lambda) - lambda * GetTotalLength();
    }

    // ---------------
    // collecting suff stats
    // ---------------

    // suffstats, per condition and per site
    // see SuffStat.hpp
    void CollectSuffStat() {
        ClearSuffStat();
        RecursiveCollectSuffStat(tree->GetRoot());
    }

    void ClearSuffStat() {
        for (int k = 0; k < Ncond; k++) {
            for (int i = 0; i < Nsite; i++) {
                suffstatarray[k][i].Clear();
            }
        }
    }

    void RecursiveCollectSuffStat(const Link* from) {
        if (from->isRoot()) {
            for (int i = 0; i < Nsite; i++) {
                phyloprocess->AddRootSuffStat(i, suffstatarray[0][i]);
            }
        } else {
            for (int i = 0; i < Nsite; i++) {
                phyloprocess->AddSuffStat(
                    i, from, suffstatarray[branchalloc[from->GetBranch()->GetIndex()]][i]);
            }
        }
        for (const Link* link = from->Next(); link != from; link = link->Next()) {
            RecursiveCollectSuffStat(link->Out());
        }
    }

    // length suff stats (count and beta, for each branch)
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
                phyloprocess->AddLengthSuffStat(i, from,
                                                branchlengthcount[from->GetBranch()->GetIndex()],
                                                branchlengthbeta[from->GetBranch()->GetIndex()]);
            }
        }
        for (const Link* link = from->Next(); link != from; link = link->Next()) {
            RecursiveCollectLengthSuffStat(link->Out());
        }
    }

    // still to be implemented
    // more compact suff stats
    // for fitness profiles: collapsed onto the 20 amino-acid states (instead of the 61 codons)
    // for nuc rate parameters: collapsed onto the 4 nucleotide states (and summed over all sites
    // and conditions)
    /*
    void CollectAASuffStat() {}
    void CollectNucSuffStat()	{}
    */

    // move cycle schedule
    // does not yet implement any monitoring (success rates, time spent, etc)
    double Move() override {
        SubMatrix::ResetDiagCount();

        int nrep0 = 3;
        int nrep = 20;

        for (int rep0 = 0; rep0 < nrep0; rep0++) {
            CollectLengthSuffStat();

            MoveBranchLength();
            MoveLambda(1.0, 10);
            MoveLambda(0.3, 10);

            CollectSuffStat();

            UpdateAll();

            for (int rep = 0; rep < nrep; rep++) {
                /* ci gissaient movebaseline, movedelta et move varsel*/
            }

            MoveRR(0.1, 1, 10);
            MoveRR(0.03, 3, 10);
            MoveRR(0.01, 3, 10);

            MoveNucStat(0.1, 1, 10);
            MoveNucStat(0.01, 1, 10);
        }

        UpdateAll();
        // phyloprocess->ResampleSub();
        // fraction of sites that are resampled
        // can be made into a parameter of the mcmc
        double frac = 1;
        phyloprocess->Move(frac);

        return 1.0;
    }

    double MoveFitness(int cond, double tuning, int nrep) {
        int ntot = 0, nacc = 0;
        auto partial_gamma_log_density = [](double alpha, double m, double x) {
            return (alpha - 1.) * log(x) - alpha / m * x;
        };

        for (int rep = 0; rep < nrep; rep++) {
            for (int i = 0; i < Nsite; i++) {
                int aa = Random::Choose(Naa);
                double bk = fitness[cond](i, aa);
                BackupSite(i);

                double loglikelihood_before =
                    partial_gamma_log_density(fitness_shape, fitness_inv_rates[aa], bk) +
                    GetSiteSuffStatLogProb(i);

                fitness[cond](i, aa) *= tuning * exp(Random::Uniform() - 0.5);
                UpdateSite(i);
                double loglikelihood_after =
                    partial_gamma_log_density(fitness_shape, fitness_inv_rates[aa],
                                              fitness[cond](i, aa)) +
                    GetSiteSuffStatLogProb(i);

                double loghastings = 0.;

                double deltalogprob = loglikelihood_after - loglikelihood_before + loghastings;

                int accepted = (log(Random::Uniform()) < deltalogprob);
                if (accepted) {
                    nacc++;
                } else {
                    fitness[cond](i, aa) = bk;
                    RestoreSite(i);
                }
                ntot++;
            }
        }
        return nacc / ntot;
    }

    double MoveFitnessShape(double tuning) {
        auto partial_gamma_log_density = [](double alpha, double m, double x) {
            double beta = alpha / m;
            return alpha * log(beta) - log(tgamma(alpha)) + (alpha - 1) * log(x);
        };

        auto fitness_log_density = [&]() {
            double loglikelihood = -fitness_shape;
            for (auto& G_k : fitness)
                for (int i = 0; i < Nsite; i++)
                    for (int aa = 0; aa < Naa; aa++) {
                        loglikelihood += partial_gamma_log_density(
                            fitness_shape, fitness_inv_rates[aa], G_k(i, aa));
                    }
            return loglikelihood;
        };

        double bk = fitness_shape;

        double loglikelihood_before = fitness_log_density();
        fitness_shape *= tuning * exp(Random::Uniform() - 0.5);
        double loglikelihood_after = fitness_log_density();

        double loghastings = 0.;
        double deltalogprob = loglikelihood_after - loglikelihood_before + loghastings;

        bool accepted = (log(Random::Uniform()) < deltalogprob);
        if (!accepted) {
            fitness_shape = bk;
        }
        return static_cast<bool>(accepted);
    }

    double MoveFitnessInvRates(double tuning) {
        auto partial_gamma_log_density = [](double alpha, double m, double x) {
            double beta = alpha / m;
            return alpha * log(beta) - beta * x;
        };

        auto fitness_log_density = [&]() {
            double loglikelihood = 0;
            for (auto& G_k : fitness)
                for (int i = 0; i < Nsite; i++)
                    for (int aa = 0; aa < Naa; aa++) {
                        loglikelihood += partial_gamma_log_density(
                            fitness_shape, fitness_inv_rates[aa], G_k(i, aa));
                    }
            return loglikelihood;
        };

        AAProfile bk = fitness_inv_rates;

        double loglikelihood_before = fitness_log_density();
        double loghastings = Random::ProfileProposeMove(fitness_inv_rates, tuning, 0);
        double loglikelihood_after = fitness_log_density();

        double deltalogprob = loglikelihood_after - loglikelihood_before + loghastings;

        bool accepted = (log(Random::Uniform()) < deltalogprob);
        if (!accepted) {
            fitness_inv_rates = bk;
        }
        return static_cast<bool>(accepted);
    }

    double MoveProbConv(double tuning) {
        auto partial_beta_log_density = [](double m, double v, double x) {
            double alpha = m / v;
            double beta = (1 - m) / v;
            return (alpha - 1) * log(x) + (beta - 1) * log(1 - x);
        };

        auto ind_conv_log_density = [&]() {
            double loglikelihood = 0;
            for (int k = 1; k < Ncond; k++)
                for (int i = 0; i < Nsite; i++)
                    for (int aa = 0; aa < Naa; aa++) {
                        loglikelihood += (ind_conv[k](i,aa)) ? log(prob_conv[k]) : log(1 - prob_conv[k]);
                    }
            return loglikelihood;
        };

        double ntot = 0, nacc = 0;

        for(int k = 1; k < Ncond; k++) {
            double bk = prob_conv[k];
            double loglikelihood_before =
                ind_conv_log_density()
                + partial_beta_log_density(prob_conv_m, prob_conv_v, prob_conv[k]);

            prob_conv[k] *= tuning * exp(Random::Uniform() - 0.5);
            double loghastings = 0.;

            double loglikelihood_after =
                ind_conv_log_density()
                + partial_beta_log_density(prob_conv_m, prob_conv_v, prob_conv[k]);

            double deltalogprob = loglikelihood_after - loglikelihood_before + loghastings;

            bool accepted = (log(Random::Uniform()) < deltalogprob);
            if (!accepted) {
                prob_conv[k] = bk;
            }
            else nacc++;
            ntot++;
        }
        return nacc / ntot;
    }

    double MoveRR(double tuning, int n, int nrep) {
        double nacc = 0;
        double ntot = 0;
        double bk[Nrr];
        for (int rep = 0; rep < nrep; rep++) {
            // CheckAll();

            for (int l = 0; l < Nrr; l++) {
                bk[l] = nucrelrate[l];
            }
            BackupAll();

            double deltalogprob = -GetSuffStatLogProb();
            double loghastings = Random::ProfileProposeMove(nucrelrate, Nrr, tuning, n);
            deltalogprob += loghastings;

            UpdateAll();

            deltalogprob += GetSuffStatLogProb();

            int accepted = (log(Random::Uniform()) < deltalogprob);
            if (accepted) {
                nacc++;
            } else {
                for (int l = 0; l < Nrr; l++) {
                    nucrelrate[l] = bk[l];
                }
                RestoreAll();
            }
            ntot++;
            // CheckAll();
        }
        return nacc / ntot;
    }

    double MoveNucStat(double tuning, int n, int nrep) {
        double nacc = 0;
        double ntot = 0;
        double bk[Nnuc];
        for (int rep = 0; rep < nrep; rep++) {
            // CheckAll();
            for (int l = 0; l < Nnuc; l++) {
                bk[l] = nucstat[l];
            }
            BackupAll();

            double deltalogprob = -GetSuffStatLogProb();
            double loghastings = Random::ProfileProposeMove(nucstat, Nnuc, tuning, n);
            deltalogprob += loghastings;

            UpdateAll();

            deltalogprob += GetSuffStatLogProb();

            int accepted = (log(Random::Uniform()) < deltalogprob);
            if (accepted) {
                nacc++;
            } else {
                for (int l = 0; l < Nnuc; l++) {
                    nucstat[l] = bk[l];
                }

                RestoreAll();
            }
            ntot++;
            // CheckAll();
        }
        return nacc / ntot;
    }


    double MoveBranchLength() {
        for (int j = 0; j < Nbranch; j++) {
            branchlength[j] =
                Random::Gamma(1.0 + branchlengthcount[j], lambda + branchlengthbeta[j]);
            if (!branchlength[j]) {
                cerr << "error: resampled branch length is 0\n";
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

    // accessors and
    // summary statistics

    int GetNsite() { return Nsite; }
    int GetNcond() { return Ncond; }

    double GetEntropy(double* profile, int dim) {
        double tot = 0;
        for (int i = 0; i < dim; i++) {
            tot -= (profile[i] < 1e-6) ? 0 : profile[i] * log(profile[i]);
        }
        return tot;
    }

    double GetTotalLength() {
        double tot = 0;
        for (int j = 0; j < Nbranch; j++) {
            tot += branchlength[j];
        }
        return tot;
    }


    void TraceHeader(std::ostream& os) override {
        os << "length\t";
        os << "statent\t";
        os << "rrent\t";
        os << "diag\n";
    }

    void Trace(ostream& os) override {
        os << GetTotalLength() << '\t';
        os << GetEntropy(nucstat, Nnuc) << '\t';
        os << GetEntropy(nucrelrate, Nrr) << '\t';
        os << SubMatrix::GetDiagCount() << '\n';
    }

    void Monitor(ostream&) override {}

    void FromStream(istream& is) override {
        is >> lambda;
        for (int i = 0; i < Nbranch; i++) {
            is >> branchlength[i];
        }
        for (int i = 0; i < Nrr; i++) {
            is >> nucrelrate[i];
        }
        for (int i = 0; i < Nnuc; i++) {
            is >> nucstat[i];
        }
    }

    void ToStream(ostream& os) override {
        os << lambda << '\n';
        for (int i = 0; i < Nbranch; i++) {
            os << branchlength[i] << '\t';
        }
        os << '\n';
        for (int i = 0; i < Nrr; i++) {
            os << nucrelrate[i] << '\t';
        }
        os << '\n';
        for (int i = 0; i < Nnuc; i++) {
            os << nucstat[i] << '\t';
        }
        os << '\n';
    }
};
