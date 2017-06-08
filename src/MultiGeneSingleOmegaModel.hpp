
#include "CodonSequenceAlignment.hpp"
#include "CodonSubMatrix.hpp"
#include "GTRSubMatrix.hpp"
#include "PhyloProcess.hpp"
#include "ProbModel.hpp"
#include "Tree.hpp"

const int Nrr = Nnuc * (Nnuc - 1) / 2;
const int Nstate = 61;

class MultiGeneSingleOmegaModel : public ProbModel {
    Tree* tree;
    FileSequenceAlignment** data;
    const TaxonSet* taxonset;
    CodonSequenceAlignment** codondata;

    int Ngene;
    int* Nsite;
    int Ntaxa;
    int Nbranch;

    double lambda;
    double* branchlength;
    int* branchlengthcount;
    double* branchlengthbeta;

    int** genebranchlengthcount;
    double** genebranchlengthbeta;

    double** nucstat;
    double** nucrelrate;
    GTRSubMatrix** nucmatrix;

    // gene specific omegas iid gamma(alpha,beta)
    double* omega;
    double alpha;
    double beta;

    MGOmegaCodonSubMatrix** codonmatrix;

    SubMatrix**** phylosubmatrix;
    SubMatrix*** rootsubmatrix;

    PhyloProcess** phyloprocess;

    PathSuffStat* suffstat;
    double* suffstatlogprob;
    double* bksuffstatlogprob;

  public:
    MultiGeneSingleOmegaModel(string genelistfile, string treefile) {
        ifstream is(genelistfile.c_str());
        is >> Ngene;
        data = new FileSequenceAlignment*[Ngene];
        codondata = new CodonSequenceAlignment*[Ngene];
        Nsite = new int[Ngene];

        std::cerr << "-- Number of sites: \n";
        for (int gene = 0; gene < Ngene; gene++) {
            string datafile;
            is >> datafile;
            data[gene] = new FileSequenceAlignment(datafile);
            codondata[gene] = new CodonSequenceAlignment(data[gene], true);

            Nsite[gene] = codondata[gene]->GetNsite();  // # columns
            std::cerr << datafile << '\t' << Nsite[gene] << '\n';
        }
        Ntaxa = codondata[0]->GetNtaxa();

        taxonset = codondata[0]->GetTaxonSet();

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
        for (int gene = 0; gene < Ngene; gene++) {
            cerr << "-- unfold\n";
            phyloprocess[gene]->Unfold();
            std::cerr << "-- mapping substitutions\n";
            phyloprocess[gene]->ResampleSub();
            std::cerr << "-- collect suffstat\n";
            CollectSuffStat(gene);
        }
        Trace(cerr);
    }

    void Allocate() {
        lambda = 10;
        branchlength = new double[Nbranch];
        for (int j = 0; j < Nbranch; j++) {
            branchlength[j] = Random::sExpo() / lambda;
        }

        branchlengthcount = new int[Nbranch];
        branchlengthbeta = new double[Nbranch];

        genebranchlengthcount = new int*[Ngene];
        genebranchlengthbeta = new double*[Ngene];
        for (int gene = 0; gene < Ngene; gene++) {
            genebranchlengthcount[gene] = new int[Nbranch];
            genebranchlengthbeta[gene] = new double[Nbranch];
        }

        nucrelrate = new double*[Ngene];
        nucstat = new double*[Ngene];
        nucmatrix = new GTRSubMatrix*[Ngene];
        omega = new double[Ngene];
        codonmatrix = new MGOmegaCodonSubMatrix*[Ngene];
        phylosubmatrix = new SubMatrix***[Ngene];
        rootsubmatrix = new SubMatrix**[Ngene];
        phyloprocess = new PhyloProcess*[Ngene];

        alpha = beta = 1.0;
        for (int gene = 0; gene < Ngene; gene++) {
            nucrelrate[gene] = new double[Nrr];
            double totrr = 0;
            for (int k = 0; k < Nrr; k++) {
                nucrelrate[gene][k] = Random::sExpo();
                totrr += nucrelrate[gene][k];
            }
            for (int k = 0; k < Nrr; k++) {
                nucrelrate[gene][k] /= totrr;
            }

            nucstat[gene] = new double[Nnuc];
            double totstat = 0;
            for (int k = 0; k < Nnuc; k++) {
                nucstat[gene][k] = Random::sGamma(1.0);
                totstat += nucstat[gene][k];
            }
            for (int k = 0; k < Nnuc; k++) {
                nucstat[gene][k] /= totstat;
            }

            nucmatrix[gene] = new GTRSubMatrix(Nnuc, nucrelrate[gene], nucstat[gene], true);
            omega[gene] = 1.0;
            codonmatrix[gene] = new MGOmegaCodonSubMatrix(
                (CodonStateSpace*)codondata[gene]->GetStateSpace(), nucmatrix[gene], omega[gene]);

            // codon matrices
            // per branch and per site
            // (array of ptrs based on condsubmatrixarray)
            phylosubmatrix[gene] = new SubMatrix**[Nbranch];
            for (int j = 0; j < Nbranch; j++) {
                phylosubmatrix[gene][j] = new SubMatrix*[Nsite[gene]];
                for (int i = 0; i < Nsite[gene]; i++) {
                    phylosubmatrix[gene][j][i] = codonmatrix[gene];
                }
            }

            rootsubmatrix[gene] = new SubMatrix*[Nsite[gene]];
            for (int i = 0; i < Nsite[gene]; i++) {
                rootsubmatrix[gene][i] = codonmatrix[gene];
            }

            // phyloprocess
            phyloprocess[gene] = new PhyloProcess(tree, codondata[gene], branchlength, 0,
                                                  phylosubmatrix[gene], 0, rootsubmatrix[gene]);

            suffstat = new PathSuffStat[Ngene];
            suffstatlogprob = new double[Ngene];
            bksuffstatlogprob = new double[Ngene];
        }
    }

    void UpdateNucMatrix(int gene) {
        nucmatrix[gene]->CopyStationary(nucstat[gene]);
        nucmatrix[gene]->CorruptMatrix();
    }

    void UpdateCodonMatrix(int gene) {
        codonmatrix[gene]->SetOmega(omega[gene]);
        codonmatrix[gene]->CorruptMatrix();
    }

    void CollectSuffStat(int gene) {
        suffstat[gene].Clear();
        RecursiveCollectSuffStat(gene, tree->GetRoot());
    }

    void RecursiveCollectSuffStat(int gene, const Link* from) {
        if (from->isRoot()) {
            for (int i = 0; i < Nsite[gene]; i++) {
                phyloprocess[gene]->AddRootSuffStat(i, suffstat[gene]);
            }
        } else {
            for (int i = 0; i < Nsite[gene]; i++) {
                phyloprocess[gene]->AddSuffStat(i, from->Out(), suffstat[gene]);
            }
        }
        for (const Link* link = from->Next(); link != from; link = link->Next()) {
            RecursiveCollectSuffStat(gene, link->Out());
        }
    }

    void CollectLengthSuffStat() {
        for (int j = 0; j < Nbranch; j++) {
            branchlengthcount[j] = 0;
            branchlengthbeta[j] = 0;
            for (int gene = 0; gene < Ngene; gene++) {
                CollectLengthSuffStat(gene);
                branchlengthcount[j] += genebranchlengthcount[gene][j];
                branchlengthbeta[j] += genebranchlengthbeta[gene][j];
            }
        }
    }

    void CollectLengthSuffStat(int gene) {
        ClearLengthSuffStat(gene);
        RecursiveCollectLengthSuffStat(gene, tree->GetRoot());
    }

    void ClearLengthSuffStat(int gene) {
        for (int j = 0; j < Nbranch; j++) {
            genebranchlengthcount[gene][j] = 0;
            genebranchlengthbeta[gene][j] = 0;
        }
    }

    void RecursiveCollectLengthSuffStat(int gene, const Link* from) {
        if (!from->isRoot()) {
            for (int i = 0; i < Nsite[gene]; i++) {
                phyloprocess[gene]->AddLengthSuffStat(
                    i, from->Out(), genebranchlengthcount[gene][from->GetBranch()->GetIndex()],
                    genebranchlengthbeta[gene][from->GetBranch()->GetIndex()]);
            }
        }
        for (const Link* link = from->Next(); link != from; link = link->Next()) {
            RecursiveCollectLengthSuffStat(gene, link->Out());
        }
    }

    void UpdateSuffStatLogProb(int gene) {
        UpdateNucMatrix(gene);
        UpdateCodonMatrix(gene);
        suffstatlogprob[gene] = codonmatrix[gene]->SuffStatLogProb(&suffstat[gene]);
    }

    double GetSuffStatLogProb(int gene) { return suffstatlogprob[gene]; }

    void BackupSuffStatLogProb(int gene) { bksuffstatlogprob[gene] = suffstatlogprob[gene]; }

    void RestoreSuffStatLogProb(int gene) { suffstatlogprob[gene] = bksuffstatlogprob[gene]; }

    double OmegaLogProb() {
        double total = 0;
        for (int gene = 0; gene < Ngene; gene++) {
            total += OmegaLogProb(gene);
        }
        return total;
    }

    double OmegaLogProb(int gene) {
        return alpha * log(beta) - Random::logGamma(alpha) + (alpha - 1) * log(omega[gene]) -
               beta * omega[gene];
    }

    double AlphaLogProb() { return 0; }

    double BetaLogProb() { return 0; }

    double LambdaLogProb() { return -lambda / 10; }

    double LengthLogProb() { return Nbranch * log(lambda) - lambda * GetTotalLength(); }

    double Move() {
        for (int gene = 0; gene < Ngene; gene++) {
            phyloprocess[gene]->ResampleSub();
        }

        int nrep = 30;

        for (int rep = 0; rep < nrep; rep++) {
            CollectLengthSuffStat();
            MoveBranchLength();
            MoveLambda(1.0, 10);
            MoveLambda(0.3, 10);

            for (int gene = 0; gene < Ngene; gene++) {
                CollectSuffStat(gene);
                UpdateSuffStatLogProb(gene);

                MoveRR(gene, 0.1, 1, 3);
                MoveRR(gene, 0.03, 3, 3);
                MoveRR(gene, 0.01, 3, 3);

                MoveNucStat(gene, 0.1, 1, 3);
                MoveNucStat(gene, 0.01, 1, 3);

                MoveOmega(gene, 0.3, 3);
                MoveOmega(gene, 0.1, 3);
                MoveOmega(gene, 0.03, 3);
            }
            MoveAlpha(1, 3);
            MoveAlpha(0.3, 3);
            MoveBeta(1, 3);
            MoveBeta(0.3, 3);
        }
        return 1.0;
    }

    double MoveRR(int gene, double tuning, int n, int nrep) {
        double nacc = 0;
        double ntot = 0;
        double bk[Nrr];
        for (int rep = 0; rep < nrep; rep++) {
            for (int l = 0; l < Nrr; l++) {
                bk[l] = nucrelrate[gene][l];
            }
            BackupSuffStatLogProb(gene);
            double deltalogprob = -GetSuffStatLogProb(gene);
            double loghastings = Random::ProfileProposeMove(nucrelrate[gene], Nrr, tuning, n);
            deltalogprob += loghastings;
            UpdateSuffStatLogProb(gene);
            deltalogprob += GetSuffStatLogProb(gene);
            int accepted = (log(Random::Uniform()) < deltalogprob);
            if (accepted) {
                nacc++;
            } else {
                for (int l = 0; l < Nrr; l++) {
                    nucrelrate[gene][l] = bk[l];
                }
                RestoreSuffStatLogProb(gene);
            }
            ntot++;
        }
        return nacc / ntot;
    }

    double MoveNucStat(int gene, double tuning, int n, int nrep) {
        double nacc = 0;
        double ntot = 0;
        double bk[Nnuc];
        for (int rep = 0; rep < nrep; rep++) {
            for (int l = 0; l < Nnuc; l++) {
                bk[l] = nucstat[gene][l];
            }
            BackupSuffStatLogProb(gene);
            double deltalogprob = -GetSuffStatLogProb(gene);
            double loghastings = Random::ProfileProposeMove(nucstat[gene], Nnuc, tuning, n);
            deltalogprob += loghastings;
            UpdateSuffStatLogProb(gene);
            deltalogprob += GetSuffStatLogProb(gene);
            int accepted = (log(Random::Uniform()) < deltalogprob);
            if (accepted) {
                nacc++;
            } else {
                for (int l = 0; l < Nnuc; l++) {
                    nucstat[gene][l] = bk[l];
                }
                RestoreSuffStatLogProb(gene);
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

    double MoveOmega(int gene, double tuning, int nrep) {
        double nacc = 0;
        double ntot = 0;
        for (int rep = 0; rep < nrep; rep++) {
            BackupSuffStatLogProb(gene);
            double deltalogprob = -OmegaLogProb(gene) - GetSuffStatLogProb(gene);
            double m = tuning * (Random::Uniform() - 0.5);
            double e = exp(m);
            omega[gene] *= e;
            UpdateSuffStatLogProb(gene);
            deltalogprob += OmegaLogProb(gene) + GetSuffStatLogProb(gene);
            deltalogprob += m;
            int accepted = (log(Random::Uniform()) < deltalogprob);
            if (accepted) {
                nacc++;
            } else {
                omega[gene] /= e;
                RestoreSuffStatLogProb(gene);
            }
            ntot++;
        }
        return nacc / ntot;
    }

    double MoveAlpha(double tuning, int nrep) {
        double nacc = 0;
        double ntot = 0;
        for (int rep = 0; rep < nrep; rep++) {
            double deltalogprob = -AlphaLogProb() - OmegaLogProb();
            double m = tuning * (Random::Uniform() - 0.5);
            double e = exp(m);
            alpha *= e;
            deltalogprob += AlphaLogProb() + OmegaLogProb();
            deltalogprob += m;
            int accepted = (log(Random::Uniform()) < deltalogprob);
            if (accepted) {
                nacc++;
            } else {
                alpha /= e;
            }
            ntot++;
        }
        return nacc / ntot;
    }

    double MoveBeta(double tuning, int nrep) {
        double nacc = 0;
        double ntot = 0;
        for (int rep = 0; rep < nrep; rep++) {
            double deltalogprob = -BetaLogProb() - OmegaLogProb();
            double m = tuning * (Random::Uniform() - 0.5);
            double e = exp(m);
            beta *= e;
            deltalogprob += BetaLogProb() + OmegaLogProb();
            deltalogprob += m;
            int accepted = (log(Random::Uniform()) < deltalogprob);
            if (accepted) {
                nacc++;
            } else {
                beta /= e;
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

    double GetLogLikelihood() {
        double total = 0;
        for (int gene = 0; gene < Ngene; gene++) {
            total += phyloprocess[gene]->GetLogProb();
        }
        return total;
    }

    double GetEntropy(double* profile, int dim) {
        double tot = 0;
        for (int i = 0; i < dim; i++) {
            tot -= (profile[i] < 1e-6) ? 0 : profile[i] * log(profile[i]);
        }
        return tot;
    }

    double GetMeanOmega() {
        double total = 0;
        for (int gene = 0; gene < Ngene; gene++) {
            total += omega[gene];
        }
        return total / Ngene;
    }

    void TraceHeader(std::ostream& os) {
        os << "#logprior\tlnL\tlength\tlambda\t";
        os << "omega\t";
        os << "alpha\tbeta\n";
    }

    void Trace(ostream& os) {
        os << GetLogPrior() << '\t';
        os << GetLogLikelihood() << '\t';
        os << GetTotalLength() << '\t';
        os << lambda << '\t';
        os << GetMeanOmega() << '\t';
        os << alpha << '\t' << beta << '\n';
    }

    void Monitor(ostream&) {}

    void FromStream(istream&) {}
    void ToStream(ostream&) {}
};
