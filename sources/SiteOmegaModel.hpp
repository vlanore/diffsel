
#include "CodonSequenceAlignment.hpp"
#include "Tree.hpp"
#include "ProbModel.hpp"
#include "GTRSubMatrix.hpp"
#include "CodonSubMatrix.hpp"
#include "PhyloProcess.hpp"
#include "IIDGamma.hpp"
#include "CodonSuffStat.hpp"
#include "CodonSubMatrixArray.hpp"

const int Nrr = Nnuc * (Nnuc-1) / 2;
const int Nstate = 61;

class SiteOmegaModel	{

	Tree* tree;
	FileSequenceAlignment* data;
	const TaxonSet* taxonset;
	CodonSequenceAlignment* codondata;

	int Nsite;
	int Ntaxa;
	int Nbranch;

	double lambda;
	BranchIIDGamma* branchlength;
	
	double alpha;
	double beta;
	IIDGamma* omegaarray;

	vector<double> nucrelrate;
	vector<double> nucstat;
	GTRSubMatrix* nucmatrix;

	MGOmegaHeterogeneousCodonSubMatrixArray* codonmatrixarray;
	// BranchHomogeneousSiteHeterogeneousArray<SubMatrix>* matrixsystem;

    PhyloProcess* phyloprocess;

	// suffstats

	PoissonSuffStatBranchArray* lengthsuffstatarray;
	GammaSuffStat lambdasuffstat;
	OmegaSuffStatArray* omegasuffstatarray;
	PathSuffStatArray* pathsuffstatarray;
	// NucPathSuffStat nucsuffstat;	

	double suffstatlogprob;
	double bksuffstatlogprob;

    public:

	SiteOmegaModel(string datafile, string treefile)	{

		data = new FileSequenceAlignment(datafile);
		codondata = new CodonSequenceAlignment(data, true);

		Nsite = codondata->GetNsite();    // # columns
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
		cerr << "-- unfold\n";
		phyloprocess->Unfold();
		cerr << phyloprocess->GetLogProb() << '\n';
		std::cerr << "-- mapping substitutions\n";
		phyloprocess->ResampleSub();
		Trace(cerr);
	}

    int GetNsite() {return data->GetNsite();}

	void Allocate()	{

		lambda = 10.0;
		branchlength = new BranchIIDGamma(tree,1.0,lambda);

		alpha = beta = 1.0;
		omegaarray = new IIDGamma(GetNsite(),alpha,beta);

		nucrelrate.assign(Nrr,0);
		double totrr = 0;
		for (int k=0; k<Nrr; k++)	{
			nucrelrate[k] = Random::sExpo();
			totrr += nucrelrate[k];
		}
		for (int k=0; k<Nrr; k++)	{
			nucrelrate[k] /= totrr;
		}

		nucstat.assign(Nnuc,0);
		double totstat = 0;
		for (int k=0; k<Nnuc; k++)	{
			nucstat[k] = Random::sGamma(1.0);
			totstat += nucstat[k];
		}
		for (int k=0; k<Nnuc; k++)	{
			nucstat[k] /= totstat;
		}
		nucmatrix = new GTRSubMatrix(Nnuc,nucrelrate,nucstat,true);

		codonmatrixarray = new MGOmegaHeterogeneousCodonSubMatrixArray((CodonStateSpace*) codondata->GetStateSpace(),nucmatrix,omegaarray);

		phyloprocess = new PhyloProcess(tree,data,branchlength,0,codonmatrixarray);

		lengthsuffstatarray = new PoissonSuffStatBranchArray(tree);
		pathsuffstatarray = new PathSuffStatArray(GetNsite());
		omegasuffstatarray = new OmegaSuffStatArray(GetNsite());
	}

	void UpdateNucMatrix()	{
		nucmatrix->CopyStationary(nucstat);
		nucmatrix->CorruptMatrix();
	}

	void UpdateCodonMatrices()	{
		codonmatrixarray->UpdateCodonMatrices();
	}
		
	void UpdateSuffStatLogProb()	{
		UpdateNucMatrix();
		UpdateCodonMatrices();
		suffstatlogprob = pathsuffstatarray->GetLogProb(codonmatrixarray);
	}

	double GetSuffStatLogProb()	{
		return suffstatlogprob;
	}

	void BackupSuffStatLogProb()	{
		bksuffstatlogprob = suffstatlogprob;
	}

	void RestoreSuffStatLogProb()	{
		suffstatlogprob = bksuffstatlogprob;
	}

	// exponential of mean 1
	double OmegaLogProb()	{
		return omegaarray->GetLogProb();
	}

    double AlphaLogProb()   {
        return -alpha;
    }

    double BetaLogProb()    {
        return -beta;
    }

	double LambdaLogProb()	{
		return -lambda / 10;
	}

	double LengthSuffStatLogProb()	{
		return lambdasuffstat.GetLogProb(1.0,lambda);
	}

	double LengthLogProb()	{
		return branchlength->GetLogProb();
	}

	void Move()	{

		phyloprocess->ResampleSub();

		int nrep = 30;

		for (int rep=0; rep<nrep; rep++)	{

			ResampleBranchLengths();
			MoveLambda();

			CollectPathSuffStat();

			MoveOmega();
            /*
            MoveAlpha();
            MoveBeta();
            */

			UpdateSuffStatLogProb();

			MoveRR(0.1,1,3);
			MoveRR(0.03,3,3);
			MoveRR(0.01,3,3);

			MoveNucStat(0.1,1,3);
			MoveNucStat(0.01,1,3);

			UpdateSuffStatLogProb();

		}
	}

	void ResampleBranchLengths()	{

		lengthsuffstatarray->Clear();
		phyloprocess->AddLengthSuffStat(lengthsuffstatarray);
		branchlength->GibbsResample(lengthsuffstatarray);
	}

	void MoveLambda()	{

		lambdasuffstat.Clear();
		branchlength->AddSuffStat(lambdasuffstat);
		MoveLambda(1.0,10);
		MoveLambda(0.3,10);
		branchlength->SetScale(lambda);
	}

	void CollectPathSuffStat()	{

		pathsuffstatarray->Clear();
		phyloprocess->AddPathSuffStat(pathsuffstatarray);
		UpdateSuffStatLogProb();
	}

	void MoveOmega()	{

		omegasuffstatarray->Clear();
		omegasuffstatarray->AddSuffStat(*codonmatrixarray,*pathsuffstatarray);
        omegaarray->GibbsResample(omegasuffstatarray);
		UpdateCodonMatrices();
	}

	double MoveRR(double tuning, int n, int nrep)	{
		double nacc = 0;
		double ntot = 0;
		double bk[Nrr];
		for (int rep=0; rep<nrep; rep++)	{
			for (int l=0; l<Nrr; l++)	{
				bk[l] = nucrelrate[l];
			}
			BackupSuffStatLogProb();
			double deltalogprob = -GetSuffStatLogProb();
			double loghastings = Random::ProfileProposeMove(nucrelrate,Nrr,tuning,n);
			deltalogprob += loghastings;
			UpdateSuffStatLogProb();
			deltalogprob += GetSuffStatLogProb();
			int accepted = (log(Random::Uniform()) < deltalogprob);
			if (accepted)	{
				nacc ++;
			}
			else	{
				for (int l=0; l<Nrr; l++)	{
					nucrelrate[l] = bk[l];
				}
				RestoreSuffStatLogProb();
			}
			ntot++;
		}
		return nacc/ntot;
	}

	double MoveNucStat(double tuning, int n, int nrep)	{
		double nacc = 0;
		double ntot = 0;
		double bk[Nnuc];
		for (int rep=0; rep<nrep; rep++)	{
			for (int l=0; l<Nnuc; l++)	{
				bk[l] = nucstat[l];
			}
			BackupSuffStatLogProb();
			double deltalogprob = -GetSuffStatLogProb();
			double loghastings = Random::ProfileProposeMove(nucstat,Nnuc,tuning,n);
			deltalogprob += loghastings;
			UpdateSuffStatLogProb();
			deltalogprob += GetSuffStatLogProb();
			int accepted = (log(Random::Uniform()) < deltalogprob);
			if (accepted)	{
				nacc ++;
			}
			else	{
				for (int l=0; l<Nnuc; l++)	{
					nucstat[l] = bk[l];
				}
				RestoreSuffStatLogProb();
			}
			ntot++;
		}
		return nacc/ntot;
	}

	double MoveLambda(double tuning, int nrep)	{

		double nacc = 0;
		double ntot = 0;
		for (int rep=0; rep<nrep; rep++)	{
			double deltalogprob = - LambdaLogProb() - LengthSuffStatLogProb();
			double m = tuning * (Random::Uniform() - 0.5);
			double e = exp(m);
			lambda *= e;
			deltalogprob += LambdaLogProb() + LengthSuffStatLogProb();
			deltalogprob += m;
			int accepted = (log(Random::Uniform()) < deltalogprob);
			if (accepted)	{
				nacc ++;
			}
			else	{
				lambda /= e;
			}
			ntot++;
		}
		return nacc/ntot;
	}

	// summary statistics

	double GetLogPrior() {
		double total = 0;
		total += LambdaLogProb();
		total += LengthLogProb();
		total += OmegaLogProb();
		return total;
	}

	double GetLogLikelihood()	{
		return phyloprocess->GetLogProb();
	}

	double GetTotalLength()	{
		double tot = 0;
		for (int j=1; j<Nbranch; j++)	{
			tot += branchlength->GetVal(j);
		}
		return tot;
	}

    double GetMeanOmega()   {
        return 0;
    }

    double GetVarOmega()    {
        return 0;
    }

	double GetEntropy(const std::vector<double>& profile, int dim) const {
		double tot = 0;
		for (int i=0; i<dim; i++)	{
			tot -= (profile[i] < 1e-6) ? 0 : profile[i]*log(profile[i]);
		}
		return tot;
	}

	void TraceHeader(std::ostream& os)  {
		os << "#logprior\tlnL\tlength\tlambda\t";
		os << "meanomega\tvaromega\talpha\tbeta\t";
		os << "statent\t";
		os << "rrent\n";
	}

	void Trace(ostream& os) {	
		os << GetLogPrior() << '\t';
		os << GetLogLikelihood() << '\t';
		os << GetTotalLength() << '\t';
		os << lambda << '\t';
		os << GetMeanOmega() << '\t';
        os << GetVarOmega() << '\t';
        os << alpha << '\t' << beta << '\t';
		os << GetEntropy(nucstat,Nnuc) << '\t';
		os << GetEntropy(nucrelrate,Nrr) << '\n';
	}

	void Monitor(ostream& os) {}

	void FromStream(istream& is) {}
	void ToStream(ostream& os) {}

};

