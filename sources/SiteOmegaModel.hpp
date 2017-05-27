
#include "CodonSequenceAlignment.hpp"
#include "Tree.hpp"
#include "ProbModel.hpp"
#include "GTRSubMatrix.hpp"
#include "CodonSubMatrix.hpp"
#include "PhyloProcess.hpp"

class MGOmegaHeterogeneousCodonSubMatrixArray : public Array<SubMatrix>	{

	public:
	MGOmegaHeterogeneousCodonSubMatrixArray(CodonStateSpace* incodonstatespace, GTRSubMatrix* innucmatrix, const Array<double>* inomegaarray) : codonstatespace(incodonstatespace), nucmatrix(innucmatrix), omegaarray(inomegaarray), matrixarray(inomegaarray->GetSize()) {
		for (int i=0; i<matrixarray.size(); i++)	{
			matrixarray[i] = new MGOmegaCodonSubMatrix(nucmatrix,omegarray->GetVal(i));
		}
	}

	~MGOmegaHeterogeneousCodonSubMatrixArray()	{
		for (int i=0; i<matrixarray.size(); i++)	{
			delete matrixarray[i];
		}
	}
		
	int GetSize() {return matrixarray.size();}

	MGOmegaCodonSubMatrix& GetMGOmegaCodonSubMatrix(int site) const {
		return *matrixarray[site];
	}

	const GTRSubMatrix& GetNucMatrix() const {return *nucmatrix;}

	const SubMatrix& GetVal(int site) const {return *matrixarray[site];}

	private:

	CodonStateSpace* codonstatespace;
	const GTRSubMatrix* nucmatrix;
	const Array<double>* omegaarray;
	vector<MGOmegaCodonSubMatrix*> matrixarray;
};

class SiteOmegaModel	{

	double lambda;
	BranchIIDGamma* branchlength;
	
	double alpha;
	double beta;
	IIDGamma* omegaarray;

	vector<double> nucrelrate;
	vector<double> nucstat;
	GTRSubMatrix* nucmatrix;

	MGOmegaCodonSubMatrixArray* codonmatrixarray;
	BranchHomogeneousSiteHeterogeneousArray<SubMatrix>* matrixsystem;

	// suffstats

	BranchPoissonSuffStatArray* lengthsuffstat;
	PathSuffStatArray* pathsuffstat;
	PoissonSuffStatArray* omegasuffstat;
	NucPathSuffStat nucsuffstat;	

	SiteOmegaModel()	{

	}

	void MakeModel()	{

		lambda = 10.0;
		branchlength = new BranchIIDGamma(tree,1.0,lambda);

		alpha = beta = 1.0;
		omegaarray = new IIDGamma(GetNsite(),alpha,beta);

		nucrelrate ...
		nucstat ...
		nucmatrix = neew GTRSubMatrix(Nnuc,nucrelrate,nucstat,true);

		codonmatrixarray = new MGOmegaHeterogeneousCodonSubMatrixArray(codonstatespace,nucmatrix,omegaarray);

		matrixsystem = new BranchHomogeneousSiteHeterogeneousArray<SubMatrix>(codonmatrixarray);

		phyloprocess = new PhyloProcess(tree,data,branchlength,0,matrixsystem,0,codonmatrixarray);

		lengthsuffstatarray = new BranchPoissonSuffStatArray(tree);
		pathsuffstatarray = new PathSuffStatArray(GetNsite());
		omegasuffstatarray = new PoissonSuffStatArray(GetNsite());
	}

	void Move()	{

		phyloprocess->ResampleSub();

		lengthsuffstatarray->Clear();
		phyloprocess->AddLengthSuffStat(lengthsuffstatarray);
		branchlength->GibbsResample(lengthsuffstatarray);
		MoveLambda();

		pathsuffstatarray->Clear();
		phyloprocess->AddPathSuffStat(pathsuffstatarray);
		omegasuffstatarray->Clear();
		omegasuffstatarray->AddSuffStat(codonstatespace,pathsuffstatarray);
		omegaarray->GibbsResample(omegasuffstatarray);

		MoveAlpha();
		MoveBeta();

		nucsuffstat->Clear();
		nucsuffstat->Add(codonstatespace,pathsuffstatarray);
		MoveNucRates();
	}

};

