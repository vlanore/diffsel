
#include "Sample.hpp"
#include "DiffSelModel.hpp"

class DiffSelSample : public Sample	{

  private:
    // Sample parameters
    string modeltype, datafile, treefile;
    int codonmodel, category, level, fixglob, fixvar; 

	public:

	string GetModelType() {return modeltype;}

	DiffSelModel* GetModel() {return (DiffSelModel*) model;}

	DiffSelSample(string filename, int inburnin, int inevery, int inuntil) : Sample(filename,inburnin,inevery,inuntil)	{
		Open();
	}

	void Open()	{

		// open <name>.param
		ifstream is((name + ".param").c_str());

		// check that file exists
		if (!is)	{
			cerr << "error : cannot find file : " << name << ".param\n";
			exit(1);
		}

		// read model type, and other standard fields
		is >> modeltype;
		is >> datafile >> treefile >> category >> level;
		is >> fixglob >> fixvar >> codonmodel;
		int tmp;
		is >> tmp;
		if (tmp) {
		    cerr << "-- Error when reading model\n";
		    exit(1);
		}
		is >> chainevery >> chainuntil >> chainsize;
		// the chain's saving frequency, upper limit and current size
		// not to be confused with the sample's subsampling frequency, upper limit and size

		// make a new model depending on the type obtained from the file
		if (modeltype == "DIFFSEL")	{
			cerr << "CREATE\n";
			model = new DiffSelModel(datafile, treefile, category, level, fixglob, fixvar, codonmodel, false);
		}
		else	{
			cerr << "error when opening file "  << name << '\n';
			exit(1);
		}

		// read model (i.e. chain's last point) from <name>.param
		cerr << "from stream\n";
		model->FromStream(is);
		// cerr << "UPDATE\n";
		// model->Update();
		// open <name>.chain, and prepare stream and stream iterator
		OpenChainFile();
		// now, size is defined (it is the total number of points with which this Sample object will make all its various posterior averages)
		// all these points can be accessed to (only once) by repeated calls to GetNextPoint()
	}

	void Read(double cutoff)	{

		// prepare the mean and variance
		
		int Nsite = GetModel()->GetNsite();
		int K = GetModel()->GetNcond();

		double*** meansel = new double**[K];
		for (int k=1; k<K; k++)	{
			meansel[k] = new double*[Nsite];
			for (int i=0; i<Nsite; i++)	{
				meansel[k][i] = new double[Naa];
				for (int a=0; a<Naa; a++)	{
					meansel[k][i][a] = 0;
				}
			}
		}

		double*** ppsel = new double**[K];
		for (int k=1; k<K; k++)	{
			ppsel[k] = new double*[Nsite];
			for (int i=0; i<Nsite; i++)	{
				ppsel[k][i] = new double[Naa];
				for (int a=0; a<Naa; a++)	{
					ppsel[k][i][a] = 0;
				}
			}
		}

		// cycle over the sample
		for (int c=0; c<size; c++)	{
			GetNextPoint();
			cerr << '.';

			for (int k=1; k<K; k++)	{
				for (int i=0; i<Nsite; i++)	{
					double* delta = GetModel()->GetDelta(k,i);

					double mean = 0;
					for (int a=0; a<Naa; a++)	{
						mean += delta[a];
					}
					mean /= Naa;

					for (int a=0; a<Naa; a++)	{
						double tmp = delta[a] - mean;
						meansel[k][i][a] += tmp;
						if (tmp > 0)	{
							ppsel[k][i][a] ++;
						}
					}
				}
			}
		}
		cerr << '\n';

		for (int k=1; k<K; k++)	{

			ostringstream s1,s2;
			s1 << GetName() << "_" << k << ".meandiffsel";
			s2 << GetName() << "_" << k << ".signdiffsel";
			ofstream os1(s1.str().c_str());
			ofstream os2(s2.str().c_str());

			for (int i=0; i<Nsite; i++)	{

				os1 << i;
				for (int a=0; a<Naa; a++)	{
					ppsel[k][i][a] /= size;
					os1 << '\t' << ppsel[k][i][a];
					
				}
				for (int a=0; a<Naa; a++)	{
					meansel[k][i][a] /= size;
					os1 << '\t' << meansel[k][i][a];
				}
				os1 << '\n';

				for (int a=0; a<Naa; a++)	{
					if (ppsel[k][i][a] > cutoff)	{
						os2 << i << '\t' << a << '\t' << ppsel[k][i][a] << '\t' << meansel[k][i][a] << '\n';
					}
				}
			}
		}
	}
};

int main(int argc, char* argv[])	{

	int burnin = 0;
	int every = 1;
	int until = -1;
	string name;

	double cutoff = 0.9;

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;

		while (i < argc)	{
			string s = argv[i];
			if (s == "-c")	{
				i++;
				cutoff = atof(argv[i]);
			}
			else if ( (s == "-x") || (s == "-extract") )	{
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				burnin = atoi(argv[i]);
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				every = atoi(argv[i]);
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				until = atoi(argv[i]);
			}
			else	{
					if (i != (argc -1))	{
					throw(0);
				}
				name = argv[i];
			}
			i++;
		}
		if (name == "")	{
			throw(0);
		}
		
	}
	catch(...)	{
		cerr << "readpb [-x <burnin> <every> <until>] <chainname> \n";
		cerr << '\n';
		exit(1);
	}

	DiffSelSample sample(name,burnin,every,until);
	sample.Read(cutoff);

}

