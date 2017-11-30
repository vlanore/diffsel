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

#ifndef SAMPLE_H
#define SAMPLE_H

#include <fstream>
#include <string>

class ProbModel;

// Sample is in charge of reading chains from file,
// and computing posterior averages or distributions
class Sample {
  public:
    /// opening chain from file
    Sample(std::string filename, int in_burnin = 0, int in_every = 1, int in_until = -1);

    virtual ~Sample() = default;

    // return the base for all file names (same as for the chain object)
    std::string GetName() { return name; }

    // get the next point (automatically accounts for subsampling, as specified by the "every"
    // parameter)
    // this point is accessed to using GetModel()
    void GetNextPoint();

    // returns a pointer to current point
    // in general this function will be overriden by derived classes
    // because derived classes manipulate models that are derived from ProbModel
    virtual ProbModel* GetModel() { return model; }

    // std::string meant as a check when opening files (model name)
    // you should give different names to chains based on different models
    virtual std::string GetModelType() = 0;

    // abstract method: is responsible for model type checking, and model creation
    // then, calls OpenChainFile()
    virtual void Open() = 0;

    // opens files, and prepare data structures
    // sets the stream iterator at the right point (i.e. discards the burnin)
    // after OpenChainFile() has been called,
    // size is defined (it is the total number of points with which this Sample object will make all
    // its various posterior averages)
    // all these points can be accessed to (only once) by repeated calls to GetNextPoint()
    void OpenChainFile();

    int size;  // sample size (calculated from parameters above)

  protected:
    std::ifstream chain_is;
    int chainevery{-1};  // chain's saving frequency
    int chainuntil{-1};  // chain's intended size of the run (number of saved points)
    int chainsize{-1};   // chain's current size
    int burnin{-1};      // burnin
    int every{-1};       // subsampling frequency
    int until{-1};       // reading chain until this point
    int currentpoint{-1};
    ProbModel* model;               // the model
    std::string name{"undefined"};  // the name of the chain in the filesystem
};

#endif  // SAMPLE_H
