#ifndef SAMPLE_H
#define SAMPLE_H

#include <fstream>
#include <string>
using namespace std;

class ProbModel;

// Sample is in charge of reading chains from file,
// and computing posterior averages or distributions
class Sample {
  public:
    /// opening chain from file
    Sample(string filename, int in_burnin = 0, int in_every = 1, int in_until = -1);

    virtual ~Sample() = default;

    // return the base for all file names (same as for the chain object)
    string GetName() { return name; }

    // get the next point (automatically accounts for subsampling, as specified by the "every"
    // parameter)
    // this point is accessed to using GetModel()
    void GetNextPoint();

    // returns a pointer to current point
    // in general this function will be overriden by derived classes
    // because derived classes manipulate models that are derived from ProbModel
    virtual ProbModel* GetModel() { return model; }

    // string meant as a check when opening files (model name)
    // you should give different names to chains based on different models
    virtual string GetModelType() = 0;

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
    ifstream chain_is;
    int chainevery{-1};  // chain's saving frequency
    int chainuntil{-1};  // chain's intended size of the run (number of saved points)
    int chainsize{-1};   // chain's current size
    int burnin{-1};      // burnin
    int every{-1};       // subsampling frequency
    int until{-1};       // reading chain until this point
    int currentpoint{-1};
    ProbModel* model;          // the model
    string name{"undefined"};  // the name of the chain in the filesystem
};

#endif  // SAMPLE_H
