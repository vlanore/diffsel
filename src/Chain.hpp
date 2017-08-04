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

#ifndef CHAIN_H
#define CHAIN_H

#include <memory>
#include <string>
#include "ProbModel.hpp"

/// Chain is a Monte Carlo Markov Chain
//  it is responsible for creating a model, applying it to data
//  running a MCMC, to obtain a sample approximately from the posterior
//  distribution
//  saving the sample onto a file, restarting a run
//
// file nomenclature:
// <chainname>.param   : current state
// <chainname>.chain   : list of all points since the beginning of the Monte
// Carlo (burnin included)
// <chainname>.trace   : trace file, each row corresponding to one point of the
// .chain file
// <chainname>.monitor : monitoring the success rate, time spent in each move,
// numerical errors, etc
// <chainname>.run     : put 0 in this file to stop the chain

class Chain {
  public:
    Chain();

    virtual ~Chain() = default;

    virtual void MakeFiles(int force = 0);
    // overwrites files if force == 1

    virtual void Monitor();
    // write the .trace and .monitor files

    virtual void SavePoint();
    // save one point in the .chain file

    virtual void Reset(int force = 0);
    // initialise model and make the files (overwrite if force == 1)

    virtual void Move();
    // perform one cycle of Monte Carlo "moves" (updates)

    virtual void Start();
    // start Monte Carlo

    virtual int GetRunningStatus();
    // returns 0 (means STOP) if one the following conditions holds true
    // 	.run file contains a 0 ("echo 0 > <chainname>.run" is the proper way to
    // stop a chain from a
    // shell)
    // 	size >= until (and until != -1)

    virtual void Run();
    // Move, Monitor amd Save while running status == 1

    virtual std::string GetModelType() = 0;
    // std::string meant as a check when opening files (model name)
    // you should give different names to chains based on different models

    virtual void New(int force = 0) = 0;
    // new chain (force = 1 : this will overwrite files)

    virtual void Open() = 0;
    // open a chain from files

    virtual void Save() = 0;
    // save a chain to files

    std::string GetName() { return name; }

    void SetEvery(int inevery) { every = inevery; }
    void SetUntil(int inuntil) { until = inuntil; }
    void SetName(std::string inname) { name = inname; }

    ProbModel *GetModel() { return model.get(); }
    int GetSize() { return size; }

  protected:
    int every;                         // saving frequency
    int until;                         // intended size of the run (number of saved points)
    int size;                          // current size
    std::unique_ptr<ProbModel> model;  // the model
    std::string name;                  // the name of the chain in the filesystem
    // all files for this chain will be of the form : <name>.<ext>
};

#endif  // CHAIN_H
