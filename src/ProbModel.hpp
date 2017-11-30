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

#ifndef PROBMODEL_H
#define PROBMODEL_H

#include <iostream>
#include <set>

/// Probabilistic model
/**
 * A model is defined as a Graphical Model:
 * - a directed acyclic graph (DAG) with N nodes
 * - the nodes of the DAG are random variables (X_n)
 * - the joint probability law is given by
 *   \prod_n P(X_n | Pa(X_n))
 *   where  Pa(X_n) represents the set of all parents of node X_n
 *
 * A model implements MCMC:
 * - GetLogProb returns the log of the probability (or probability density)
 * mentioned above
 * - Sample draws a model configuration from this joint probability
 * - Move resample the model's current configuration conditional on the data */

class ProbModel {
  public:
    ProbModel() = default;
    virtual ~ProbModel() = default;

    /// obtain the set ("state") of all the nodes of the DAG by a recursive
    /// traversal from the root
    /// nodes to the tips

    // returns the log of the probability (or probability density) mentioned above

    virtual void MakeScheduler() {}

    virtual double Move() { return 1; }
    virtual void Update() {}

    // save model configuration to stream
    virtual void ToStream(std::ostream &os) = 0;
    // get model configuration from stream
    virtual void FromStream(std::istream &is) = 0;

    // monitoring the run
    virtual void Trace(std::ostream & /*unused*/) {}
    virtual void TraceHeader(std::ostream & /*unused*/) {}
    virtual void Monitor(std::ostream &) {}
    // MCScheduler scheduler;
};

#endif  // PROBMODEL_H
