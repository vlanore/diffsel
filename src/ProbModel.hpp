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
    ProbModel() {}
    ~ProbModel() {}

    /// obtain the set ("state") of all the nodes of the DAG by a recursive
    /// traversal from the root
    /// nodes to the tips

    // returns the log of the probability (or probability density) mentioned above

    virtual void MakeScheduler() {}

    virtual double Move() {return 1;}
    virtual void Update() {}

    // save model configuration to stream
    virtual void ToStream(std::ostream &os) = 0;
    // get model configuration from stream
    virtual void FromStream(std::istream &is) = 0;

    // monitoring the run
    virtual void Trace(std::ostream & /*unused*/) {}
    virtual void TraceHeader(std::ostream & /*unused*/) {}
    virtual void Monitor(std::ostream &os) {}
    // MCScheduler scheduler;
};

#endif  // PROBMODEL_H
