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

#ifndef SITEPATH_H
#define SITEPATH_H

#include <map>
#include <string>
#include "PathSuffStat.hpp"
#include "StateSpace.hpp"
#include "SubMatrix.hpp"


class Plink {
    friend class BranchSitePath;

  public:
    Plink() : next(nullptr), prev(nullptr), state(0), rel_time(0) {}
    Plink(int instate, double inrel_time)
        : next(nullptr), prev(nullptr), state(instate), rel_time(inrel_time) {}
    ~Plink() { Splice(); }

    Plink *Prev() { return prev; }
    Plink *Next() { return next; }

    bool IsFirst() { return prev == nullptr; }
    bool IsLast() { return next == nullptr; }

    void Insert(Plink *link) {
        link->next = next;
        if (next != nullptr) {
            next->prev = link;
        }
        link->prev = this;
        next = link;
    }

    void Splice() {
        if (prev != nullptr) {
            prev->next = next;
        }
        if (next != nullptr) {
            next->prev = prev;
        }
        prev = next = nullptr;
    }

    void SetState(int instate) { state = instate; }
    void SetRelativeTime(double inrel_time) { rel_time = inrel_time; }
    double GetRelativeTime() { return rel_time; }
    int GetState() { return state; }

  private:
    Plink *next;
    Plink *prev;

    int state;
    double rel_time;
};

class BranchSitePath {
  public:
    BranchSitePath();
    explicit BranchSitePath(int state);
    virtual ~BranchSitePath();

    Plink *Init() { return init; }
    Plink *Last() { return last; }

    int GetNsub() { return nsub; }
    int GetInitState() { return init->GetState(); }
    int GetFinalState() { return last->GetState(); }
    double GetRelativeTime(Plink *link) { return link->GetRelativeTime(); }

    void AddCounts(int **paircounts, int *statecounts);
    void AddLengthSuffStat(int &count, double &beta, double factor, SubMatrix *mat);
    void AddSuffStat(PathSuffStat &suffstat, double factor);

    void Reset(int state);
    void Append(int instate, double reltimelength) {
        last->SetRelativeTime(reltimelength);
        auto link = new Plink(instate, 0);
        last->Insert(link);
        last = link;
        nsub++;
    }

    void BKReset(int state);
    void BKAppend(int instate, double reltimelength) {
        bklast->SetRelativeTime(reltimelength);
        auto link = new Plink(instate, 0);
        bklast->Insert(link);
        bklast = link;
        bknsub++;
    }
    void BackupPath();
    void RestorePath();

    // pulley: the beginning of path p, up to absolute time point t,
    // is inverted and transferred at the base of this path
    // path p is modified accordingly
    // void Prefix(BranchSitePath *p, BranchSitePath *root, double abstime);

  protected:
    Plink *init;
    Plink *last;
    int nsub;

    Plink *bkinit;
    Plink *bklast;
    int bknsub;
};

#endif  // SITEPATH_H
