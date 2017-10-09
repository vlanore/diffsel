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
#include "StateSpace.hpp"
using namespace std;

#include "PathSuffStat.hpp"
#include "SubMatrix.hpp"


class Plink {
    friend class BranchSitePath;

  public:
    Plink();
    Plink(int instate, double inrel_time);
    ~Plink();
    Plink *Prev();
    Plink *Next();

    bool IsFirst();
    bool IsLast();

    void Splice();
    void Insert(Plink *link);

    void SetState(int instate);
    int GetState();

    void SetRelativeTime(double inrel_time);
    double GetRelativeTime();

  private:
    Plink *next;
    Plink *prev;

    int state;
    double rel_time;
};

class BranchSitePath {
  public:
    BranchSitePath();
    BranchSitePath(int state);
    virtual ~BranchSitePath();

    /*
    virtual double GetTotalTime() = 0;
    virtual void SetTotalTime(double intime) = 0;
    */

    Plink *Init();
    Plink *Last();
    // StateSpace *GetStateSpace();
    int GetNsub();
    /*
    std::string GetState(Plink *link);
    std::string GetCharInitState();
    std::string GetCharFinalState();
    */
    int GetInitState();
    int GetFinalState();
    // double GetAbsoluteTime(Plink *link);
    double GetRelativeTime(Plink *link) { return link->GetRelativeTime(); }

    void AddCounts(int **paircounts, int *statecounts);

    void AddLengthSuffStat(int &count, double &beta, double factor, SubMatrix *mat);
    void AddSuffStat(PathSuffStat &suffstat, double factor);

    /*
    void SetTimesRelativeToAbsolute();
    void SetTimesAbsoluteToRelative();
    double CheckTotalTime();
    */

    void Reset(int state);
    void Append(int instate, double reltimelength);

    void BKReset(int state);
    void BKAppend(int instate, double reltimelength);
    void BackupPath();
    void RestorePath();

    // pulley: the beginning of path p, up to absolute time point t,
    // is inverted and transferred at the base of this path
    // path p is modified accordingly
    // void Prefix(BranchSitePath *p, BranchSitePath *root, double abstime);

    // std::string ToString(bool redundant = false);

  protected:
    Plink *init;
    Plink *last;
    int nsub;

    Plink *bkinit;
    Plink *bklast;
    int bknsub;
};

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* Inline definitions
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
//	* Plink
//-------------------------------------------------------------------------

inline Plink::Plink() : next(nullptr), prev(nullptr), state(0), rel_time(0) {}
inline Plink::Plink(int instate, double inrel_time)
    : next(nullptr), prev(nullptr), state(instate), rel_time(inrel_time) {}
inline Plink::~Plink() { Splice(); }

inline Plink *Plink::Prev() { return prev; }
inline Plink *Plink::Next() { return next; }

inline bool Plink::IsFirst() { return prev == nullptr; }
inline bool Plink::IsLast() { return next == nullptr; }

inline void Plink::Insert(Plink *link) {
    link->next = next;
    if (next != nullptr) {
        next->prev = link;
    }
    link->prev = this;
    next = link;
}

inline void Plink::Splice() {
    if (prev != nullptr) {
        prev->next = next;
    }
    if (next != nullptr) {
        next->prev = prev;
    }
    prev = next = nullptr;
}

inline void Plink::SetState(int instate) { state = instate; }
inline void Plink::SetRelativeTime(double inrel_time) { rel_time = inrel_time; }
inline double Plink::GetRelativeTime() { return rel_time; }
inline int Plink::GetState() { return state; }

inline Plink *BranchSitePath::Init() { return init; }
inline Plink *BranchSitePath::Last() { return last; }

//-------------------------------------------------------------------------
//	* BranchSitePath
//-------------------------------------------------------------------------

inline void BranchSitePath::Append(int instate, double reltimelength) {
    last->SetRelativeTime(reltimelength);
    auto link = new Plink(instate, 0);
    last->Insert(link);
    last = link;
    nsub++;
}

inline void BranchSitePath::BKAppend(int instate, double reltimelength) {
    bklast->SetRelativeTime(reltimelength);
    auto link = new Plink(instate, 0);
    bklast->Insert(link);
    bklast = link;
    bknsub++;
}

inline int BranchSitePath::GetNsub() { return nsub; }

/*
inline std::string BranchSitePath::GetState(Plink *link) {
    return GetStateSpace()->GetState(link->GetState());
}
inline std::string BranchSitePath::GetCharInitState() { return GetState(init); }
inline std::string BranchSitePath::GetCharFinalState() { return GetState(last); }
*/
inline int BranchSitePath::GetInitState() { return init->GetState(); }
inline int BranchSitePath::GetFinalState() { return last->GetState(); }

/*
inline double BranchSitePath::GetAbsoluteTime(Plink *link) {
    return link->GetRelativeTime() * GetTotalTime();
}
*/
#endif  // SITEPATH_H
