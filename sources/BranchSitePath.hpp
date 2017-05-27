#ifndef SITEPATH_H
#define SITEPATH_H

#include <string>
#include "StateSpace.hpp"
#include <map>
using namespace std;

#include "SuffStat.hpp"
#include "PoissonSuffStat.hpp"
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

    const Plink *Prev() const;
    const Plink *Next() const;

    bool IsFirst() const;
    bool IsLast() const;

    void Splice();
    void Insert(Plink *link);

    void SetState(int instate);
    int GetState() const;

    void SetRelativeTime(double inrel_time);
    double GetRelativeTime() const;

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

    const Plink *Init() const;
    const Plink *Last() const;

    Plink *Init();
    Plink *Last();
    // StateSpace *GetStateSpace();
    int GetNsub() const;
    /*
    std::string GetState(Plink *link);
    std::string GetCharInitState();
    std::string GetCharFinalState();
    */
    int GetInitState() const;
    int GetFinalState() const;
    // double GetAbsoluteTime(Plink *link);
    double GetRelativeTime(const Plink *link) const { return link->GetRelativeTime(); }

    // void AddCounts(int **paircounts, int *statecounts);

    void AddLengthSuffStat(PoissonSuffStat& suffstat, double factor, const SubMatrix& mat) const;
    void AddPathSuffStat(PathSuffStat& suffstat, double factor) const;

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

inline const Plink *Plink::Prev() const { return prev; }
inline const Plink *Plink::Next() const { return next; }

inline bool Plink::IsFirst() const { return prev == nullptr; }
inline bool Plink::IsLast() const { return next == nullptr; }

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
inline double Plink::GetRelativeTime() const { return rel_time; }
inline int Plink::GetState() const { return state; }

inline Plink *BranchSitePath::Init() { return init; }
inline Plink *BranchSitePath::Last() { return last; }

inline const Plink *BranchSitePath::Init() const { return init; }
inline const Plink *BranchSitePath::Last() const { return last; }

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

inline int BranchSitePath::GetNsub() const { return nsub; }

/*
inline std::string BranchSitePath::GetState(Plink *link) {
    return GetStateSpace()->GetState(link->GetState());
}
inline std::string BranchSitePath::GetCharInitState() { return GetState(init); }
inline std::string BranchSitePath::GetCharFinalState() { return GetState(last); }
*/
inline int BranchSitePath::GetInitState() const { return init->GetState(); }
inline int BranchSitePath::GetFinalState() const { return last->GetState(); }

/*
inline double BranchSitePath::GetAbsoluteTime(Plink *link) {
    return link->GetRelativeTime() * GetTotalTime();
}
*/
#endif  // SITEPATH_H
