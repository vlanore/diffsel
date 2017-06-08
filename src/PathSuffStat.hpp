
#ifndef PATHSUFFSTAT_H
#define PATHSUFFSTAT_H

#include <map>
using namespace std;

class PathSuffStat {
  public:
    PathSuffStat() {}
    ~PathSuffStat() {}

    void Clear() {
        rootcount.clear();
        paircount.clear();
        waitingtime.clear();
    }

    map<int, int> rootcount;
    map<pair<int, int>, int> paircount;
    map<int, double> waitingtime;
};


#endif
