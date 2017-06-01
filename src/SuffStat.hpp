
#ifndef SUFFSTAT_H
#define SUFFSTAT_H

#include <map>
using namespace std;

class SuffStat	{

	public:

	SuffStat() {}
	~SuffStat() {}

	void Clear()	{
		rootcount.clear();
		paircount.clear();
		waitingtime.clear();
	}

	map<int,int> rootcount;
	map<pair<int,int>,int> paircount;
	map<int,double> waitingtime;
};


#endif

