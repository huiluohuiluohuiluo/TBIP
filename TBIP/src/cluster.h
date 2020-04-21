#ifndef CLUSTER_H_
#define CLUSTER_H_
#include <list>
#include <set>

using namespace std;
class Cluster {
public:
	set<int> cEdges;
//	set<int> cInfluEdges;
	double score;

//	Cluster(set<int> cEdges, set<int> cInfluEdges, double score);
	Cluster(set<int> cEdges, double score);

};

#endif /* CLUSTER_H_ */
