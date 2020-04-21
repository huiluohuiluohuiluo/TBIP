#include <iostream>
#include <list>
#include <set>

#include "cluster.h"
using namespace std;

Cluster::Cluster(set<int> cEdges, double score){
	this->cEdges = cEdges;
	this->score = score;
}

//Cluster::Cluster(set<int> cEdges, set<int> cInfluEdges, double score){
//	this->cEdges = cEdges;
//	this->cInfluEdges = cInfluEdges;
//	this->score = score;
//}

