#include <iostream>

#include "edge.h"
using namespace std;

Edge::Edge(int edgeId, int sourceNode, int desNode, double dis){
	this->edgeId = edgeId;
	this->sourceNode = sourceNode;
	this->desNode = desNode;
	this->dis = dis;
}
