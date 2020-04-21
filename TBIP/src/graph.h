#ifndef GRAPH_H_
#define GRAPH_H_
#include <list>
#include <set>

#include "path.h"
#include "edge.h"

using namespace std;

class Graph {
public:
	int V;    		// No. of vertices
	list<int> *adj; // Pointer to an array containing adjacency

	Graph(int V);   // Constructor
	void addEdge(int v, int w);
	set<Path> BFS(int edgeId, double disPara, map<int, int> edge_sourceNode_map,
			map<int, int> edge_desNode_map, map<string, int> nodes_edge_map,
			map<int, double> edge_dis_map);

	map<int, double> getEdgeDensity_afterDiffusion(
			map<int, double> edge_density_map, map<string, double> edge_pro_map,
			map<int, set<Path>> edge_influEdges);
};

#endif /* GRAPH_H_ */
