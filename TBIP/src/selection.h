#ifndef SELECTION_H_
#define SELECTION_H_
#include <iostream>
#include <map>
#include <vector>

#include "graph.h"
using namespace std;

vector<int> selection(Graph graph, int K, int M, double disPara, double lambda,
		double theta, map<int, set<int>> edge_influEdges,
		map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, vector<Edge> edges,
		string input_file, string input_node_partition, string output_file);

#endif /* SELECTION_H_ */
