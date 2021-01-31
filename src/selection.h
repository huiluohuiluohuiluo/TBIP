#ifndef SELECTION_H_
#define SELECTION_H_
#include <iostream>
#include <map>
#include <vector>

#include "graph.h"
using namespace std;

vector<int> selection(Graph graph, int K, int M, double lambda,
		double theta, map<int, double> edge_dis_map,
		map<int, set<int>> edge_influEdges, map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, string input_file,
		string output_file);

#endif /* SELECTION_H_ */
