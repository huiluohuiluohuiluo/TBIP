#ifndef GREEDY_H_
#define GREEDY_H_
#include <iostream>
#include <map>
#include <vector>

#include "graph.h"

using namespace std;

vector<int> greedy(Graph graph, int K, double lambda,
		map<int, set<int>> edge_influEdges, map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, string output_file);

vector<int> greedyByInfluVal(Graph graph, int K, double lambda,
		map<int, set<int>> edge_influEdges, map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, string output_file);

vector<int> greedyByMarginInfluVal(Graph graph, int K, double lambda,
		map<int, set<int>> edge_influEdges, map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, string output_file);

vector<int> greedyByMarginNumber(Graph graph, int K,
		double lambda, map<int, set<int>> edge_influEdges,
		map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, string output_file);

#endif /* GREEDY_H_ */
