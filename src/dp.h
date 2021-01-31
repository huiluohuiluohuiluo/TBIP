#ifndef DP_H_
#define DP_H_
#include <iostream>
#include <map>
#include <vector>

#include "graph.h"
using namespace std;

vector<int> dp(Graph graph, int K, int M, double lambda, double theta,
		map<int, double> edge_dis_map, map<int, set<int>> edge_influEdges,
		map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, string input_file,
		string output_file);

vector<int> dpByInfluVal(Graph graph, int K, int M, double lambda, double theta,
		map<int, set<int>> edge_influEdges, map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, vector<Edge> edges,
		string input_file, string input_node_partition, string output_file);

vector<int> dpByMarginInfluVal(Graph graph, int K, int M, double lambda,
		double theta, map<int, set<int>> edge_influEdges,
		map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, vector<Edge> edges,
		string input_file, string input_node_partition, string output_file);

#endif /* DP_H_ */
