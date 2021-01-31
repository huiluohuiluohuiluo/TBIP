/*
 * topK.h
 *
 *  Created on: 9 Feb 2020
 *      Author: hui
 */

#ifndef TOPK_H_
#define TOPK_H_
#include <iostream>
#include <map>
#include <vector>

#include "graph.h"

using namespace std;

void topK(int K, double lambda, map<int, set<int>> edge_influEdges,
		map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, string output_file);

void topK_minCut(int K, int M, double lambda,
		map<int, set<int>> edge_influEdges, map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, vector<Edge> edges,
		string input_node_partition, string output_file);

#endif /* TOPK_H_ */
