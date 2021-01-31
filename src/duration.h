/*
 * duration.h
 *
 *  Created on: 20 Aug 2020
 *      Author: hui
 */

#ifndef DURATION_H_
#define DURATION_H_

#include <iostream>
#include <map>
#include <vector>

#include "graph.h"

using namespace std;

set<int> *durationSel(Graph graph, int K, int M, int batchSize,
		int durationSize, double lambda, map<int, double> edge_dis_map,
		map<int, set<int>> edge_influEdges,
		map<int, set<int>> r_edge_influEdges, map<int, double> edge_density_map,
		map<int, double> new_edge_density_map,
		map<int, map<int, double>> edge_pro_map,
		map<int, double> edge_residual_map, string input_partition_file,
		string output_file, string parType, int pruneFlag, int algFlag,
		string city, string input_partition_minCut_file,
		map<int, Edge> edges_map, string input_influEdges_file);

set<int> *greedy(int K, int edgesNum, map<int, set<int>> edgesAndInfluEdges);

set<int> *greedySelSampling(int K, int edgesNum,
		map<int, set<int>> edgesAndInfluEdges, string city);

set<int> *dp(int K, int edgesNum, int M, map<int, set<int>> edgesAndInfluEdges,
		string inputFile_partition, string parType);

set<int> *selection(int K, int edgesNum, int M,
		map<int, set<int>> edgesAndInfluEdges, string inputFile_partition,
		string parType);

set<int> *partition(int K, int edgesNum, int M,
		map<int, set<int>> edgesAndInfluEdges, string input_partition_file,
		string parType);

set<int> *greedySelForPartition(int k, set<int> edges, int edgesNum,
		map<int, set<int>> edge_influEdges, set<int> s_hat, set<int> i_hat);

set<int> *topK_minCut(int K, int edgesNum, string city,
		map<int, set<int>> edgesAndInfluEdges, string input_partition_file,
		map<int, Edge> edges_map, map<int, double> edge_density_map);

set<int> getMinCutsEdges(int edgesNum, map<int, Edge> edges,
		string input_node_partition);

#endif /* DURATION_H_ */
