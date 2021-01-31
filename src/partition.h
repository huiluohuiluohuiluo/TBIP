#ifndef PARTITION_H_
#define PARTITION_H_

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iterator>
#include <map>
#include <set>
#include <bits/stdc++.h>

#include "time.h"
#include "edge.h"
#include "fileProcess.h"

using namespace std;

list<set<int>> getThetaPartition(double theta, double lambda, string input_file,
		map<int, set<int>> edge_influEdges, map<int, double> edge_density_map,
		map<int, double> new_edge_density_map);

list<set<int>> getDensityPartition(int M, vector<Edge> edges,
		string input_node_partition);

cluster_info getDensityPartition(int M, int edgesNum,
		string input_edge_partition);

list<set<int>> getKMeansPartition_list(string input_file,
		map<int, set<int>> old_edge_influEdges);

map<int, int> getEdgePartition(int M, vector<Edge> edges,
		string input_node_partition);

map<int, int> getEdgePartition_list(string input_file,
		map<int, set<int>> old_edge_influEdges);

list<int> getPartitionEdges(int M, vector<Edge> edges,
		string input_node_partition);

double getJaccardDistance(set<int> cEdges1, set<int> cEdges2);

map<int, set<int>> getKMeansPartition(int M, double lambda, string input_file,
		map<int, set<int>> old_edge_influEdges,
		map<int, double> edge_density_map,
		map<int, double> new_edge_density_map);

cluster_info getKMeansPartition(int M, string input_kMeansPartition,
		map<int, set<int>> edgesAndInfluEdges);

cluster_info getPartitionMinCut(int M, string input_partition);

cluster_info getPartitionGrid(int M, string input_partition);

/**
 * given a set of edges, compute their influence scores
 */
double getEdgesInfluScore(double lambda, set<int> edges,
		map<int, set<int>> edge_influEdges, map<int, double> edge_density_map,
		map<int, double> new_edge_density_map);

/**
 * given an edge, compute their influence scores
 */
double getEdgeInfluScore(double lambda, set<int> influEdges,
		map<int, double> edge_density_map,
		map<int, double> new_edge_density_map);

/**
 * given a set of edges, compute their influence edges
 */
set<int> getEdgesInfluEdges(set<int> edges, map<int, set<int>> edge_influEdges);

double getOverlapRatio(double lambda, double theta, set<int> cEdges1,
		set<int> cEdges2, map<int, set<int>> edge_influEdges,
		map<int, double> edge_density_map,
		map<int, double> new_edge_density_map);

list<set<int>> getSubClusters(set<int> edges);

#endif /* PARTITION_H_ */
