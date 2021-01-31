#ifndef FILEPROCESS_H_
#define FILEPROCESS_H_
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iterator>
#include <set>

#include "edge.h"
#include "path.h"
#include "graph.h"
using namespace std;

struct cluster_info {
	map<int, set<int>> mClusters;
	map<int, int> edgeAndPartition;
};

void split(string& s, vector<string>& sv, char flag);
map<int, Edge> readEdgeBasicFile(int edgesNum, string input_file);
set<int> readEdgeSmallFile(string input_file);
map<int, double> readEdgeDensityFile(string input_file);
map<int, double> readEdgeSpeedFile(string input_file);
map<string, double> readEdgeProFile(string input_file);
map<int, map<int, double>> readEdgeNewProFile(string input_file);
map<int, double> readEdgeResidualFile(int edgesNum, string input_file);
map<string, double> getEdgeInfluPaths(string input_file_paths,
		string input_file_edges, Graph graph, int edgesNum, double disPara,
		map<int, int> edge_sourceNode_map, map<int, int> edge_desNode_map,
		map<string, double> edge_pro_map, map<string, int> nodes_edge_map,
		map<int, double> edge_dis_map, map<int, double> edge_density_map);
map<int, map<int, double>> getEdgeInfluPaths(string input_edge_influPaths,
		string input_edge_influEdges, string input_edge_residual_new,
		Graph graph, int spreadWindow, map<int, int> edge_sourceNode_map,
		map<int, int> edge_desNode_map, map<string, double> edge_pro_map,
		map<string, int> nodes_edge_map, map<int, double> edge_dis_map,
		map<int, double> edge_residual_map, map<int, double> edge_time_map);
map<int, set<int>> getEdgeInfluEdges(string input_file, int edgesNum);
map<int, set<int>> getReverseEdgeInfluEdges(int edgesNum, string input_file);
void writeOutputFile(string output_file, string line);
map<int, int> getNodePartition(string input_node_partition);
map<int, int> getEdgePartition(string input_edge_partition);
void getRoadforPartitionOnDensity(vector<Edge> edges,
		map<int, double> edge_density_map, map<string, double> edge_pro_map,
		string output_road_for_partition);
void outputEdgesPartition(int partitionNum, string input_edge_partition,
		string output_partition);

#endif /* FILEPROCESS_H_ */
