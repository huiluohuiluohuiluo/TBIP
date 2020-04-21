#include <iostream>
#include <string>
#include <map>
#include <list>
#include <vector>
#include <set>
#include <algorithm>

#include "time.h"
#include "graph.h"
#include "fileProcess.h"
#include "greedy.h"
#include "topK.h"
#include "dp.h"
#include "selection.h"
using namespace std;

int K = 50;
double disPara = 0.005;
//double disPara = -1; // disPara is to control the spatial influence range centred at a specific edge: 194m->0.0087(5-hop)
double lambda = 450; // lambda is to judge whether an edge is a traffic bottleneck
string city = "chengdu"; 	// which city
int hour = 14;      	// which hour
double theta = 0.05;		// partition ratio
int M = 10;	// the number of partitions
string pType = "max";
string parType = "Density";
//string baseFile = "../dataset/"; // path for server
string baseFile = "dataset/";

//0: top-k from all edges
//1: top-k from minimum cuts
//2: greedy
//3: greedyByInfluVal
//4: greedyByMarginInfluVal
//5: dynamic programming
//6: dpByInfluVal
//7: dpByMarginInfluVal
//8: selection
int algorithm = 44;

int main() {
	string input_edge_basic = baseFile + city + "/road.txt";
	string input_edge_density = baseFile + city + "/edge_density/"
			+ to_string(hour) + ".txt";
	string input_edge_pro = baseFile + city + "/edge_pro/" + to_string(hour)
			+ ".txt";
	string input_edge_influPaths = baseFile + city + "/edge_influ/"
			+ to_string(disPara) + "_paths.txt";
	string input_edge_influEdges = baseFile + city + "/edge_influ/"
			+ to_string(disPara) + "_edges.txt";
	string input_partition = baseFile + city + "/partition/" + to_string(theta)
			+ "_" + to_string(lambda) + "_" + to_string(hour) + ".txt";
	string input_kMeansPartition = baseFile + city + "/partition/kMeans_"
			+ to_string(theta) + "_" + to_string(lambda) + "_" + to_string(hour)
			+ "_" + pType + ".txt";
	string input_densityPartition = baseFile + city + "/partition/" + parType
			+ "_" + to_string(theta) + "_" + to_string(lambda) + "_"
			+ to_string(hour) + ".txt";
	string input_node_partition = baseFile + city + "/partition/"
			+ to_string(hour) + ".txt.part." + to_string(M);
	string input_edge_partition = input_densityPartition + ".part."
			+ to_string(M);
	string input_edges_partition = input_edge_partition + "_Edges";
	string output_file = baseFile + city + "/output.txt";

	vector<Edge> edges;
	map<int, double> edge_dis_map;
	map<int, int> edge_sourceNode_map;
	map<int, int> edge_desNode_map;
	map<string, int> nodes_edge_map;

	map<int, double> edge_density_map;
	map<string, double> edge_pro_map;

//	read files as input
	edges = readEdgeBasicFile(input_edge_basic);
	int edgesNum = edges.size();

	for (int i = 0; i < edgesNum; i++) {
		Edge edge = edges[i];
		edge_dis_map[edge.edgeId] = edge.dis;
		edge_sourceNode_map[edge.edgeId] = edge.sourceNode;
		edge_desNode_map[edge.edgeId] = edge.desNode;
		nodes_edge_map[to_string(edge.sourceNode) + ","
				+ to_string(edge.desNode)] = edge.edgeId;
	}
	int nodesNum;
	if (city == "chengdu") {
		nodesNum = 4326;
	} else if (city == "xian") {
		nodesNum = 2086;
	}

	edge_density_map = readEdgeDensityFile(input_edge_density);
	edge_pro_map = readEdgeProFile(input_edge_pro);

//	contruct the graph and insert edges
	Graph graph(nodesNum);
	for (int i = 0; i < edgesNum; i++) {
		Edge edge = edges[i];
		graph.addEdge(edge.sourceNode, edge.desNode);
	}

//	time counter
	vector<int> seed_edges;
	clock_t time, time0, time1, time2;
	time = clock();

//	find the influenced edges and paths for each edge
	map<int, set<Path>> edge_influPaths;
	map<int, set<int>> edge_influEdges;

	edge_influPaths = getEdgeInfluPaths(input_edge_influPaths,
			input_edge_influEdges, graph, edgesNum, disPara,
			edge_sourceNode_map, edge_desNode_map, nodes_edge_map, edge_dis_map,
			edge_density_map);
	edge_influEdges = getEdgeInfluEdges(input_edge_influEdges);
	cout << "finish" << endl;

	time0 = clock();

	map<int, double> new_edge_density_map = graph.getEdgeDensity_afterDiffusion(
			edge_density_map, edge_pro_map, edge_influPaths);
	time1 = clock();

	if (algorithm == 0) {
		topK(K, lambda, edge_influEdges, edge_density_map, new_edge_density_map,
				output_file);
	} else if (algorithm == 1) {
		topK_minCut(K, M, lambda, edge_influEdges, edge_density_map,
				new_edge_density_map, edges, input_node_partition, output_file);
	} else if (algorithm == 2) {
		seed_edges = greedy(graph, K, lambda, edge_influEdges, edge_density_map,
				new_edge_density_map, output_file);
	} else if (algorithm == 3) {
		seed_edges = greedyByInfluVal(graph, K, lambda, edge_influEdges,
				edge_density_map, new_edge_density_map, output_file);
	} else if (algorithm == 4) {
		seed_edges = greedyByMarginInfluVal(graph, K, lambda, edge_influEdges,
				edge_density_map, new_edge_density_map, output_file);
	} else if (algorithm == 44) {
		seed_edges = greedyByMarginNumber(graph, K, lambda, edge_influEdges,
				edge_density_map, new_edge_density_map, output_file);
	} else if (algorithm == 5) {
		seed_edges = dp(graph, K, M, lambda, theta, edge_influEdges,
				edge_density_map, new_edge_density_map, edges,
				input_edges_partition, input_node_partition, output_file);
	} else if (algorithm == 6) {
		seed_edges = dpByInfluVal(graph, K, M, lambda, theta, edge_influEdges,
				edge_density_map, new_edge_density_map, edges,
				input_edges_partition, input_node_partition, output_file);
	} else if (algorithm == 7) {
		seed_edges = dpByMarginInfluVal(graph, K, M, lambda, theta,
				edge_influEdges, edge_density_map, new_edge_density_map, edges,
				input_edges_partition, input_node_partition, output_file);
	} else if (algorithm == 8) {
		seed_edges = selection(graph, K, M, disPara, lambda, theta,
				edge_influEdges, edge_density_map, new_edge_density_map, edges,
				input_edges_partition, input_node_partition, output_file);
	} else if (algorithm == -1) {
//		getRoadforPartitionOnDensity(edges, edge_density_map, edge_pro_map,
//				input_densityPartition);
		outputEdgesPartition(M, input_edge_partition, input_edges_partition);

	}

	time2 = clock();
	double duration0 = (double) (time0 - time) / CLOCKS_PER_SEC;
	double duration1 = (double) (time1 - time0) / CLOCKS_PER_SEC;
	double duration2 = (double) (time2 - time1) / CLOCKS_PER_SEC;

	string parameters = to_string(K) + " " + to_string(M) + " "
			+ to_string(disPara) + " " + to_string(lambda) + " " + city + " "
			+ to_string(hour);
	writeOutputFile(output_file, parameters);

	string runtime_o = "Running time: " + to_string(duration0) + " "
			+ to_string(duration1) + " " + to_string(duration2) + " seconds"
			+ "\n" + "***********";
	writeOutputFile(output_file, runtime_o);
	cout << runtime_o << endl;

	return 0;
}
