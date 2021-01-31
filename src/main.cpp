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
#include "globalSel.h"
#include "greedySel.h"
#include "branchAndbound.h"
#include "duration.h"

using namespace std;

//basic parameters:
double eplison = 0.01;
double lambda = 111000; // assume there is one vehicle per meter, lambda is to judge whether an edge is a traffic bottleneck, which will used to multiply with the road length
string city = "xian"; 	// which city
int hour = 11;      	// which hour
int spreadWindow = 20; // indicate the time period for the traffic spread window, and the information will be updated every spread window, the unit is second
int monitorWindow = 3600; // indicate the time period for the monitor window, and select K seed edges every monitor window, the unit is second
int batchSize = monitorWindow / spreadWindow;
int durationSize = 15; // durationSize is to define the time interval a bottleneck may last (duration time = durationSize * spreadWindow
int algFlag = 5;
//0: greedy; 1: dp; 2: selection; 3: partition 4: greedy+sampling 5: topMinCut
int pruneFlag = 0;
int K = 15;

string baseFile = "dataset/"; // path for server
//string baseFile = "dataset/";

//optional parameters for different methods:
//double theta = 0.05;		// partition ratio
int M = 100;	// the number of partitions
string pType = "max";
string parType = "kmeans";
//string parType = "grid";
//string parType = "mincut";

int main(int argc, char* argv[]) {
//	int dataId = atoi(argv[1]);
//	algFlag = atoi(argv[2]);	//0: greedy; 1: dp; 2: selection; 3: partition
//	eplison = atoi(argv[3]) / 1000.0;
//	lambda = atoi(argv[4]) * 111000;
//	hour = atoi(argv[5]);
//	spreadWindow = atoi(argv[6]);
//	monitorWindow = atoi(argv[7]);
//	durationSize = atoi(argv[8]);
//	pruneFlag = atoi(argv[9]);
//	M = atoi(argv[10]);
//	int parTypeId = atoi(argv[11]);
//
//	batchSize = monitorWindow / spreadWindow;
//
//	if (dataId == 0) {
//		city = "chengdu";
//	} else if (dataId == 1) {
//		city = "xian";
//	} else if (dataId == 2) {
//		city = "porto";
//	}
//
//	if (parTypeId == 0) {
//		parType = "kmeans";
//	} else if (parTypeId == 1) {
//		parType = "grid";
//	}
//	comment above

	string parameterStr = to_string(lambda) + "_" + to_string(hour) + "_"
			+ to_string(spreadWindow) + "_" + to_string(monitorWindow) + "_"
			+ to_string(durationSize) + "_" + to_string(M) + "_" + pType + "_"
			+ parType;
	string input_edge_basic = baseFile + city + "/road.txt";
	string input_edge_density = baseFile + city + "/edge_density/"
			+ to_string(hour) + ".txt";
	string input_edge_pro = baseFile + city + "/edge_pro/" + to_string(hour)
			+ ".txt";
	string input_edge_residual = baseFile + city + "/edge_residual/"
			+ to_string(hour) + ".txt";
	string input_edge_influPaths = baseFile + city + "/edge_influ/"
			+ to_string(spreadWindow) + "_" + to_string(hour) + "_paths.txt";
	string input_edge_influEdges = baseFile + city + "/edge_influ/"
			+ to_string(spreadWindow) + "_" + to_string(hour) + "_edges.txt";
	string input_edge_speed = baseFile + city + "/edge_speed/" + to_string(hour)
			+ ".txt";

	string input_partition_file;
	string input_partition_minCut_file;
	if (parType == "kmeans") {
		input_partition_file = baseFile + city + "/partition/" + parameterStr
				+ ".txt";
	} else {
		if (city == "chengdu") {
			input_partition_file = baseFile + city
					+ "/partition/road_partition_6135_" + to_string(M) + ".txt";
		} else if (city == "xian") {
			input_partition_file = baseFile + city
					+ "/partition/road_partition_5045_" + to_string(M) + ".txt";
		} else if (city == "porto") {
			input_partition_file = baseFile + city
					+ "/partition/road_partition_108384_" + to_string(M)
					+ ".txt";
		}
	}

	if (city == "chengdu") {
		input_partition_minCut_file = baseFile + city + "/partition/" + city
				+ "_" + to_string(hour) + ".txt.part.100";
	} else if (city == "xian") {
		input_partition_minCut_file = baseFile + city + "/partition/" + city
				+ "_" + to_string(hour) + ".txt.part.100";
	} else if (city == "porto") {
		input_partition_minCut_file = baseFile + city + "/partition/" + city
				+ "_" + to_string(hour) + ".txt.part.2000";
	}
	string input_influEdges_file = baseFile + city + "/edge_influ/realtime_"
			+ to_string(lambda) + "_" + to_string(hour) + "_"
			+ to_string(spreadWindow) + "_" + to_string(monitorWindow) + "_"
			+ to_string(durationSize);
	string output_file = baseFile + city + "/output.txt";

	map<int, Edge> edges;
	map<int, double> edge_dis_map;
	map<int, int> edge_sourceNode_map;
	map<int, int> edge_desNode_map;
	map<string, int> nodes_edge_map;

	map<int, double> edge_density_map;
	map<string, double> edge_pro_map;
	map<int, double> edge_residual_map;
	map<int, double> edge_speed_map;
	map<int, double> edge_time_map;
	int nodesNum;
	int edgesNum;

//	read files as input
	if (city == "chengdu") {
		nodesNum = 4326;
		edgesNum = 6135;
	} else if (city == "xian") {
		nodesNum = 2086;
		edgesNum = 5045;
	} else if (city == "porto") { // #traj: 1565595 # time range: -70 to 207 (278 totally)
		nodesNum = 60287;
		edgesNum = 108571;
	}
	edges = readEdgeBasicFile(edgesNum, input_edge_basic);
	edge_density_map = readEdgeDensityFile(input_edge_density);
	edge_speed_map = readEdgeSpeedFile(input_edge_speed);

	map<int, Edge>::iterator iter_edge;

	for (iter_edge = edges.begin(); iter_edge != edges.end(); iter_edge++) {
		int i = iter_edge->first;
		Edge edge = edges[i];
		if (edge.edgeId != -1) {
			edge_dis_map[edge.edgeId] = edge.dis;
			edge_sourceNode_map[edge.edgeId] = edge.sourceNode;
			edge_desNode_map[edge.edgeId] = edge.desNode;
			nodes_edge_map[to_string(edge.sourceNode) + ","
					+ to_string(edge.desNode)] = edge.edgeId;
		}
		if (edge_density_map.find(i) == edge_density_map.end()) {
			edge_density_map[i] = 0;
		}
		if (edge_speed_map.find(i) == edge_speed_map.end()
				|| edge_speed_map[i] == 0) {
			edge_speed_map[i] = 30; // if there is not speed information, then we set the default speed as 30km/h
		}
		edge_time_map[i] = edge_dis_map[i] / edge_speed_map[i] * 3600 * 111;

	}

//	K = eplison * edgesNum;
//	contruct the graph and insert edges
	Graph graph(nodesNum);
	for (int i = 0; i < edgesNum; i++) {
		Edge edge = edges[i];
		if (edge.edgeId != -1 && edge.sourceNode != -1 && edge.desNode != -1) {
			graph.addEdge(edge.sourceNode, edge.desNode);
		}
	}

	edge_pro_map = readEdgeProFile(input_edge_pro);
	edge_residual_map = readEdgeResidualFile(edgesNum, input_edge_residual);

//	graph.getRoadforPartitionWithoutWeight(edges, output_road_for_partition);
	vector<int> seed_edges;
	clock_t time1, time2;

//	find the influenced edges and paths for each edge
	map<int, set<int>> edge_influEdges;
	map<int, set<int>> r_edge_influEdges; // for traffic spread

	map<int, map<int, double>> new_edge_pro_map;
	new_edge_pro_map = getEdgeInfluPaths(input_edge_influPaths,
			input_edge_influEdges, input_edge_residual, graph, spreadWindow,
			edge_sourceNode_map, edge_desNode_map, edge_pro_map, nodes_edge_map,
			edge_dis_map, edge_residual_map, edge_time_map);

	edge_influEdges = getEdgeInfluEdges(input_edge_influEdges, edgesNum);

	r_edge_influEdges = getReverseEdgeInfluEdges(edgesNum,
			input_edge_influEdges);

	cout << "haha: " << edge_density_map.size() << " "
			<< edge_residual_map.size() << " " << edge_influEdges.size()
			<< endl;

	map<int, double> new_edge_density_map =
			graph.getEdgeDensity_afterDiffusion_Pro(edge_density_map,
					new_edge_pro_map, edge_residual_map, r_edge_influEdges);

	time1 = clock();
	durationSel(graph, K, M, batchSize, durationSize, lambda, edge_dis_map,
			edge_influEdges, r_edge_influEdges, edge_density_map,
			new_edge_density_map, new_edge_pro_map, edge_residual_map,
			input_partition_file, output_file, parType, pruneFlag, algFlag,
			city, input_partition_minCut_file, edges, input_influEdges_file);

	time2 = clock();
	double duration22 = (double) (time2 - time1) / CLOCKS_PER_SEC;

	string parameters2 = to_string(K) + " " + to_string(eplison) + " "
			+ to_string(batchSize) + " " + to_string(spreadWindow) + " "
			+ to_string(lambda) + " " + city + " " + to_string(hour) + " "
			+ to_string(durationSize) + " " + to_string(M) + " " + pType + " "
			+ parType;
	writeOutputFile(output_file, parameters2);

	string runtime_o2 = "Running time: " + to_string(duration22) + " seconds"
			+ "\n" + "***********";
	writeOutputFile(output_file, runtime_o2);
	cout << runtime_o2 << endl;
	return 0;
}
