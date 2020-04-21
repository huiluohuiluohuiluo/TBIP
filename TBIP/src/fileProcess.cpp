#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iterator>
#include <map>
#include <set>

#include "edge.h"
#include "fileProcess.h"
#include "path.h"
#include "graph.h"
using namespace std;

void split(string& s, vector<string>& sv, char flag) {
	sv.clear();
	istringstream iss(s);
	string temp;

	while (getline(iss, temp, flag)) {
		sv.push_back(temp);
//		stoi(str) -> int
//		atoi(str.c_str()) -> str
	}
}

vector<Edge> readEdgeBasicFile(string input_file) {
	string line;
	vector<Edge> edges;
	ifstream file(input_file);

	if (!file) {
		cout << "The file can not open!" << endl;
	} else {
		while (getline(file, line)) {
			vector<string> line_vec;
			split(line, line_vec, ',');
			int edgeId = stoi(line_vec[0]);
			int sourceNode = stoi(line_vec[1]);
			int desNode = stoi(line_vec[2]);
			istringstream iss(line_vec[3]);
			double dis;
			iss >> dis;
//			cout << edgeId << " " << sourceNode << " " << desNode << " " << dis << endl;
			Edge e(edgeId, sourceNode, desNode, dis);
			edges.push_back(e);
		}
	}
//	cout << edges.size() << endl;
	return edges;
}

map<int, double> readEdgeDensityFile(string input_file) {
	string line;
	ifstream file(input_file);
	map<int, double> edge_density_map;

	if (!file) {
		cout << "The file can not open!" << endl;
	} else {
		while (getline(file, line)) {
			vector<string> line_vec;
			split(line, line_vec, ',');
			int edgeId = stoi(line_vec[0]);
			istringstream iss(line_vec[1]);
			double density;
			iss >> density;
			edge_density_map[edgeId] = density;
		}
	}
	return edge_density_map;
}

map<string, double> readEdgeProFile(string input_file) {
	string line;
	ifstream file(input_file);
	map<string, double> edge_pro_map;

	if (!file) {
		cout << "The file can not open!" << endl;
	} else {
		while (getline(file, line)) {
			vector<string> line_vec;
			split(line, line_vec, ';');
//			int source = stoi(line_vec[0]);
//			int des = stoi(line_vec[1]);
			vector<string> line_vec_1;
			split(line_vec[1], line_vec_1, ',');
			string key = line_vec[0] + ',' + line_vec_1[0];
			istringstream iss(line_vec_1[1]);
			double pr;
			iss >> pr;
			edge_pro_map[key] = pr;
		}
	}
	return edge_pro_map;
}

/**
 * input_edge_influPaths: the influenced path for each edge
 * input_edge_influEdges: the influenced edges for each edge
 */
map<int, set<Path>> getEdgeInfluPaths(string input_file_paths,
		string input_file_edges, Graph graph, int edgesNum, double disPara,
		map<int, int> edge_sourceNode_map, map<int, int> edge_desNode_map,
		map<string, int> nodes_edge_map, map<int, double> edge_dis_map,
		map<int, double> edge_density_map) {
	string line;
	ifstream file(input_file_paths);
	map<int, set<Path>> edge_influPaths;
	map<int, set<int>> edge_influEdges;
//	write into files
	ofstream outfile_paths(input_file_paths, ios::app);
	ofstream outfile_edges(input_file_edges, ios::app);

	if (!file) {
		cout << "Create edge influence paths file ..." << endl;
		for (int i = 0; i < edgesNum; i++) {
			set<Path> influPaths;
			set<int> influEdges;
			if (edge_density_map[i] != 0) {
				influPaths = graph.BFS(i, disPara, edge_sourceNode_map,
						edge_desNode_map, nodes_edge_map, edge_dis_map);
				set<Path>::iterator iter_set;
				for (iter_set = influPaths.begin();
						iter_set != influPaths.end(); iter_set++) {
					Path path = *iter_set;
					list<int>::iterator iter_list;
					for (iter_list = path.edges.begin();
							iter_list != path.edges.end(); iter_list++) {
						influEdges.insert(*iter_list);
					}
				}
			}
			cout << "edgeId is: " << i << endl;
			edge_influPaths[i] = influPaths;
			set<Path>::iterator iter_set_paths;
			string line_paths = to_string(i) + ";";
			for (iter_set_paths = influPaths.begin();
					iter_set_paths != influPaths.end(); iter_set_paths++) {
				Path path = *iter_set_paths;
				list<int>::iterator iter_list;
				line_paths += to_string(path.dis) + "-";
				for (iter_list = path.edges.begin();
						iter_list != path.edges.end(); iter_list++) {
					line_paths += to_string(*iter_list) + ",";
				}
				line_paths.pop_back();
				line_paths += "*";
			}
			line_paths.pop_back();
			outfile_paths << line_paths << endl;

			edge_influEdges[i] = influEdges;
			set<int>::iterator iter_set_edges;
			string line_edges = to_string(i) + ";";
			for (iter_set_edges = influEdges.begin();
					iter_set_edges != influEdges.end(); iter_set_edges++) {
				if (i != *iter_set_edges) {
					line_edges += to_string(*iter_set_edges) + ",";
				}
			}
			line_edges.pop_back();
			outfile_edges << line_edges << endl;
		}
		cout << "finish ..." << endl;

//		for (int i = 0; i < edgesNum; i++) {
//			set<Path> influPaths = edge_influPaths[i];
//			set<Path>::iterator iter_set;
//			string line = to_string(i) + ";";
//			for (iter_set = influPaths.begin(); iter_set != influPaths.end();
//					iter_set++) {
//				Path path = *iter_set;
//				list<int>::iterator iter_list;
//				line += to_string(path.dis) + "-";
//				for (iter_list = path.edges.begin();
//						iter_list != path.edges.end(); iter_list++) {
//					line += to_string(*iter_list) + ",";
//				}
//				line.pop_back();
//				line += "*";
//			}
//			line.pop_back();
//			outfile_paths << line << endl;
//		}

//		for (int i = 0; i < edgesNum; i++) {
//			set<int> influEdges = edge_influEdges[i];
//			set<int>::iterator iter_set;
//			string line = to_string(i) + ";";
//			for (iter_set = influEdges.begin(); iter_set != influEdges.end();
//					iter_set++) {
//				if (i != *iter_set) {
//					line += to_string(*iter_set) + ",";
//				}
//			}
//			line.pop_back();
//			outfile_edges << line << endl;
//		}
	} else {
		while (getline(file, line)) {
			vector<string> line_vec;
			split(line, line_vec, ';');
			int edgeId = stoi(line_vec[0]);
			vector<string> line_vec_1;
			set<Path> influPaths;
			if (line_vec.size() > 1) {
				split(line_vec[1], line_vec_1, '*');
				int pathLen = line_vec_1.size();
				for (int i = 0; i < pathLen; i++) {
					string path_str = line_vec_1[i];
					vector<string> line_vec_2;
					double dis;
					list<int> i_edges;
					if (path_str.size() != 0) {
						split(path_str, line_vec_2, '-');
						istringstream iss(line_vec_2[0]);
						iss >> dis;
						vector<string> line_vec_3;
						split(line_vec_2[1], line_vec_3, ',');
						int edgesLen = line_vec_3.size();
						for (int j = 0; j < edgesLen; j++) {
							int i_edgeId = stoi(line_vec_3[j]);
							i_edges.push_back(i_edgeId);
						}
					}
					Path path(dis, i_edges);
					influPaths.insert(path);
				}
			}
			edge_influPaths[edgeId] = influPaths;
		}
	}
	return edge_influPaths;
}

map<int, set<int>> getEdgeInfluEdges(string input_file) {
	string line;
	ifstream file(input_file);
	map<int, set<int>> edge_influEdges;

	if (!file) {
		cout << "The file can not open!" << endl;
	} else {
		while (getline(file, line)) {
			vector<string> line_vec;
			split(line, line_vec, ';');
			int edgeId = stoi(line_vec[0]);
			vector<string> line_vec_1;
			set<int> influEdges;
			if (line_vec.size() > 1) {
				split(line_vec[1], line_vec_1, ',');
				int edgesLen = line_vec_1.size();
				for (int i = 0; i < edgesLen; i++) {
					influEdges.insert(stoi(line_vec_1[i]));
				}
			}
			edge_influEdges[edgeId] = influEdges;
		}
	}
	return edge_influEdges;
}

void writeOutputFile(string output_file, string line) {
	ofstream output(output_file, ios::app);
	output << line << endl;
}

void getRoadforPartition(vector<Edge> edges, map<int, double> edge_density_map,
		string output_road_for_partition) {
	int edgesNum = edges.size();
//	cout << edgesNum << " " << edge_density_map.size() << endl;
	map<int, string> vertexAndNeighbors;
	map<string, double> edgesAndDensity_old;
	map<string, double> edgesAndDensity_new;

	vector<Edge>::iterator iter_vec;
	for (iter_vec = edges.begin(); iter_vec != edges.end(); iter_vec++) {
		Edge edge = *iter_vec;
		int edgeId = edge.edgeId;
		int sourceNode = edge.sourceNode;
		int desNode = edge.desNode;
		double density = edge_density_map[edgeId];
		string key = to_string(sourceNode) + ";" + to_string(desNode);
		if (density == 0) {
			density = 1;
		}
		edgesAndDensity_old[key] = density;
	}

	map<string, double>::iterator iter_map_old;
	for (iter_map_old = edgesAndDensity_old.begin();
			iter_map_old != edgesAndDensity_old.end(); iter_map_old++) {
		string map_key = iter_map_old->first;
		double map_value = iter_map_old->second;
		edgesAndDensity_new[map_key] = map_value;
		vector<string> key_vec;
		split(map_key, key_vec, ';');
		string new_map_key = key_vec[1] + ";" + key_vec[0];
		if (edgesAndDensity_old.find(new_map_key)
				== edgesAndDensity_old.end()) {
			edgesAndDensity_new[new_map_key] = 1;
		} else {
			edgesAndDensity_new[new_map_key] = edgesAndDensity_old[new_map_key];
		}
	}
	cout << edgesAndDensity_old.size() << " " << edgesAndDensity_new.size()
			<< endl;

	map<string, double>::iterator iter_map_new;
	for (iter_map_new = edgesAndDensity_new.begin();
			iter_map_new != edgesAndDensity_new.end(); iter_map_new++) {
		vector<string> key_vec_new;
		string map_key_new = iter_map_new->first;
		split(map_key_new, key_vec_new, ';');
		int key = stoi(key_vec_new[0]);
		string value;
		if (vertexAndNeighbors.find(key) == vertexAndNeighbors.end()) {
			value = key_vec_new[1] + " 1 ";
		} else {
			value = vertexAndNeighbors[key] + key_vec_new[1] + " "
					+ to_string(int(iter_map_new->second)) + " ";
		}
		vertexAndNeighbors[key] = value;
	}

//	int nodesNum = vertexAndNeighbors.size();
	edgesNum = edgesAndDensity_new.size() / 2;
	int nodesNum = 4326;
	cout << "test: " << vertexAndNeighbors.size() << " "
			<< edgesAndDensity_new.size() << endl;
//	cout << vertexAndNeighbors[0] << endl;
//	cout << vertexAndNeighbors[731] << endl;
//	cout << vertexAndNeighbors[nodesNum] << endl;
	string line_start = to_string(nodesNum) + " " + to_string(edgesNum)
			+ " 001";
	writeOutputFile(output_road_for_partition, line_start);
	for (int i = 1; i <= nodesNum; i++) {
		writeOutputFile(output_road_for_partition, vertexAndNeighbors[i]);
	}
	cout << "The number of nodes: " << nodesNum << endl;
	cout << "The number of edges: " << edgesNum << endl;

	cout << "finish ... " << endl;
}

map<int, int> getNodePartition(string input_node_partition) {
	map<int, int> nodeAndPartition;
	string line;
	ifstream file(input_node_partition);
	int i = 1;
	if (!file) {
		cout << "The file can not open!" << endl;
	} else {
		while (getline(file, line)) {
			int partitionId = stoi(line);
			nodeAndPartition[i] = partitionId;
			i++;
		}
	}
	return nodeAndPartition;
}

void outputEdgesPartition(int partitionNum, string input_edge_partition,
		string output_partition) {
	map<int, int> edgeAndPartition;
	string line;
	ifstream file(input_edge_partition);
	int i = 1;
	if (!file) {
		cout << "The file can not open!" << endl;
	} else {
		while (getline(file, line)) {
			int partitionId = stoi(line);
			edgeAndPartition[i] = partitionId;
			i++;
		}
	}

	map<int, set<int>> partitionAndEdges;
	for (int j = 0; j < partitionNum; j++) {
		set<int> edges_set;
		partitionAndEdges[j] = edges_set;
	}
	map<int, int>::iterator iter_map;
	for (iter_map = edgeAndPartition.begin();
			iter_map != edgeAndPartition.end(); iter_map++) {
		int key = iter_map->first;
		int value = iter_map->second;
		partitionAndEdges[value].insert(key);
	}

	map<int, set<int>>::iterator iter_map_new;
	for (iter_map_new = partitionAndEdges.begin();
			iter_map_new != partitionAndEdges.end(); iter_map_new++) {
		int partitionId = iter_map_new->first;
		set<int> edges = iter_map_new->second;
		set<int>::iterator iter_set;
		string line = "";
		for (iter_set = edges.begin(); iter_set != edges.end(); iter_set++) {
			line += to_string(*iter_set) + ";";
		}
		line.pop_back();
		writeOutputFile(output_partition, line);
	}
	cout << "finish output ..." << endl;
}

void getRoadforPartitionOnDensity(vector<Edge> edges,
		map<int, double> edge_density_map, map<string, double> edge_pro_map,
		string output_road_for_partition) {
	int connectsNum = 0;
//	int edgesNum = edges.size();
	int edgesNum = 6135;
	cout << edgesNum << " " << edge_density_map.size() << endl;

	map<int, string> vertexAndNeighbors;
	map<int, int> edge_density_map_new;
	map<int, set<string>> edgeAndNeighborEdges;

	vector<Edge>::iterator iter_vec;
	for (int i = 0; i < edgesNum; i++) {
		int edgeId = i;
		map<int, double>::iterator iter_map = edge_density_map.find(edgeId);
		int density;
		if (iter_map == edge_density_map.end()) {
			density = 1;
		} else {
			density = (int) (edge_density_map[edgeId]);
			if (density == 0) {
				density = 1;
			}
		}
		edge_density_map_new[edgeId] = density;
		set<string> nEdges;
		edgeAndNeighborEdges[edgeId] = nEdges;
	}

	cout << edgeAndNeighborEdges.size() << " " << edge_density_map_new.size()
			<< " " << edge_pro_map.size() << endl;

	map<string, double>::iterator iter_map;
	for (iter_map = edge_pro_map.begin(); iter_map != edge_pro_map.end();
			iter_map++) {
		string edgeAndEdge = iter_map->first;
		vector<string> edgeAndEdge_vec;
		split(edgeAndEdge, edgeAndEdge_vec, ',');
		int edge_s = stoi(edgeAndEdge_vec[0]);
		int edge_d = stoi(edgeAndEdge_vec[1]);
		int pro = iter_map->second * 100;
//		cout << to_string(edge_s) << " " << to_string(edge_d) << " " << to_string(pro) << endl;
		string key_new = to_string(edge_d) + "," + to_string(edge_s);
		map<string, double>::iterator iter_map_new = edge_pro_map.find(key_new);
		string edgeDPro_new = to_string(edge_s) + ";";
		if (iter_map_new == edge_pro_map.end()) {
			edgeDPro_new += "1";
		} else {
			edgeDPro_new += to_string(edge_pro_map[key_new]);
		}
		edgeAndNeighborEdges[edge_d].insert(edgeDPro_new);

		string edgeDPro = to_string(edge_d) + ";" + to_string(pro);
		edgeAndNeighborEdges[edge_s].insert(edgeDPro);
		connectsNum++;
	}
	int sum = 0;
	string line_start = to_string(edgesNum) + " " + to_string(connectsNum)
			+ " 011";
	writeOutputFile(output_road_for_partition, line_start);

	for (int i = 0; i < edgesNum; i++) {
		int edgeId = i;
		set<string> edgeDPro_set = edgeAndNeighborEdges[edgeId];
		string line = to_string(edge_density_map_new[edgeId]) + " ";
		if (edgeDPro_set.size() != 0) {
			set<string>::iterator iter_set;
			for (iter_set = edgeDPro_set.begin();
					iter_set != edgeDPro_set.end(); iter_set++) {
				string edgeDPro = *iter_set;
				vector<string> edgeDPro_vec;
				split(edgeDPro, edgeDPro_vec, ';');
				int edge_d = stoi(edgeDPro_vec[0]);
				int pro = stoi(edgeDPro_vec[1]);
				if (pro == 0) {
					pro = 1;
				}
				line += to_string(edge_d) + " " + to_string(pro) + " ";
				sum++;
			}
		}
		writeOutputFile(output_road_for_partition, line);
	}
	cout << "sum: " << sum << endl;
//	int nodesNum = 4326;
//	cout << "The number of nodes: " << nodesNum << endl;
	cout << "The number of edges: " << edgesNum << endl;
	cout << "finish ... " << endl;
}

//int main() {
//	string city = "chengdu";
//	string hour = "1";
//	string input_road = "dataset/" + city + "/road.txt";
//	vector<Edge> edges = readEdgeBasicFile(input_road);
//	string input_edge_density = "dataset/" + city + "/edge_density/" + hour
//			+ ".txt";
//	map<int, double> edge_density_map = readEdgeDensityFile(input_edge_density);
//	string output_road_for_partition = "dataset/" + city + "/partition/" + hour
//			+ ".txt";
//
//	getRoadforPartition(edges, edge_density_map, output_road_for_partition);
//}
