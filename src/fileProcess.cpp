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

map<int, Edge> readEdgeBasicFile(int edgesNum, string input_file) {
	string line;
	map<int, Edge> edges;
	ifstream file(input_file);
	int maxEdgeId = 0;
	int maxNode = 0;
	int minEdgeId = 1000;
	int minNode = 1000;
	if (!file) {
		cout << "The file can not open!" << endl;
	} else {
		while (getline(file, line)) {
			vector<string> line_vec;
			split(line, line_vec, ',');
			int edgeId = stoi(line_vec[0]);
			int sourceNode = stoi(line_vec[1]);
			int desNode = stoi(line_vec[2]);
			if (edgeId > maxEdgeId) {
				maxEdgeId = edgeId;
			}
			if (sourceNode > maxNode) {
				maxNode = sourceNode;
			}
			if (desNode > maxNode) {
				maxNode = desNode;
			}
			if (edgeId < minEdgeId) {
				minEdgeId = edgeId;
			}
			if (sourceNode < minNode) {
				minNode = sourceNode;
			}
			if (desNode < minNode) {
				minNode = desNode;
			}
			istringstream iss(line_vec[3]);
			double dis;
			iss >> dis;
//			cout << edgeId << " " << sourceNode << " " << desNode << " " << dis << endl;
			Edge e;
			e.edgeId = edgeId;
			e.sourceNode = sourceNode;
			e.desNode = desNode;
			e.dis = dis;
			edges[edgeId] = e;
		}
	}
//	cout << edges.size() << endl;   // xian: 4969
	for (int i = 0; i < edgesNum; i++) {
		if (edges.find(i) == edges.end()) {
			Edge e;
			e.edgeId = -1;
			e.sourceNode = -1;
			e.desNode = -1;
			e.dis = -1;
			edges[i] = e;
		}
	}
//	cout << "max: " << maxEdgeId << "    " << maxNode << endl;
//	cout << "min: " << minEdgeId << "    " << minNode << endl;
	return edges;
}

set<int> readEdgeSmallFile(string input_file) {
	string line;
	ifstream file(input_file);
	set<int> edgeSmall;

	if (!file) {
		cout << "The file can not open!" << endl;
	} else {
		while (getline(file, line)) {
			int edgeId = stoi(line);
			edgeSmall.insert(edgeId);
		}
	}
	return edgeSmall;
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
			edge_density_map[edgeId - 1] = density; // because the density file of edgeId starts from 1
		}
	}
	return edge_density_map;
}

map<int, double> readEdgeSpeedFile(string input_file) {
	string line;
	ifstream file(input_file);
	map<int, double> edge_speed_map;

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
			edge_speed_map[edgeId - 1] = density; // because the density file of edgeId starts from 1
		}
	}
	return edge_speed_map;
}

map<string, double> readEdgeProFile(string input_edge_influPaths) {
	string line;
	ifstream file(input_edge_influPaths);
	map<string, double> edge_pro_map;

	if (!file) {
		cout << "The file can not open!" << endl;
	} else {
		while (getline(file, line)) {
			vector<string> line_vec;
			split(line, line_vec, ';');
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

map<int, map<int, double>> readEdgeNewProFile(string input_edge_influPaths) {
	string line;
	ifstream file(input_edge_influPaths);
	map<int, map<int, double>> edge_pro_map;

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
			//			string key = line_vec[0] + ',' + line_vec_1[0];
			int edge1 = stoi(line_vec[0]);
			int edge2 = stoi(line_vec_1[0]);
			istringstream iss(line_vec_1[1]);
			double pr;
			iss >> pr;
			edge_pro_map[edge1][edge2] = pr;
		}
	}

	return edge_pro_map;
}

map<int, double> readEdgeResidualFile(int edgesNum, string input_file) {
	string line;
	ifstream file(input_file);
	map<int, double> edge_residual_map;

	if (!file) {
		cout << "The file can not open!" << endl;
	} else {
		while (getline(file, line)) {
			vector<string> line_vec;
			split(line, line_vec, ';');
			int key = stoi(line_vec[0]);
			istringstream iss(line_vec[1]);
			double residual;
			iss >> residual;
			edge_residual_map[key] = residual;
		}
	}
	for (int i = 0; i < edgesNum; i++) {
		if (edge_residual_map.find(i) == edge_residual_map.end()) {
			edge_residual_map[i] = 1;
		}
	}
	return edge_residual_map;
}

/**
 * input_edge_influPaths: the influenced path for each edge
 * input_edge_influEdges: the influenced edges for each edge
 *
 * get the update pro between any two edges
 */
map<string, double> getEdgeInfluPaths(string input_file_paths,
		string input_file_edges, Graph graph, int edgesNum, double disPara,
		map<int, int> edge_sourceNode_map, map<int, int> edge_desNode_map,
		map<string, double> edge_pro_map, map<string, int> nodes_edge_map,
		map<int, double> edge_dis_map, map<int, double> edge_density_map) {
	string line;
	ifstream file(input_file_paths);
	map<int, set<Path>> edge_influPaths;
	map<int, set<int>> edge_influEdges;
	map<string, double> new_edge_pro_map = edge_pro_map;
	map<string, int> path_time; // to check whether the same path to a des node has been computed on the ratio repeatedly
	//	write into files
	ofstream outfile_paths(input_file_paths, ios::app);
	ofstream outfile_edges(input_file_edges, ios::app);

	if (!file) {
		cout << "Create edge influence paths file ..." << endl;
		for (int i = 0; i < edgesNum; i++) {
			set<Path> influPaths;
			set<int> influEdges;
			influPaths = graph.BFS(i, disPara, edge_sourceNode_map,
					edge_desNode_map, nodes_edge_map, edge_dis_map);
			set<Path>::iterator iter_set;
			cout << "edgeId is: " << i << " " << influPaths.size() << endl;
			for (iter_set = influPaths.begin(); iter_set != influPaths.end();
					iter_set++) {
				Path path = *iter_set;
				list<int>::iterator iter_list;
				list<int>::iterator iter_list_next = path.edges.begin();
				double pro = 0.0;
				double temp_pro = 1.0;
				if (path.edges.size() > 2) {
					string pathStr = to_string(*iter_list_next) + ",";
					for (iter_list = path.edges.begin();
							iter_list != path.edges.end(); iter_list++) {
						influEdges.insert(*iter_list);
						iter_list_next++;
						pathStr += to_string(*iter_list_next) + ",";
						if (iter_list == path.edges.begin()) {
							path_time[pathStr] = 1;
						}

						if (iter_list_next != path.edges.end()) {
							string temp_key = to_string(*iter_list) + ','
									+ to_string(*iter_list_next);
							if (edge_pro_map.find(temp_key)
									== edge_pro_map.end()) {
								temp_pro = 0.0;
							} else {
								temp_pro *= edge_pro_map[temp_key];
							}

							string key = to_string(i) + ','
									+ to_string(*iter_list_next);
							if (new_edge_pro_map.find(key)
									!= new_edge_pro_map.end()) {
								pro = new_edge_pro_map[key];
							} else {
								pro = 0.0;
							}

							if (path_time.find(pathStr) == path_time.end()) {
								new_edge_pro_map[key] = pro + temp_pro;
								path_time[pathStr] = 1;
							}
						}
					}
				}
			}

			edge_influEdges[i] = influEdges;
		}
		cout << "finish 1..." << new_edge_pro_map.size() << endl;

		//		output into the "edge_influPaths" file
		map<string, double>::iterator iter_map;
		iter_map = new_edge_pro_map.begin();
		string line_paths;
		while (iter_map != new_edge_pro_map.end()) {
			vector<string> key_vec;
			string key = iter_map->first;
			split(key, key_vec, ',');
			line_paths = key_vec[0] + ";" + key_vec[1] + ","
					+ to_string(iter_map->second);
			outfile_paths << line_paths << endl;
			iter_map++;
		}
		cout << "finish 1 ..." << endl;

		//		output into the "edge_influEdges" file
		map<int, set<int>> new_edge_influEdges;
		for (int i = 0; i < edgesNum; i++) {
			set<int> influEdges = edge_influEdges[i];
			set<int>::iterator iter_set_edges;
			string line_edges = to_string(i) + ";";
			for (iter_set_edges = influEdges.begin();
					iter_set_edges != influEdges.end(); iter_set_edges++) {
				new_edge_influEdges[*iter_set_edges].insert(i);
			}
		}

		for (int i = 0; i < edgesNum; i++) {
			set<int> influEdges = new_edge_influEdges[i];
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
	} else {
		//		new_edge_pro_map = readEdgeProFile(input_file_paths);
		//		map<string, double>::iterator iter_m = new_edge_pro_map.begin();
		//		int sum = 0;
		//		while (iter_m != new_edge_pro_map.end()) {
		//			if (iter_m->second == 0.0) {
		//				sum++;
		//			}
		//			iter_m++;
		//		}
		//		cout << sum << " " << new_edge_pro_map.size() << endl;
	}
	return new_edge_pro_map;
}

/**
 * input_edge_influPaths: the influenced path for each edge
 * input_edge_influEdges: the influenced edges for each edge
 *
 * get the update pro between any two edges
 */
map<int, map<int, double>> getEdgeInfluPaths(string input_edge_influPaths,
		string input_edge_influEdges, string input_edge_residual_new,
		Graph graph, int spreadWindow, map<int, int> edge_sourceNode_map,
		map<int, int> edge_desNode_map, map<string, double> edge_pro_map,
		map<string, int> nodes_edge_map, map<int, double> edge_dis_map,
		map<int, double> edge_residual_map, map<int, double> edge_time_map) {
	ifstream file(input_edge_influPaths);
	map<int, set<int>> edge_influEdges;
	map<int, map<int, double>> new_edge_pro_map;
	set<string> path_time;
	set<string> path_prev_time;
	//	write into files
	ofstream outfile_paths(input_edge_influPaths, ios::app);
	ofstream outfile_edges(input_edge_influEdges, ios::app);
	int edgesNum = edge_residual_map.size();

	if (!file) {
		cout << "Create edge influence paths file ..." << endl;
		for (int i = 0; i < edgesNum; i++) {
//		int test_edgeId = 1;
//		for (int i = test_edgeId; i < test_edgeId + 1; i++) {
			set<string> influPaths;
			set<int> influEdges;
			influPaths = graph.BFS_Str(i, spreadWindow, edge_sourceNode_map,
					edge_desNode_map, nodes_edge_map, edge_dis_map,
					edge_time_map);
			set<string>::iterator iter_set;
			for (iter_set = influPaths.begin(); iter_set != influPaths.end();
					iter_set++) {
				string path = *iter_set;
				vector<string> path_vec;
				split(path, path_vec, ',');
				string pathStr = path_vec[0] + ",";
				string pathStr_prev = "";
				for (int j = 0; j < path_vec.size() - 1; j++) {
					pathStr += path_vec[j + 1] + ",";
					pathStr_prev += path_vec[j] + ",";
					influEdges.insert(stoi(path_vec[j + 1]));

					if (path_time.find(pathStr) == path_time.end()) { // the path was not computed before
						path_time.insert(pathStr);
					}
					if (path_prev_time.find(pathStr_prev)
							== path_prev_time.end()) {
						path_prev_time.insert(pathStr_prev);
					}
				}
			}
			edge_influEdges[i] = influEdges;
			cout << "influEdges: " << i << " " << influEdges.size() << endl;
		}
		cout << "finish 1..." << endl;

//		This is to update the residual pro because of the spatial range, not all the edges which have the pro ratio are within the spatial range
		//		map<string, double> pathPro;
		set<string> pathDone;
		//		map<string, int>::iterator iter_map_str = path_time.begin();
		set<string>::iterator iter_set_str = path_time.begin();
		while (iter_set_str != path_time.end()) {
			string path = *iter_set_str;
			vector<string> line_vec;
			split(path, line_vec, ',');
			double temp_pro = 1.0;
			string pathStr = line_vec[0] + ",";
			string pathStr_pre = "";
			int pathLen = line_vec.size();
//			cout << "PATH: " << path << " " << pathLen << endl;
			for (int i = 1; i < pathLen; i++) {
				int edgeId_cur = stoi(line_vec[i - 1]);
				int edgeId_next = stoi(line_vec[i]);
				string temp_key = to_string(edgeId_cur) + ','
						+ to_string(edgeId_next);
				pathStr += line_vec[i] + ",";
				pathStr_pre += line_vec[i - 1] + ",";
				temp_pro *= edge_pro_map[temp_key]
						* (1 - edge_residual_map[edgeId_cur]);
				if (pathDone.find(pathStr) == pathDone.end()) { // this path hasn't been computed
//					cout << "enter: " << pathStr << " ** " << pathStr_pre
//							<< endl;
					double residual_pro;
					if (path_prev_time.find(pathStr) == path_prev_time.end()) { // the path is not the subset of other paths
						residual_pro = temp_pro;
//						cout << "temp_pro1: " << residual_pro << endl;
					} else {
						residual_pro = temp_pro
								* edge_residual_map[edgeId_next];
//						cout << "temp_pro2: " << residual_pro << endl;
					}
//					if (stoi(line_vec[0]) == test_edgeId) {
//						cout << pathStr << " \t" << stoi(line_vec[0]) << "\t"
//								<< edgeId_next << endl;
//						cout << edge_pro_map[temp_key] << " "
//								<< edge_residual_map[edgeId_cur] << " "
//								<< edge_residual_map[edgeId_next] << " "
//								<< temp_pro << endl;
//					}
					new_edge_pro_map[stoi(line_vec[0])][edgeId_next] +=
							residual_pro;
					pathDone.insert(pathStr);
					//					pathPro[pathStr] = temp_pro
					//							* (1 - edge_residual_map[edgeId_next]);
				}
				//				else {
				//					temp_pro = pathPro[pathStr];
				//				}
			}
			iter_set_str++;
		}

		cout << "finish 2..." << endl;

		//		output into the "edge_influPaths" file
		map<int, map<int, double>>::iterator iter_map;
		iter_map = new_edge_pro_map.begin();
		string line_paths;
		while (iter_map != new_edge_pro_map.end()) {
			//			vector<string> key_vec;
			int edge1 = iter_map->first;
			//			cout << edge1 << endl;
			map<int, double> temp_map = iter_map->second;
			map<int, double>::iterator temp_iter_map = temp_map.begin();
			//			cout << edge1 << endl;
			while (temp_iter_map != temp_map.end()) {
				//				cout << "enter" << endl;
				int edge2 = temp_iter_map->first;
				//				cout << edge2 << endl;
				double pro = temp_iter_map->second;
				//				cout << pro << endl;
				line_paths = to_string(edge1) + ";" + to_string(edge2) + ","
						+ to_string(pro);
				//				cout << line_paths << endl;
				outfile_paths << line_paths << endl;
				temp_iter_map++;
			}
			iter_map++;
		}
		cout << "finish 4 ..." << endl;

//		output into the "edge_influEdges" file
		map<int, set<int>> new_edge_influEdges; // in-neighbor
		for (int i = 0; i < edgesNum; i++) {
			set<int> influEdges = edge_influEdges[i];
			set<int>::iterator iter_set_edges;
			//			cout << i << ": ";
			for (iter_set_edges = influEdges.begin();
					iter_set_edges != influEdges.end(); iter_set_edges++) {
				new_edge_influEdges[*iter_set_edges].insert(i);
				//				cout << *iter_set_edges << " ";
			}
			//			cout << endl;
		}

		for (int i = 0; i < edgesNum; i++) {
			set<int> influEdges = new_edge_influEdges[i];
			set<int>::iterator iter_set_edges;
			string line_edges = to_string(i) + ";";
			//			cout << i << ": ";
			for (iter_set_edges = influEdges.begin();
					iter_set_edges != influEdges.end(); iter_set_edges++) {
				line_edges += to_string(*iter_set_edges) + ",";
				//				cout << *iter_set_edges << " ";
			}
			//			cout << endl;
			line_edges.pop_back();
			outfile_edges << line_edges << endl;
		}
		//		cout << "finish ..." << endl;
	} else {
		new_edge_pro_map = readEdgeNewProFile(input_edge_influPaths);
//		cout << edge_pro_map["1003,2172"] << endl;
//		cout << new_edge_pro_map[1003][2172] << endl;
	}
	return new_edge_pro_map;
}

/**
 * input_edge_influPaths: the influenced path for each edge
 * input_edge_influEdges: the influenced edges for each edge
 */
map<int, set<Path>> getEdgeInfluPaths_backup(string input_file_paths,
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

map<int, set<int>> getEdgeInfluEdges(string input_file, int edgesNum) {
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
			set<int> influEdges;
			vector<string> line_vec_1;
			if (line_vec.size() > 1) {
				split(line_vec[1], line_vec_1, ',');
				int edgesLen = line_vec_1.size();
				for (int i = 0; i < edgesLen; i++) {
					if (stoi(line_vec_1[i]) != edgeId) {
						influEdges.insert(stoi(line_vec_1[i]));
					}
				}
			}
			edge_influEdges[edgeId] = influEdges;
		}
	}
	return edge_influEdges;
}

map<int, set<int>> getReverseEdgeInfluEdges(int edgesNum, string input_file) {
	string line;
	ifstream file(input_file);
	map<int, set<int>> r_edge_influEdges;
	bool edgeEffect[edgesNum];
	for (int i = 0; i < edgesNum; i++) {
		edgeEffect[i] = false;
	}

	if (!file) {
		cout << "The file can not open!" << endl;
	} else {
		while (getline(file, line)) {
			vector<string> line_vec;
			split(line, line_vec, ';');
			int edgeId = stoi(line_vec[0]);
			vector<string> line_vec_1;
			if (line_vec.size() > 1) {
				split(line_vec[1], line_vec_1, ',');
				int edgesLen = line_vec_1.size();
				for (int i = 0; i < edgesLen; i++) {
					if (stoi(line_vec_1[i]) == true) {
						r_edge_influEdges[stoi(line_vec_1[i])].insert(edgeId);
					}
				}
			}
		}
	}
	return r_edge_influEdges;
}

void writeOutputFile(string output_file, string line) {
	ofstream output(output_file, ios::app);
	output << line << endl;
}

void getRoadforPartition(vector<Edge> edges, map<int, double> edge_density_map,
		string output_road_for_partition, int nodesNum) {
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
		int neighbor = stoi(key_vec_new[1]) + 1;
		if (vertexAndNeighbors.find(key) == vertexAndNeighbors.end()) {
			value = to_string(neighbor) + " 1 ";
		} else {
			value = vertexAndNeighbors[key] + to_string(neighbor) + " "
					+ to_string(int(iter_map_new->second)) + " ";
		}
		vertexAndNeighbors[key] = value;
	}

	//	int nodesNum = vertexAndNeighbors.size();
	int edgesNum = edgesAndDensity_new.size() / 2;
//	int nodesNum = 4326;
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
}

map<int, int> getNodePartition(string input_node_partition) {
	map<int, int> nodeAndPartition;
	string line;
	ifstream file(input_node_partition);
	int i = 1;
	if (!file) {
		cout << "The file can not open!" << endl;
		cout << input_node_partition << endl;
	} else {
		while (getline(file, line)) {
			int partitionId = stoi(line);
			nodeAndPartition[i] = partitionId;
			i++;
		}
	}
	return nodeAndPartition;
}

map<int, int> getEdgePartition(string input_edge_partition) {
	map<int, int> edgeAndPartition;
	string line;
	ifstream file(input_edge_partition);
	int i = 0;
	if (!file) {
		cout << "The file can not open!" << endl;
	} else {
		while (getline(file, line)) {
			int partitionId = stoi(line);
			edgeAndPartition[i] = partitionId;
			i++;
		}
	}
	return edgeAndPartition;
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
		string output_road_for_partition, int edgesNum) {
	int connectsNum = 0;
	//	int edgesNum = edges.size();
//	int edgesNum = 6135;
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
			edgeDPro_new += "1"; //edit
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
		string line;
		line = to_string(edge_density_map_new[edgeId]) + " ";
		if (edgeDPro_set.size() != 0) {
			set<string>::iterator iter_set;
			for (iter_set = edgeDPro_set.begin();
					iter_set != edgeDPro_set.end(); iter_set++) {
				string edgeDPro = *iter_set;
				vector<string> edgeDPro_vec;
				split(edgeDPro, edgeDPro_vec, ';');
				int edge_d = stoi(edgeDPro_vec[0]) + 1; //gpmetis can only start from 1
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
	cout << "The number of edges: " << edgesNum << endl;
	cout << "finish ... " << endl;
}

//int main() {
//	string city = "porto";
//	string hour = "11";
//	string input_road = "dataset/" + city + "/road.txt";
//	int nodesNum;
//	int edgesNum;
//	if (city == "chengdu") {
//		nodesNum = 4326;
//		edgesNum = 6135;
//	} else if (city == "xian") {
//		nodesNum = 2086;
//		edgesNum = 5045;
//	} else if (city == "porto") { // #traj: 1565595 # time range: -70 to 207 (278 totally)
//		nodesNum = 60287;
//		edgesNum = 108571;
//	}
//
//	map<int, Edge> edges_map = readEdgeBasicFile(edgesNum, input_road);
//	map<int, Edge>::iterator iter_map;
//	vector<Edge> edges;
//	for (iter_map = edges_map.begin(); iter_map != edges_map.end();
//			iter_map++) {
//		edges.push_back(iter_map->second);
//	}
//
//	string input_edge_density = "dataset/" + city + "/edge_density/" + hour
//			+ ".txt";
//	map<int, double> edge_density_map = readEdgeDensityFile(input_edge_density);
//	string output_road_for_partition = "dataset/" + city + "/partition/" + city
//			+ "_" + hour + ".txt";
//	string input_edge_pro = "dataset/" + city + "/edge_pro/" + hour + ".txt";
//
////	getRoadforPartition(edges, edge_density_map, output_road_for_partition,
////			nodesNum);
//	map<string, double> edge_pro_map;
//	edge_pro_map = readEdgeProFile(input_edge_pro);
//	getRoadforPartitionOnDensity(edges, edge_density_map, edge_pro_map,
//			output_road_for_partition, edgesNum);
//	cout << "finish ... " << endl;
//}
