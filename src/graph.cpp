#include <iostream>
#include <string>
#include <map>
#include <list>
#include <vector>
#include <set>
#include <algorithm>

#include "path.h"
#include "pathStr.h"
#include "graph.h"
#include "edge.h"
#include "fileProcess.h"
using namespace std;

Graph::Graph(int V) {
	this->V = V;
//	adj = new list<int> [V];
}

void Graph::addEdge(int v, int w) {
	adj[v].push_back(w); // Add w to v’s list.
}

/** This function is used to compute the influenced edges by the edgeId
 * edgeId: the edgeId
 * disPara: control the spatial diffusion range
 * edge_sourceNode_map: key->edgeId(int), value->sourceNode(int)
 * edge_desNode_map: key->edgeId(int), value->desNode(int)
 * nodes_edge_map: key->node1,node2(string), value->edgeId(int)
 * edge_dis_map: key->edgeId(int), value->dis(double)
 * edge_density_map: key->edgeId(int), value->traffic volume(int)
 * edge_pro_map: key->edgeId1,edgeId2(string), pro->probability(double)
 */
set<string> Graph::BFS_Str(int edgeId, int spreadWindow,
		map<int, int> edge_sourceNode_map, map<int, int> edge_desNode_map,
		map<string, int> nodes_edge_map, map<int, double> edge_dis_map,
		map<int, double> edge_time_map) {
	set<Path> results;
	list<Path> paths;
	list<int> edges;

	edges.push_back(edgeId);
	Path path(0, edges);
	paths.push_back(path);

	list<int>::iterator it;
	clock_t time1, time2, time3;
	time1 = clock();
	while (!paths.empty()) {
//		cout << paths.size() << endl;
		Path p = paths.front();
		paths.pop_front();
		double t_dis = p.dis; // dis is actually saved as time
		list<int> t_edges = p.edges;
		int t_end = t_edges.back();
		int desNode = edge_desNode_map[t_end];
		for (it = adj[desNode].begin(); it != adj[desNode].end(); ++it) {
			string key = to_string(desNode) + ',' + to_string(*it);
			if (nodes_edge_map.find(key) != nodes_edge_map.end()) {
				int new_edgeId = nodes_edge_map[key];
				double new_time = edge_time_map[new_edgeId];
//				cout << "new_dis: " << new_dis << " " << new_edgeId << endl;
				list<int>::iterator t_it = find(t_edges.begin(), t_edges.end(),
						new_edgeId);
				if (t_it == t_edges.end()) {
					if (t_dis + new_time > spreadWindow) {
						if (p.edges.size() > 1) {
							results.insert(p);
						}
					} else {
						Path t_p = p;
						t_p.dis += new_time;
						t_p.edges.push_back(new_edgeId);
						paths.push_back(t_p);
					}
				}
			}
		}
		if (adj[desNode].size() == 0 && p.edges.size() > 1) { // there is no further path to go through
			results.insert(p);
		}
	}

	time2 = clock();
	double duration0 = (double) (time2 - time1) / CLOCKS_PER_SEC;
//	cout << duration0 << endl;

	set<Path>::iterator iter_set = results.begin();

//	do the complements: because of the edge_pro maybe not all direct-connected edges are included in the path
	set<Path> new_results = results;
	set<string> pathTime;
	for (iter_set = results.begin(); iter_set != results.end(); iter_set++) {
		Path path = *iter_set;
		list<int> edges = path.edges;
		list<int>::iterator iter_list = edges.begin();
		list<int>::iterator iter_list_next = edges.begin();
		string pathStr = "";
		for (iter_list = edges.begin(); iter_list != edges.end(); iter_list++) {
			iter_list_next++;
			if (iter_list_next != edges.end()) {
				pathStr += to_string(*iter_list) + ",";
				string path = pathStr + to_string(*iter_list_next);
				pathTime.insert(path);
				int edgeId_cur = *iter_list;
//				int edgeId_next = *iter_list_next;
				int desNode = edge_desNode_map[edgeId_cur];
				for (it = adj[desNode].begin(); it != adj[desNode].end();
						++it) {
					int edgeId_neighbor = nodes_edge_map[to_string(desNode)
							+ ',' + to_string(*it)];
					path = pathStr + to_string(edgeId_neighbor);
					pathTime.insert(path);
				}
			}
		}
	}
	time3 = clock();
	double duration1 = (double) (time3 - time2) / CLOCKS_PER_SEC;
//	cout << duration1 << endl;
	return pathTime;
}

/** This function is used to compute the influenced edges by the edgeId
 * edgeId: the edgeId
 * disPara: control the spatial diffusion range
 * edge_sourceNode_map: key->edgeId(int), value->sourceNode(int)
 * edge_desNode_map: key->edgeId(int), value->desNode(int)
 * nodes_edge_map: key->node1,node2(string), value->edgeId(int)
 * edge_dis_map: key->edgeId(int), value->dis(double)
 * edge_density_map: key->edgeId(int), value->traffic volume(int)
 * edge_pro_map: key->edgeId1,edgeId2(string), pro->probability(double)
 */
set<Path> Graph::BFS(int edgeId, double disPara,
		map<int, int> edge_sourceNode_map, map<int, int> edge_desNode_map,
		map<string, int> nodes_edge_map, map<int, double> edge_dis_map) {
	set<Path> results;
	list<Path> paths;
	list<int> edges;

	edges.push_back(edgeId);
	Path path(0.0, edges);
	paths.push_back(path);

	list<int>::iterator it;

	while (!paths.empty()) {
		Path p = paths.front();
		paths.pop_front();
		double t_dis = p.dis;
		list<int> t_edges = p.edges;
		int t_end = t_edges.back();
//		int sourceNode = edge_sourceNode_map[t_end];
//		cout << "enter 1 ... " << endl;
		int desNode = edge_desNode_map[t_end];
//		bool flag = false;
//		cout << "enter 2 ... " << endl;
		for (it = adj[desNode].begin(); it != adj[desNode].end(); ++it) {
//			flag = true;
			int new_edgeId = nodes_edge_map[to_string(desNode) + ','
					+ to_string(*it)];

			double new_dis = edge_dis_map[new_edgeId];

			list<int>::iterator t_it = find(t_edges.begin(), t_edges.end(),
					new_edgeId); // iterator is used to save the idx

			if (t_it == t_edges.end()) {
				if (disPara != -1) { // which means that we have a bounded graph with the spatial diffusion range
					if (t_dis + new_dis > disPara) {
						if (p.edges.size() > 1) {
							results.insert(p);
						}
					} else {
						Path t_p = p;
						t_p.dis += new_dis;
						t_p.edges.push_back(new_edgeId);
						paths.push_back(t_p);
					}
				} else {
					Path t_p = p;
					t_p.dis += new_dis;
					t_p.edges.push_back(new_edgeId);
					paths.push_back(t_p);
				}
			}
		}
		if (adj[desNode].size() == 0 && p.edges.size() > 1) { // there is no further path to go through
			results.insert(p);
		}
	}

	return results;
}

map<int, double> Graph::getEdgeDensity_afterDiffusion_Pro(
		map<int, double> edge_density_map,
		map<int, map<int, double>> edge_pro_map,
		map<int, double> edge_residual_map,
		map<int, set<int>> r_edge_influEdges) {
	int edgesNum = edge_density_map.size();
	map<int, double> new_edge_density_map;
	for (int i = 0; i < edgesNum; i++) {
		set<int> influEdges = r_edge_influEdges[i];
		set<int>::iterator iter_set;
		for (iter_set = influEdges.begin(); iter_set != influEdges.end();
				iter_set++) {
			int influEdge = *iter_set;
			double pro = edge_pro_map[i][influEdge];
			new_edge_density_map[influEdge] += pro * edge_density_map[i];
		}
	}
	for (int i = 0; i < edgesNum; i++) {
		if (r_edge_influEdges[i].size() > 0) {
			new_edge_density_map[i] += edge_residual_map[i]
					* edge_density_map[i];
		} else {
			new_edge_density_map[i] += edge_density_map[i];
		}
	}

	return new_edge_density_map;
}

map<int, double> Graph::getEdgeDensity_afterDiffusion(
		map<int, double> edge_density_map, map<string, double> edge_pro_map,
		map<int, double> edge_residual_map,
		map<int, set<Path>> edge_influEdges) {
	int edgesNum = edge_density_map.size();
	map<int, double> new_edge_density_map = edge_density_map;

	for (int i = 0; i < edgesNum; i++) {
		set<Path> paths = edge_influEdges[i];
		set<Path>::iterator iter_set;
		string key;
		for (iter_set = paths.begin(); iter_set != paths.end(); iter_set++) {
			Path path = *iter_set;
			list<int> edges = path.edges;
			double pro = 1.0;
			list<int>::iterator iter_list = edges.begin();
			list<int>::iterator iter_list_next = edges.begin();
			iter_list_next++;
//			key = to_string(*iter_list) + "," + to_string(*iter_list_next);
//			update the density of current edge i
//			new_edge_density_map[i] -= edge_pro_map[key] * edge_density_map[i];
//			cout << "----" << new_edge_density_map[i] << " "
//					<< edge_pro_map[key] << " " << double(edge_density_map[i])
//					<< " " << edge_pro_map[key] * double(edge_density_map[i])
//					<< endl;
			for (iter_list = edges.begin(); iter_list_next != edges.end();
					iter_list++) {
				iter_list_next = iter_list;
				iter_list_next++;
				key = to_string(*iter_list) + "," + to_string(*iter_list_next);
				pro *= edge_pro_map[key];
//				update the density of influenced edges
//				cout << "++++ " << new_edge_density_map[*iter_list_next] << " "
//						<< pro << " " << double(edge_density_map[i]) << " "
//						<< pro * double(edge_density_map[i]) << endl;
				new_edge_density_map[*iter_list_next] += pro
						* edge_density_map[i];
//				cout << "++++ after " << new_edge_density_map[*iter_list_next] << endl;
			}
		}
	}
	return new_edge_density_map;
}

void Graph::getRoadforPartitionWithoutWeight(map<int, Edge> edges,
		string output_road_for_partition) {
	int edgesNum = edges.size();
	//	cout << edgesNum << " " << edge_density_map.size() << endl;
	map<int, set<int>> edgesAndNeighbors;
	set<int> nodes;
	map<string, int> edgesAndId;
	map<string, int> nodesAndEdgeId;
	for (int i = 0; i < edgesNum; i++) {
		Edge edge = edges[i];
		int edgeId = edge.edgeId;
		int sourceNode = edge.sourceNode;
		int desNode = edge.desNode;
		string key = to_string(sourceNode) + "," + to_string(desNode);
		nodesAndEdgeId[key] = edgeId;
	}

	int nodesNum = edgesNum;
	int edgesNumNew = 0;
	for (int i = 0; i < edgesNum; i++) {
		Edge edge = edges[i];
		int edgeId = edge.edgeId;
		int sourceNode = edge.sourceNode;
		int desNode = edge.desNode;
		nodes.insert(sourceNode);
		nodes.insert(desNode);
		list<int>::iterator it;
		for (it = adj[desNode].begin(); it != adj[desNode].end(); ++it) {
			string nextEdge = to_string(desNode) + ',' + to_string(*it);
			int nextEdgeId = nodesAndEdgeId[nextEdge];
			edgesAndNeighbors[edgeId].insert(nextEdgeId);
			edgesNumNew++;
		}
	}

	edgesNumNew = edgesNumNew / 2;
	string line_start = to_string(nodesNum) + " " + to_string(edgesNumNew);

	int sum = 0;
	writeOutputFile(output_road_for_partition, line_start);
	for (int i = 0; i < edgesNum; i++) {
		set<int> neighbors = edgesAndNeighbors[i];
		set<int>::iterator iter_set = neighbors.begin();
		string line;
		for (iter_set = neighbors.begin(); iter_set != neighbors.end();
				iter_set++) {
			line += to_string(*iter_set + 1) + " ";
			sum++;
		}
//		cout << line << endl;
		writeOutputFile(output_road_for_partition, line);
	}
	cout << "The number of nodes: " << nodesNum << endl;
	cout << "The number of edges: " << edgesNum << endl;
	cout << sum << endl;
	cout << edgesNumNew << endl;
	cout << "finish ... " << endl;
}
