#include <iostream>
#include <string>
#include <map>
#include <list>
#include <vector>
#include <set>
#include <algorithm>

#include "path.h"
#include "graph.h"
#include "edge.h"
#include "fileProcess.h"
using namespace std;

Graph::Graph(int V) {
	this->V = V;
	adj = new list<int> [V];
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
		bool flag = false;
//		cout << "enter 2 ... " << endl;
		for (it = adj[desNode].begin(); it != adj[desNode].end(); ++it) {
			flag = true;
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
		if (flag == false && p.edges.size() > 1) {
			results.insert(p);
//			cout << "results.size() is: " << results.size() << " "
//					<< p.edges.size() << endl;
		}
	}
	return results;
}

map<int, double> Graph::getEdgeDensity_afterDiffusion(
		map<int, double> edge_density_map, map<string, double> edge_pro_map,
		map<int, set<Path>> edge_influEdges) {

	int edgesNum = edge_density_map.size();
	map<int, double> new_edge_density_map;
	for (int i = 0; i < edgesNum; i++) {
		new_edge_density_map[i] = 0.0;
	}

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

