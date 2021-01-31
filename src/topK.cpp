#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <algorithm>
#include <queue>

#include "time.h"
#include "edge.h"
#include "fileProcess.h"
#include "partition.h"
using namespace std;

void topK(int K, double lambda, map<int, set<int>> edge_influEdges,
		map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, string output_file) {
	map<int, double> edgeAndScore;
	int edgesNum = edge_density_map.size();
	double result_influScore = 0.0;

	for (int i = 0; i < edgesNum; i++) {
		set<int> influEdges = edge_influEdges[i];
		set<int>::iterator iter;
		double influScore = 0;
		for (iter = influEdges.begin(); iter != influEdges.end(); iter++) {
			if (new_edge_density_map[*iter] > lambda
					&& edge_density_map[*iter] >= 0
					&& new_edge_density_map[*iter] > edge_density_map[*iter]) {
				influScore += new_edge_density_map[*iter]
						- edge_density_map[*iter];
			}
		}
		edgeAndScore[i] = influScore;
	}

	vector<pair<int, double>> vtMap;
	for (auto it = edgeAndScore.begin(); it != edgeAndScore.end(); it++) {
		vtMap.push_back(make_pair(it->first, it->second));
	}

	sort(vtMap.begin(), vtMap.end(),
			[](const pair<int, double> &x, const pair<int, double> &y) -> int {
				return x.second > y.second;
			});

//	for (auto it = vtMap.begin(); it != vtMap.end(); it++)
//	     cout << it->first << ':' << it->second << '\n';

	set<int> result_influEdges;
	int j = 0;
	auto it = vtMap.begin();
	while (j < K) {
		int edge = it->first;
		cout << edge << " ";
		set<int> influEdges = edge_influEdges[edge];
		set<int>::iterator iter;
		for (iter = influEdges.begin(); iter != influEdges.end(); iter++) {
			if (new_edge_density_map[*iter] > lambda
					&& edge_density_map[*iter] >= 0
					&& new_edge_density_map[*iter] > edge_density_map[*iter]) {
				result_influEdges.insert(*iter);
			}
		}
		it++;
		j++;
	}

	set<int>::iterator iter_result;
	for (iter_result = result_influEdges.begin();
			iter_result != result_influEdges.end(); iter_result++) {
		result_influScore += new_edge_density_map[*iter_result]
				- edge_density_map[*iter_result];
	}
	string runtime_t = "TopK: " + to_string(result_influEdges.size()) + " "
			+ to_string(result_influScore);
	writeOutputFile(output_file, runtime_t);
	cout << runtime_t << endl;
}

void topK_minCut(int K, int M, double lambda,
		map<int, set<int>> edge_influEdges, map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, vector<Edge> raw_edges,
		string input_node_partition, string output_file) {
	clock_t time0, time1, time2, time3, time4;
	map<int, double> edgeAndInfluScore;

	int edgesNum = edge_influEdges.size();
	cout << "The number of edges is: " << edgesNum << endl;

	for (int i = 0; i < edgesNum; i++) {
		set<int> influEdges = edge_influEdges[i];
		set<int>::iterator iter;
		double influScore = 0;
		for (iter = influEdges.begin(); iter != influEdges.end(); iter++) {
			if (new_edge_density_map[*iter] > lambda
					&& edge_density_map[*iter] >= 0
					&& new_edge_density_map[*iter] > edge_density_map[*iter]) {
				influScore += new_edge_density_map[*iter]
						- edge_density_map[*iter];
			}
		}
		edgeAndInfluScore[i] = influScore;
	}

//	get the edge cuts after partitions
	list<int> edges = getPartitionEdges(M, raw_edges, input_node_partition);
//	cout << "The number of cuts is: " << edges.size() << endl;

	time0 = clock();
	double duration0, duration1 = 0.0;

	list<int>::iterator iter_list;
	vector<int> result_seedEdges;
	set<int> result_influEdges;
	double result_influScore = 0;
	while (result_seedEdges.size() < K) {
		int maxEdge;
		double maxInfluScore = 0;
		time1 = clock();
		for (iter_list = edges.begin(); iter_list != edges.end(); iter_list++) {
			int tempEdge = *iter_list;
			if (edgeAndInfluScore[tempEdge] > maxInfluScore) { // use the upper bound to prune
				double tempInfluScore = 0;
				if (edge_influEdges[tempEdge].size() != 0) {
					set<int>::iterator iter_set;
					time2 = clock();
					for (iter_set = edge_influEdges[tempEdge].begin();
							iter_set != edge_influEdges[tempEdge].end();
							iter_set++) {
						if (new_edge_density_map[*iter_set] > lambda
								&& edge_density_map[*iter_set] >= 0
								&& new_edge_density_map[*iter_set]
										> edge_density_map[*iter_set]) {
							set<int>::iterator iter_vec = find(
									result_influEdges.begin(),
									result_influEdges.end(), *iter_set); // iterator is used to save the idx

							if (iter_vec == result_influEdges.end()) {
								tempInfluScore +=
										new_edge_density_map[*iter_set]
												- edge_density_map[*iter_set];
							}
						}
					}
					if (tempInfluScore > maxInfluScore) {
						maxEdge = tempEdge;
						maxInfluScore = tempInfluScore;
					}
					time3 = clock();
					duration0 += (double) (time3 - time2) / CLOCKS_PER_SEC;
				}
			}
		}
		time4 = clock();
		duration1 += (double) (time4 - time1) / CLOCKS_PER_SEC;

		result_seedEdges.push_back(maxEdge);
		set<int>::iterator iter;
		for (iter = edge_influEdges[maxEdge].begin();
				iter != edge_influEdges[maxEdge].end(); iter++) {
			if (new_edge_density_map[*iter] > lambda
					&& edge_density_map[*iter] >= 0
					&& new_edge_density_map[*iter] > edge_density_map[*iter]) {
				result_influEdges.insert(*iter);
			}
		}
		result_influScore += maxInfluScore;
		edges.remove(maxEdge);
	}

	string runtime_tc = "TopK_minCut: " + to_string(duration0) + " "
			+ to_string(duration1) + " seconds" + "\n"
			+ to_string(result_influEdges.size()) + " "
			+ to_string(result_influScore);

	writeOutputFile(output_file, runtime_tc);
	cout << runtime_tc << endl;
//	return result_seedEdges;
}

