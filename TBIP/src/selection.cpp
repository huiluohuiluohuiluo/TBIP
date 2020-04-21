#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <algorithm>

#include "time.h"
#include "edge.h"
#include "graph.h"
#include "fileProcess.h"
#include "partition.h"
using namespace std;

map<int, double> edgeAndEstScore;
list<int> max_list;

class Cmpare {
public:
//	map<int, double> edgeAndEstScore;
//	map<int, double> edgeAndScore;
	bool operator()(int a, int b) {
		return edgeAndEstScore[a] > edgeAndEstScore[b];
	}
};

vector<int> selection(Graph graph, int K, int M, double disPara, double lambda,
		double theta, map<int, set<int>> old_edge_influEdges,
		map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, vector<Edge> vec_edges,
		string input_file, string input_node_partition, string output_file) {
	clock_t time0, time1, time2, time3;
	double duration1 = 0.0;
	time0 = clock();
//  Step1: get the partition info by kMeans
	map<int, set<int>> mClusters = getKMeansPartition(M, lambda, input_file,
			old_edge_influEdges, edge_density_map, new_edge_density_map);
	cout << "finish loading the clusters info: " << mClusters.size() << endl;

//	Step2: start the selection process
	vector<int> result_seedEdges;
	set<int> result_influEdges;
	double result_influScore = 0.0;
	map<int, double> clusterAndEstScore;
	map<int, list<int>> sorted_lists;
	map<int, int> edgeAndCluster;

	int old_edgesNum = old_edge_influEdges.size();
	ifstream file(input_file);
	map<int, set<int>> edge_influEdges;
	set<int> edges;
	for (int i = 0; i < old_edgesNum; i++) {
		if (old_edge_influEdges[i].size() > 0) { // filter some edges that can not influence other edges
			set<int> old_influEdges = old_edge_influEdges[i];
			set<int> influEdges;
			set<int>::iterator iter_influEdges;
			double influScore = 0.0;
			for (iter_influEdges = old_influEdges.begin();
					iter_influEdges != old_influEdges.end();
					iter_influEdges++) {
				if (new_edge_density_map[*iter_influEdges] > lambda
						&& edge_density_map[*iter_influEdges] >= 0
						&& new_edge_density_map[*iter_influEdges]
								> edge_density_map[*iter_influEdges]) {
					influEdges.insert(*iter_influEdges);
					influScore += new_edge_density_map[*iter_influEdges]
							- edge_density_map[*iter_influEdges];
				}
			}
			if (influEdges.size() > 0) {
				edges.insert(i);
				edge_influEdges[i] = influEdges;
				edgeAndEstScore[i] = influScore;
			}
		}
	}
	int edgesNum = edges.size();
	cout << "The number of edges is: " << edgesNum << endl;

	for (int m = 0; m < M; m++) {
		list<int> sorted_list;
		set<int> edgesInCluster = mClusters[m];
		set<int>::iterator iter_set_edgesInCluster;
		double estScoreInCluster = -1;
		for (iter_set_edgesInCluster = edgesInCluster.begin();
				iter_set_edgesInCluster != edgesInCluster.end();
				iter_set_edgesInCluster++) {
			int edgeIdInCluster = *iter_set_edgesInCluster;
			edgeAndCluster[edgeIdInCluster] = m;
			sorted_list.push_back(edgeIdInCluster);
			if (edgeAndEstScore[edgeIdInCluster] > estScoreInCluster) {
				estScoreInCluster = edgeAndEstScore[edgeIdInCluster];
			}
		}
		sorted_list.sort(Cmpare());
		sorted_lists[m] = sorted_list;
		clusterAndEstScore[m] = estScoreInCluster;
	}
	while (result_seedEdges.size() < K) {
//		select the cluster with maximal estimated marginal gain (line 6)
		double estScoreCluster_max = -1;
		int clusterId_max;
		for (int m = 0; m < M; m++) {
			double temp_estScoreCluster_max = clusterAndEstScore[m];
			if (temp_estScoreCluster_max > estScoreCluster_max) {
				estScoreCluster_max = temp_estScoreCluster_max;
				clusterId_max = m;
			}
		}

//		select the edge with maximal actual marginal gain in a maximal partition (line 7)
		set<int> edgesInCluster = mClusters[clusterId_max];
		set<int>::iterator iter_set_edgesInCluster;
		double edgeActScore_max = 0.0;
		int edgeId_max;
		for (iter_set_edgesInCluster = edgesInCluster.begin();
				iter_set_edgesInCluster != edgesInCluster.end();
				iter_set_edgesInCluster++) {
			int edgeIdInCluster = *iter_set_edgesInCluster;
			set<int> influEdges = edge_influEdges[edgeIdInCluster];
			set<int>::iterator iter_set_influEdge;
			double edgeActScore_temp = 0.0;
			for (iter_set_influEdge = influEdges.begin();
					iter_set_influEdge != influEdges.end();
					iter_set_influEdge++) {
				set<int>::iterator iter_set_find = find(
						result_influEdges.begin(), result_influEdges.end(),
						*iter_set_influEdge);
				if (iter_set_find == result_influEdges.end()) {
					edgeActScore_temp +=
							new_edge_density_map[*iter_set_influEdge]
									- edge_density_map[*iter_set_influEdge];
				}
			}
//			edgeAndActScore[edgeIdInCluster] = edgeActScore_temp;
			edgeAndEstScore[edgeIdInCluster] = edgeActScore_temp;
			if (edgeActScore_max < edgeActScore_temp) {
				edgeActScore_max = edgeActScore_temp;
				edgeId_max = edgeIdInCluster;
			}
		}

//		get other possible edges with high estimated scores (lines 8-10)
		for (int j = 0; j < M; j++) {
			if (j != clusterId_max) {
				list<int> sorted_list = sorted_lists[j];
				list<int>::iterator iter_list_s;
				for (iter_list_s = sorted_list.begin();
						iter_list_s != sorted_list.end(); iter_list_s++) {
					if (edgeAndEstScore[*iter_list_s] > edgeActScore_max) {
						max_list.push_back(*iter_list_s);
					}
				}
			}
		}
		max_list.sort(Cmpare());

//		get the actual maximal edge (lines 11-15)
		list<int>::iterator iter_list_m;
		for (iter_list_m = max_list.begin(); iter_list_m != max_list.end();
				iter_list_m++) {
			int est_max_edgeId = *iter_list_m;
			if (edgeActScore_max < edgeAndEstScore[est_max_edgeId]) {
				set<int> influEdges = edge_influEdges[est_max_edgeId];
				set<int>::iterator iter_set_influEdge;
				double edgeActScore_temp = 0.0;
				for (iter_set_influEdge = influEdges.begin();
						iter_set_influEdge != influEdges.end();
						iter_set_influEdge++) {
					set<int>::iterator iter_set_find = find(
							result_influEdges.begin(), result_influEdges.end(),
							*iter_set_influEdge);
					if (iter_set_find == result_influEdges.end()) {
						edgeActScore_temp +=
								new_edge_density_map[*iter_set_influEdge]
										- edge_density_map[*iter_set_influEdge];
					}
				}
//				update the upper bound of other edges with large scores (lines 16-18)
				edgeAndEstScore[est_max_edgeId] = edgeActScore_temp;
				if (edgeActScore_temp > edgeActScore_max) {
					edgeActScore_max = edgeActScore_temp;
					edgeId_max = est_max_edgeId;
					clusterId_max = edgeAndCluster[edgeId_max];
				}
			} else {
				break;
			}
		}

//		update the upper bound of each partition (not mentioned in the code skeleton)
		for (int m = 0; m < M; m++) {
			sorted_lists[m].sort(Cmpare());
			clusterAndEstScore[m] = edgeAndEstScore[*sorted_lists[m].begin()];
		}

//		push into seedEdges set
		result_seedEdges.push_back(edgeId_max);
//		cout << result_seedEdges.size() << ": " << clusterId_max << ": "
//				<< edgeId_max << ": " << result_influEdges.size() << ": "
//				<< edge_influEdges[edgeId_max].size() << endl;
		mClusters[clusterId_max].erase(edgeId_max);
		sorted_lists[clusterId_max].remove(edgeId_max);
		result_influScore += edgeActScore_max;
		set<int> influEdges = edge_influEdges[edgeId_max];
		set<int>::iterator iter_set_influEdge;
		for (iter_set_influEdge = influEdges.begin();
				iter_set_influEdge != influEdges.end(); iter_set_influEdge++) {
//			set<int>::iterator iter_set_find = find(result_influEdges.begin(),
//					result_influEdges.end(), *iter_set_influEdge);
//			if (iter_set_find == result_influEdges.end()) {
			result_influEdges.insert(*iter_set_influEdge);
//			}
		}

	} // end while loop to get K seed edges

	double score = 0.0;
	set<int> result_influEdges1;
	vector<int>::iterator iter_result;

	for (iter_result = result_seedEdges.begin();
			iter_result != result_seedEdges.end(); iter_result++) {
		int edgeId = *iter_result;
		set<int> influEdges = edge_influEdges[edgeId];
		set<int>::iterator iter_influEdges;
		for (iter_influEdges = influEdges.begin();
				iter_influEdges != influEdges.end(); iter_influEdges++) {
			int influEdgeId = *iter_influEdges;
			if (new_edge_density_map[influEdgeId] > lambda
					&& edge_density_map[influEdgeId] >= 0
					&& new_edge_density_map[influEdgeId]
							> edge_density_map[influEdgeId]) {
				set<int>::iterator iter_vec = find(result_influEdges1.begin(),
						result_influEdges1.end(), influEdgeId); // iterator is used to save the idx
				if (iter_vec == result_influEdges1.end()) {
					score += new_edge_density_map[influEdgeId]
							- edge_density_map[influEdgeId];
					result_influEdges1.insert(influEdgeId);
				}
			}
		}
	}
	cout << "actual score result is: " << score << endl;

	time1 = clock();
	double duration = (double) (time1 - time0) / CLOCKS_PER_SEC;
	string runtime_s = "Selection: " + to_string(duration) + " "
			+ to_string(duration1) + " seconds" + "\n"
			+ to_string(result_influEdges.size()) + " "
			+ to_string(result_influEdges1.size()) + " "
			+ to_string(result_influScore) + " " + to_string(score);

	writeOutputFile(output_file, runtime_s);
	cout << runtime_s << endl;
	return result_seedEdges;

}
