#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iterator>
#include <map>
#include <set>
#include <bits/stdc++.h>
#include <algorithm>

#include "time.h"
#include "edge.h"
#include "fileProcess.h"
#include "partition.h"

using namespace std;

vector<int> dp(Graph graph, int K, int M, double lambda, double theta,
		map<int, set<int>> edge_influEdges, map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, vector<Edge> edges,
		string input_file, string input_node_partition, string output_file) {
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

//	Method 1: follow the Ping's work
//	list<set<int>> mClusters = getThetaPartition(theta, lambda, input_file,
//			edge_influEdges, edge_density_map, new_edge_density_map);
//	M = mClusters.size();

//	Method 2: partition the graph by the traffic density/volume
//	list<set<int>> mClusters = getDensityPartition(M, edges,
//			input_node_partition);
//	map<int, int> edgeAndPartition = getEdgePartition(M, edges,
//			input_node_partition);

//	Method 3: partition the graph by the k-means method
	list<set<int>> mClusters = getKMeansPartition_list(input_file,
			edge_influEdges);
	map<int, int> edgeAndPartition = getEdgePartition_list(input_file,
			edge_influEdges);

	clock_t time0, time1, time2, time3;
	double duration1 = 0.0;
	time0 = clock();

	vector<int> result_seedEdges;
	set<int> result_influEdges;
	double result_influScore = 0;
//	score_arr[m][k] denotes the influence score of current found k seed edges from the first m partitions
	double score_arr[M + 1][K + 1];
//	bit_arr[m][k] means the partition id where the k-th seed edge resides in which cluster of the first m clusters
	int bit_arr[M + 1][K + 1];
	for (int k = 1; k < K + 1; k++) {
		score_arr[0][k] = 0;
		bit_arr[0][k] = 0;
	}
	for (int m = 1; m < M + 1; m++) {
		score_arr[m][0] = 0;
	}
	for (int k = 1; k < K + 1; k++) {
		list<set<int>>::iterator iter_list;
		int m = 0;
		int edgePos[M + 1];
		for (int m = 1; m < M + 1; m++) {
			edgePos[m] = 0;
		}
		for (iter_list = mClusters.begin(); iter_list != mClusters.end();
				iter_list++) {
			m++;
			set<int> edgesInCluster = *iter_list; // edges in a cluster
			set<int>::iterator key_iter_edge;
			int maxEdge_inner;
			double maxInfluScore_inner = 0;
			time2 = clock();
			for (key_iter_edge = edgesInCluster.begin();
					key_iter_edge != edgesInCluster.end(); key_iter_edge++) {
				int edgeId = *key_iter_edge;
				double tempInfluScore_inner = 0;
				// judge whether edgeId has been considered as a seed edge
				vector<int>::iterator iter_vec_result = find(
						result_seedEdges.begin(), result_seedEdges.end(),
						edgeId);
				if (iter_vec_result == result_seedEdges.end()) {
					set<int> influEdges = edge_influEdges[edgeId];
					set<int>::iterator iter_influEdges;
					for (iter_influEdges = influEdges.begin();
							iter_influEdges != influEdges.end();
							iter_influEdges++) {
						int influEdgeId = *iter_influEdges;
						if (edgeAndPartition[edgeId]
								== edgeAndPartition[influEdgeId]) {
							if (new_edge_density_map[influEdgeId] > lambda
									&& edge_density_map[influEdgeId] >= 0
									&& new_edge_density_map[influEdgeId]
											> edge_density_map[influEdgeId]) {
								set<int>::iterator iter_vec = find(
										result_influEdges.begin(),
										result_influEdges.end(), influEdgeId); // iterator is used to save the idx
								if (iter_vec == result_influEdges.end()) {
									tempInfluScore_inner +=
											new_edge_density_map[influEdgeId]
													- edge_density_map[influEdgeId];
									//								result_influEdges.insert(influEdgeId);
								}
							}
						}
					}
					if (tempInfluScore_inner > maxInfluScore_inner) {
						maxEdge_inner = edgeId;
						maxInfluScore_inner = tempInfluScore_inner;
					}
				}
			} // end the third loop: find the edge in a cluster to generate the maximal marginal gain
//			compute the actual score the selected maxEdge_inner
			double actual_maxInfluScore_inner = 0;
			// judge whether edgeId has been considered as a seed edge
			vector<int>::iterator iter_vec_result = find(
					result_seedEdges.begin(), result_seedEdges.end(),
					maxEdge_inner);
			if (iter_vec_result == result_seedEdges.end()) {
				set<int> influEdges = edge_influEdges[maxEdge_inner];
				set<int>::iterator iter_influEdges;
				for (iter_influEdges = influEdges.begin();
						iter_influEdges != influEdges.end();
						iter_influEdges++) {
					int influEdgeId = *iter_influEdges;
					if (new_edge_density_map[influEdgeId] > lambda
							&& edge_density_map[influEdgeId] >= 0
							&& new_edge_density_map[influEdgeId]
									> edge_density_map[influEdgeId]) {
						set<int>::iterator iter_vec = find(
								result_influEdges.begin(),
								result_influEdges.end(), influEdgeId); // iterator is used to save the idx
						if (iter_vec == result_influEdges.end()) {
							actual_maxInfluScore_inner +=
									new_edge_density_map[influEdgeId]
											- edge_density_map[influEdgeId];
							//								result_influEdges.insert(influEdgeId);
						}
					}
				}
			}

			time3 = clock();
			duration1 += (double) (time3 - time2) / CLOCKS_PER_SEC;
			if (score_arr[m - 1][k]
					>= score_arr[M][k - 1] + actual_maxInfluScore_inner) {
				score_arr[m][k] = score_arr[m - 1][k];
				bit_arr[m][k] = bit_arr[m - 1][k]; // record which partition
				edgePos[m] = edgePos[m - 1];	   // record which edge
			} else {
				score_arr[m][k] = score_arr[M][k - 1]
						+ actual_maxInfluScore_inner;
				bit_arr[m][k] = m;
				edgePos[m] = maxEdge_inner;
			}
			if (m == M) {
//				int partitionId = bit_arr[M][k];
//				cout << "judge: " << k << " " << maxEdge_inner << " "
//						<< edgePos[M] << endl;
				int edge_seed = edgePos[M];
				result_seedEdges.push_back(edge_seed);
				set<int> influEdges = edge_influEdges[edge_seed];
				set<int>::iterator iter_influEdges;
				for (iter_influEdges = influEdges.begin();
						iter_influEdges != influEdges.end();
						iter_influEdges++) {
					int influEdgeId = *iter_influEdges;
					if (new_edge_density_map[influEdgeId] > lambda
							&& edge_density_map[influEdgeId] >= 0
							&& new_edge_density_map[influEdgeId]
									> edge_density_map[influEdgeId]) {
						result_influEdges.insert(influEdgeId);
					}
				}
				cout << edge_seed << " " << result_influEdges.size() << endl;
			}
		} // end the second loop
	} // end the first loop

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

	result_influScore = score_arr[M][K];
	time1 = clock();
	double duration = (double) (time1 - time0) / CLOCKS_PER_SEC;
	string runtime_d = "DP: " + to_string(duration) + " " + to_string(duration1)
			+ " seconds" + "\n" + to_string(result_influEdges.size()) + " "
			+ to_string(result_influEdges1.size()) + " "
			+ to_string(result_influScore) + " " + to_string(score);

	writeOutputFile(output_file, runtime_d);
	cout << runtime_d << endl;

	return result_seedEdges;
}

vector<int> dpByInfluVal(Graph graph, int K, int M, double lambda, double theta,
		map<int, set<int>> edge_influEdges, map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, vector<Edge> edges,
		string input_file, string input_node_partition, string output_file) {
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

//	Method 1: follow the Ping's work
//	list<set<int>> mClusters = getThetaPartition(theta, lambda, input_file,
//			edge_influEdges, edge_density_map, new_edge_density_map);
//	M = mClusters.size();

//	Method 2: partition the graph by the traffic density/volume
//	list<set<int>> mClusters = getDensityPartition(M, edges,
//			input_node_partition);
//	map<int, int> edgeAndPartition = getEdgePartition(M, edges,
//			input_node_partition);

//	Method 3: partition the graph by the k-means method
	list<set<int>> mClusters = getKMeansPartition_list(input_file,
			edge_influEdges);
	map<int, int> edgeAndPartition = getEdgePartition_list(input_file,
			edge_influEdges);

	clock_t time0, time1, time2, time3;
	double duration1 = 0.0;
	time0 = clock();

	vector<int> result_seedEdges;
	set<int> result_influEdges;
	double result_influScore = 0;
//	score_arr[m][k] denotes the influence score of current found k seed edges from the first m partitions
	double score_arr[M + 1][K + 1];
//	bit_arr[m][k] means the partition id where the k-th seed edge resides in which cluster of the first m clusters
	int bit_arr[M + 1][K + 1];
	for (int k = 1; k < K + 1; k++) {
		score_arr[0][k] = 0;
		bit_arr[0][k] = 0;
	}
	for (int m = 1; m < M + 1; m++) {
		score_arr[m][0] = 0;
	}
	for (int k = 1; k < K + 1; k++) {
		list<set<int>>::iterator iter_list;
		int m = 0;
		int edgePos[M + 1];
		for (int m = 1; m < M + 1; m++) {
			edgePos[m] = 0;
		}
		for (iter_list = mClusters.begin(); iter_list != mClusters.end();
				iter_list++) {
			m++;
			set<int> edgesInCluster = *iter_list;
			set<int>::iterator key_iter_edge;
			int maxEdge_inner;
			double maxInfluScore_inner = 0;
			time2 = clock();
			for (key_iter_edge = edgesInCluster.begin();
					key_iter_edge != edgesInCluster.end(); key_iter_edge++) {
				int edgeId = *key_iter_edge;
				if (edgeAndInfluScore[edgeId] > maxInfluScore_inner) { // use the upper bound to prune
					double tempInfluScore_inner = 0;
					// judge whether edgeId has been considered as a seed edge
					vector<int>::iterator iter_vec_result = find(
							result_seedEdges.begin(), result_seedEdges.end(),
							edgeId);
					if (iter_vec_result == result_seedEdges.end()) {
						set<int> influEdges = edge_influEdges[edgeId];
						set<int>::iterator iter_influEdges;
						for (iter_influEdges = influEdges.begin();
								iter_influEdges != influEdges.end();
								iter_influEdges++) {
							int influEdgeId = *iter_influEdges;
							if (edgeAndPartition[edgeId]
									== edgeAndPartition[influEdgeId]) {
								if (new_edge_density_map[influEdgeId] > lambda
										&& edge_density_map[influEdgeId] >= 0
										&& new_edge_density_map[influEdgeId]
												> edge_density_map[influEdgeId]) {
									set<int>::iterator iter_vec = find(
											result_influEdges.begin(),
											result_influEdges.end(),
											influEdgeId); // iterator is used to save the idx
									if (iter_vec == result_influEdges.end()) {
										tempInfluScore_inner +=
												new_edge_density_map[influEdgeId]
														- edge_density_map[influEdgeId];
//									result_influEdges.insert(influEdgeId);
									}
								}
							}
						}
						if (tempInfluScore_inner > maxInfluScore_inner) {
							maxEdge_inner = edgeId;
							maxInfluScore_inner = tempInfluScore_inner;
						}
					}
				}
			} // end the third loop: find the edge in a cluster to generate the maximal marginal gain
//			compute the actual score the selected maxEdge_inner
			double actual_maxInfluScore_inner = 0;
			// judge whether edgeId has been considered as a seed edge
			vector<int>::iterator iter_vec_result = find(
					result_seedEdges.begin(), result_seedEdges.end(),
					maxEdge_inner);
			if (iter_vec_result == result_seedEdges.end()) {
				set<int> influEdges = edge_influEdges[maxEdge_inner];
				set<int>::iterator iter_influEdges;
				for (iter_influEdges = influEdges.begin();
						iter_influEdges != influEdges.end();
						iter_influEdges++) {
					int influEdgeId = *iter_influEdges;
					if (new_edge_density_map[influEdgeId] > lambda
							&& edge_density_map[influEdgeId] >= 0
							&& new_edge_density_map[influEdgeId]
									> edge_density_map[influEdgeId]) {
						set<int>::iterator iter_vec = find(
								result_influEdges.begin(),
								result_influEdges.end(), influEdgeId); // iterator is used to save the idx
						if (iter_vec == result_influEdges.end()) {
							actual_maxInfluScore_inner +=
									new_edge_density_map[influEdgeId]
											- edge_density_map[influEdgeId];
//								result_influEdges.insert(influEdgeId);
						}
					}
				}
			}

			time3 = clock();
			duration1 += (double) (time3 - time2) / CLOCKS_PER_SEC;
			if (score_arr[m - 1][k]
					>= score_arr[M][k - 1] + actual_maxInfluScore_inner) {
				score_arr[m][k] = score_arr[m - 1][k];
				bit_arr[m][k] = bit_arr[m - 1][k]; // record which partition
				edgePos[m] = edgePos[m - 1];	   // record which edge
			} else {
				score_arr[m][k] = score_arr[M][k - 1]
						+ actual_maxInfluScore_inner;
				bit_arr[m][k] = m;
				edgePos[m] = maxEdge_inner;
			}
			if (m == M) {
				int edge_seed = edgePos[M];
				result_seedEdges.push_back(edge_seed);
				set<int> influEdges = edge_influEdges[edge_seed];
				set<int>::iterator iter_influEdges;
				for (iter_influEdges = influEdges.begin();
						iter_influEdges != influEdges.end();
						iter_influEdges++) {
					int influEdgeId = *iter_influEdges;
					if (new_edge_density_map[influEdgeId] > lambda
							&& edge_density_map[influEdgeId] >= 0
							&& new_edge_density_map[influEdgeId]
									> edge_density_map[influEdgeId]) {
						result_influEdges.insert(influEdgeId);
					}
				}
				cout << edge_seed << " " << result_influEdges.size() << endl;
			}
		} // end the second loop
	} // end the first loop

	double score = 0.0;
	set<int> result_influEdges1;
	vector<int>::iterator iter_result;

	for (iter_result = result_seedEdges.begin();
			iter_result != result_seedEdges.end(); iter_result++) {
		int edgeId = *iter_result;
//		cout << "edgeId is: " << edgeId << endl;
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

	result_influScore = score_arr[M][K];
	time1 = clock();
	double duration = (double) (time1 - time0) / CLOCKS_PER_SEC;
	string runtime_d = "DPByInfluVal: " + to_string(duration) + " "
			+ to_string(duration1) + " seconds" + "\n"
			+ to_string(result_influEdges.size()) + " "
			+ to_string(result_influEdges1.size()) + " "
			+ to_string(result_influScore) + " " + to_string(score);

	writeOutputFile(output_file, runtime_d);
	cout << runtime_d << endl;

	return result_seedEdges;
}

vector<int> dpByMarginInfluVal(Graph graph, int K, int M, double lambda,
		double theta, map<int, set<int>> edge_influEdges,
		map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, vector<Edge> edges,
		string input_file, string input_node_partition, string output_file) {
	map<int, double> edgeAndInfluScore;
	int edgesNum = edge_influEdges.size();
	bool flags[edgesNum];
	for (int i = 0; i < edgesNum; i++) {
		flags[i] = false;
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

//	Method 1: follow the Ping's work
//	list<set<int>> mClusters = getThetaPartition(theta, lambda, input_file,
//			edge_influEdges, edge_density_map, new_edge_density_map);
//	M = mClusters.size();

//	Method 2: partition the graph by the traffic density/volume
//	list<set<int>> mClusters = getDensityPartition(M, edges,
//			input_node_partition);
//	map<int, int> edgeAndPartition = getEdgePartition(M, edges,
//			input_node_partition);

//	Method 3: partition the graph by the k-means method
	list<set<int>> mClusters = getKMeansPartition_list(input_file,
			edge_influEdges);
	map<int, int> edgeAndPartition = getEdgePartition_list(input_file,
			edge_influEdges);

	clock_t time0, time1, time2, time3;
	double duration1 = 0.0;
	time0 = clock();

	vector<int> result_seedEdges;
	set<int> result_influEdges;
	double result_influScore = 0;
//	score_arr[m][k] denotes the influence score of current found k seed edges from the first m partitions
	double score_arr[M + 1][K + 1];
//	bit_arr[m][k] means the partition id where the k-th seed edge resides in which cluster of the first m clusters
	int bit_arr[M + 1][K + 1];
	for (int k = 1; k < K + 1; k++) {
		score_arr[0][k] = 0;
		bit_arr[0][k] = 0;
	}
	for (int m = 1; m < M + 1; m++) {
		score_arr[m][0] = 0;
	}
	for (int k = 1; k < K + 1; k++) {
		list<set<int>>::iterator iter_list;
		int m = 0;
		int edgePos[M + 1];
		for (int m = 1; m < M + 1; m++) {
			edgePos[m] = 0;
		}
		for (iter_list = mClusters.begin(); iter_list != mClusters.end();
				iter_list++) {
			m++;
			set<int> edgesInCluster = *iter_list;	// edges in a cluster
			set<int>::iterator key_iter_edge;
			int maxEdge_inner;
			double maxInfluScore_inner = 0;
			time2 = clock();
			for (key_iter_edge = edgesInCluster.begin();
					key_iter_edge != edgesInCluster.end(); key_iter_edge++) {
				int edgeId = *key_iter_edge;
				if (edgeAndInfluScore[edgeId] > maxInfluScore_inner) { // use the upper bound to prune
					double tempInfluScore_inner = 0;
					// judge whether edgeId has been considered as a seed edge
					vector<int>::iterator iter_vec_result = find(
							result_seedEdges.begin(), result_seedEdges.end(),
							edgeId);
					if (iter_vec_result == result_seedEdges.end()) {
						set<int> influEdges = edge_influEdges[edgeId];
						set<int>::iterator iter_influEdges;
						for (iter_influEdges = influEdges.begin();
								iter_influEdges != influEdges.end();
								iter_influEdges++) {
							int influEdgeId = *iter_influEdges;
							if (edgeAndPartition[edgeId]
									== edgeAndPartition[influEdgeId]) {
								if (new_edge_density_map[influEdgeId] > lambda
										&& edge_density_map[influEdgeId] >= 0
										&& new_edge_density_map[influEdgeId]
												> edge_density_map[influEdgeId]) {
									set<int>::iterator iter_vec = find(
											result_influEdges.begin(),
											result_influEdges.end(),
											influEdgeId); // iterator is used to save the idx
									if (iter_vec == result_influEdges.end()) {
										tempInfluScore_inner +=
												new_edge_density_map[influEdgeId]
														- edge_density_map[influEdgeId];
//										result_influEdges.insert(influEdgeId);
									}
								}
							}
						}
						if (tempInfluScore_inner > maxInfluScore_inner) {
							maxEdge_inner = edgeId;
							maxInfluScore_inner = tempInfluScore_inner;
						}
					}
				}
			} // end the third loop: find the edge in a cluster to generate the maximal marginal gain
//			compute the actual score the selected maxEdge_inner
			double actual_maxInfluScore_inner = 0;
			// judge whether edgeId has been considered as a seed edge
			vector<int>::iterator iter_vec_result = find(
					result_seedEdges.begin(), result_seedEdges.end(),
					maxEdge_inner);
			if (iter_vec_result == result_seedEdges.end()) {
				set<int> influEdges = edge_influEdges[maxEdge_inner];
				set<int>::iterator iter_influEdges;
				for (iter_influEdges = influEdges.begin();
						iter_influEdges != influEdges.end();
						iter_influEdges++) {
					int influEdgeId = *iter_influEdges;
					if (new_edge_density_map[influEdgeId] > lambda
							&& edge_density_map[influEdgeId] >= 0
							&& new_edge_density_map[influEdgeId]
									> edge_density_map[influEdgeId]) {
						set<int>::iterator iter_vec = find(
								result_influEdges.begin(),
								result_influEdges.end(), influEdgeId); // iterator is used to save the idx
						if (iter_vec == result_influEdges.end()) {
							actual_maxInfluScore_inner +=
									new_edge_density_map[influEdgeId]
											- edge_density_map[influEdgeId];
//								result_influEdges.insert(influEdgeId);
						}
					}
				}
			}
			time3 = clock();
			duration1 += (double) (time3 - time2) / CLOCKS_PER_SEC;
			if (score_arr[m - 1][k]
					>= score_arr[M][k - 1] + actual_maxInfluScore_inner) {
				score_arr[m][k] = score_arr[m - 1][k];
				bit_arr[m][k] = bit_arr[m - 1][k]; // record which partition
				edgePos[m] = edgePos[m - 1];	   // record which edge
			} else {
				score_arr[m][k] = score_arr[M][k - 1]
						+ actual_maxInfluScore_inner;
				bit_arr[m][k] = m;
				edgePos[m] = maxEdge_inner;
			}
			if (m == M) {
				int edge_seed = edgePos[M];
				result_seedEdges.push_back(edge_seed);
				cout << edge_seed << endl;
				flags[edge_seed] = true;
				set<int> influEdges = edge_influEdges[edge_seed];
				set<int>::iterator iter_influEdges;
				for (iter_influEdges = influEdges.begin();
						iter_influEdges != influEdges.end();
						iter_influEdges++) {
					int influEdgeId = *iter_influEdges;
					if (new_edge_density_map[influEdgeId] > lambda
							&& edge_density_map[influEdgeId] >= 0
							&& new_edge_density_map[influEdgeId]
									> edge_density_map[influEdgeId]) {
						result_influEdges.insert(influEdgeId);
					}
				}
//				cout << edge_seed << " " << result_influEdges.size() << endl;

				// update the array edgeAndInfluScore[] to get the new marginal gain for each edge
				vector<Edge>::iterator iter_vec_temp;
				for (iter_vec_temp = edges.begin();
						iter_vec_temp != edges.end(); iter_vec_temp++) {
					int tempEdge = iter_vec_temp->edgeId;
					if (flags[tempEdge] == false
							&& edgeAndInfluScore[tempEdge]
									> maxInfluScore_inner) {
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
							edgeAndInfluScore[tempEdge] = tempInfluScore;
						}
					}
				}
			}
		} // end the second loop
	} // end the first loop

//	compute the real influence score
	double score = 0.0;
	set<int> result_influEdges1;
	vector<int>::iterator iter_result;

	for (iter_result = result_seedEdges.begin();
			iter_result != result_seedEdges.end(); iter_result++) {
		int edgeId = *iter_result;
		cout << "edgeId is: " << edgeId << endl;
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

	result_influScore = score_arr[M][K];
	time1 = clock();
	double duration = (double) (time1 - time0) / CLOCKS_PER_SEC;
	string runtime_d = "DPByMarginInfluVal: " + to_string(duration) + " "
			+ to_string(duration1) + " seconds" + "\n"
			+ to_string(result_influEdges.size()) + " "
			+ to_string(result_influEdges1.size()) + " "
			+ to_string(result_influScore) + " " + to_string(score);
//	string runtime_d = "DPByMarginInfluVal: " + to_string(duration) + " "
//			+ to_string(duration1) + " seconds" + "\n"
//			+ to_string(result_influEdges.size()) + " "
//			+ to_string(result_influScore);

	writeOutputFile(output_file, runtime_d);
	cout << runtime_d << endl;

	return result_seedEdges;
}
