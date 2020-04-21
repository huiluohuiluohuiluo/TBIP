#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <algorithm>

#include "time.h"
#include "edge.h"
#include "graph.h"
#include "fileProcess.h"
using namespace std;

vector<int> greedy(Graph graph, int K, double lambda,
		map<int, set<int>> edge_influEdges, map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, string output_file) {
	clock_t time0, time1, time2, time3, time4;
	map<int, double> edgeAndInfluScore;
//	int edgesNum = edge_density_map.size();
	int edgesNum = edge_influEdges.size();
//	cout << "check: " << edge_influEdges.size() << " "
//			<< edge_density_map.size() << endl;

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

	time0 = clock();
	double duration0, duration1 = 0.0;
	list<int> edges;
	for (int i = 0; i < edgesNum; i++) {
		edges.push_back(i);
	}
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
//			if (edgeAndInfluScore[tempEdge] > maxInfluScore) { // use the upper bound to prune
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
							tempInfluScore += new_edge_density_map[*iter_set]
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
//			}
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

	string runtime_g = "Greedy: " + to_string(duration0) + " "
			+ to_string(duration1) + " seconds" + "\n"
			+ to_string(result_influEdges.size()) + " "
			+ to_string(result_influScore);

	writeOutputFile(output_file, runtime_g);
	cout << runtime_g << endl;
	return result_seedEdges;
}

vector<int> greedyByInfluVal(Graph graph, int K, double lambda,
		map<int, set<int>> edge_influEdges, map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, string output_file) {
	clock_t time0, time1, time2, time3, time4;
	map<int, double> edgeAndInfluScore;
//	int edgesNum = edge_density_map.size();
	int edgesNum = edge_influEdges.size();

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

	time0 = clock();
	double duration0, duration1 = 0.0;
	list<int> edges;
	for (int i = 0; i < edgesNum; i++) {
		edges.push_back(i);
	}
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

	string runtime_g = "GreedyByInfluVal: " + to_string(duration0) + " "
			+ to_string(duration1) + " seconds" + "\n"
			+ to_string(result_influEdges.size()) + " "
			+ to_string(result_influScore);

	writeOutputFile(output_file, runtime_g);
	cout << runtime_g << endl;
	return result_seedEdges;
}

vector<int> greedyByMarginInfluVal(Graph graph, int K, double lambda,
		map<int, set<int>> edge_influEdges, map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, string output_file) {
	clock_t time0, time1, time2, time3, time4;
	map<int, double> edgeAndInfluScore;
	int edgesNum = edge_influEdges.size();

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

	time0 = clock();
	double duration0, duration1 = 0.0;
	list<int> edges;
	for (int i = 0; i < edgesNum; i++) {
		edges.push_back(i);
	}
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

		result_seedEdges.push_back(maxEdge);
		cout << maxEdge << endl;

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

		// update the array edgeAndInfluScore[] to get the new marginal gain for each edge
		for (iter_list = edges.begin(); iter_list != edges.end(); iter_list++) {
			int tempEdge = *iter_list;
			if (edgeAndInfluScore[tempEdge] > maxInfluScore) {

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
					time3 = clock();
					duration0 += (double) (time3 - time2) / CLOCKS_PER_SEC;
				}
			}
		}
	}

	time4 = clock();
	duration1 = (double) (time4 - time1) / CLOCKS_PER_SEC;

	string runtime_g = "GreedyByMarginInfluVal: " + to_string(duration0) + " "
			+ to_string(duration1) + " seconds" + "\n"
			+ to_string(result_influEdges.size()) + " "
			+ to_string(result_influScore);

	writeOutputFile(output_file, runtime_g);
	cout << runtime_g << endl;
	return result_seedEdges;
}

vector<int> greedyByMarginNumber(Graph graph, int K, double lambda,
		map<int, set<int>> edge_influEdges, map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, string output_file) {
	clock_t time0, time1, time2, time3, time4;
	map<int, double> edgeAndInfluScore;
	int edgesNum = edge_influEdges.size();

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

	time0 = clock();
	double duration0, duration1 = 0.0;
	list<int> edges;
	for (int i = 0; i < edgesNum; i++) {
		edges.push_back(i);
	}
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
//			if (edgeAndInfluScore[tempEdge] > maxInfluScore) { // use the upper bound to prune
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
//							tempInfluScore += new_edge_density_map[*iter_set]
//									- edge_density_map[*iter_set];
							tempInfluScore++;
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
//			}
		}

		result_seedEdges.push_back(maxEdge);
		cout << maxEdge << endl;

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

		// update the array edgeAndInfluScore[] to get the new marginal gain for each edge
//		for (iter_list = edges.begin(); iter_list != edges.end(); iter_list++) {
//			int tempEdge = *iter_list;
//			if (edgeAndInfluScore[tempEdge] > maxInfluScore) {
//				double tempInfluScore = 0;
//				if (edge_influEdges[tempEdge].size() != 0) {
//					set<int>::iterator iter_set;
//					time2 = clock();
//					for (iter_set = edge_influEdges[tempEdge].begin();
//							iter_set != edge_influEdges[tempEdge].end();
//							iter_set++) {
//						if (new_edge_density_map[*iter_set] > lambda
//								&& edge_density_map[*iter_set] >= 0
//								&& new_edge_density_map[*iter_set]
//										> edge_density_map[*iter_set]) {
//							set<int>::iterator iter_vec = find(
//									result_influEdges.begin(),
//									result_influEdges.end(), *iter_set); // iterator is used to save the idx
//
//							if (iter_vec == result_influEdges.end()) {
//								tempInfluScore +=
//										new_edge_density_map[*iter_set]
//												- edge_density_map[*iter_set];
//							}
//						}
//					}
//					edgeAndInfluScore[tempEdge] = tempInfluScore;
//					time3 = clock();
//					duration0 += (double) (time3 - time2) / CLOCKS_PER_SEC;
//				}
//			}
//		}
	}

	time4 = clock();
	duration1 = (double) (time4 - time1) / CLOCKS_PER_SEC;

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

	string runtime_g = "GreedyByMarginInfluVal: " + to_string(duration0) + " "
			+ to_string(duration1) + " seconds" + "\n"
			+ to_string(result_influEdges.size()) + " "
			+ to_string(result_influScore);

	writeOutputFile(output_file, runtime_g);
	cout << runtime_g << endl;
	return result_seedEdges;
}

map<int, vector<int>> getInfluencedEdges(vector<Edge> edges) {
	map<int, vector<int>> edgeInfluEdges;
	int edgeSize = edges.size();
	cout << edgeSize << endl;
	for (int i = 0; i < edgeSize; i++) {
		int edgeId = edges[i].edgeId;
		int sourceNode = edges[i].sourceNode;
		int desNode = edges[i].desNode;
//		double dis = edges[i].dis;
	}
	return edgeInfluEdges;
}
