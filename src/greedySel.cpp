#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <algorithm>
#include <math.h>

#include "time.h"
#include "edge.h"
#include "graph.h"
#include "fileProcess.h"
#include "globalSel.h"
#include "greedySel.h"

using namespace std;

set<int> *greedySel(Graph graph, int k, double lambda,
		map<int, double> edge_dis_map, map<int, set<int>> edge_influEdges,
		map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, set<int> s_hat, set<int> i_hat) {
//	clock_t time0, time1, time2, time3, time4;
	map<int, double> edgeAndInfluNum;
	int edgesNum = edge_influEdges.size();
//	time0 = clock();
	bool edges_flag[edgesNum];
	bool influEdges_flag[edgesNum];
//	double duration0, duration1 = 0.0;
	list<int> edges;
	for (int i = 0; i < edgesNum; i++) {
		edges.push_back(i);
		edges_flag[i] = false;
		influEdges_flag[i] = false;
	}

	set<int>::iterator iter_s_hat;
	set<int>::iterator iter_i_hat;
	for (iter_s_hat = s_hat.begin(); iter_s_hat != s_hat.end(); iter_s_hat++) {
		edges_flag[*iter_s_hat] = true;
	}
	for (iter_i_hat = i_hat.begin(); iter_i_hat != i_hat.end(); iter_i_hat++) {
		influEdges_flag[*iter_i_hat] = true;
	}

	list<int>::iterator iter_list;

	set<int> result_seedEdges;
	set<int> result_influEdges;
	set<int> result_influEdgesInSomeBatches = i_hat;
	int result_influScore = 0;
	while (result_seedEdges.size() < k) {
		int maxEdge;
		double maxInfluScore = 0;
//		time1 = clock();
		for (iter_list = edges.begin(); iter_list != edges.end(); iter_list++) {
			int tempEdge = *iter_list;
			if (edges_flag[tempEdge] == false) { // check whether the current tempEdge has been included in the previous batches
				double tempInfluScore = 0;
				if (edge_influEdges[tempEdge].size() != 0) {
					set<int>::iterator iter_set;
					for (iter_set = edge_influEdges[tempEdge].begin();
							iter_set != edge_influEdges[tempEdge].end();
							iter_set++) {
						if (new_edge_density_map[*iter_set]
								> lambda * edge_dis_map[*iter_set]
								&& edge_density_map[*iter_set] >= 0
								&& new_edge_density_map[*iter_set]
										> edge_density_map[*iter_set]
								&& influEdges_flag[*iter_set] == false) {
							tempInfluScore++;
						}
					}
					if (tempInfluScore > maxInfluScore) {
						maxEdge = tempEdge;
						maxInfluScore = tempInfluScore;
					}
				}
			}
		}

//		can be improved
		result_seedEdges.insert(maxEdge);
		edges_flag[maxEdge] = true;
		set<int>::iterator iter;
		for (iter = edge_influEdges[maxEdge].begin();
				iter != edge_influEdges[maxEdge].end(); iter++) {
			if (new_edge_density_map[*iter] > lambda * edge_dis_map[*iter]
					&& edge_density_map[*iter] >= 0
					&& new_edge_density_map[*iter] > edge_density_map[*iter]
					&& influEdges_flag[*iter] == false) {
				result_influEdges.insert(*iter);
				result_influEdgesInSomeBatches.insert(*iter);
				influEdges_flag[*iter] = true;
			}
		}
//		cout << maxEdge << "   " << maxInfluScore << endl;
		result_influScore += maxInfluScore;
		edges.remove(maxEdge);
	}
	set<int> *result = new set<int> [3];
	result[0] = result_seedEdges;
	result[1] = result_influEdges; // the influenced edges in the i-th batches
	result[2] = result_influEdgesInSomeBatches; // the influenced edges in the first i batches
	cout << result_influScore << " ; actual: " << result_influEdges.size()
			<< endl;
	return result;
}

set<int> *greedySelProgressive(Graph graph, int k, double lambda,
		map<int, double> edge_dis_map, map<int, set<int>> edge_influEdges,
		map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, set<int> s_hat, set<int> i_hat) {
	map<int, double> edgeAndInfluNum;
	int edgesNum = edge_influEdges.size();
	bool edges_flag[edgesNum];
	bool influEdges_flag[edgesNum];
	list<int> edges;
	for (int i = 0; i < edgesNum; i++) {
		edges.push_back(i);
		edges_flag[i] = false;
		influEdges_flag[i] = false;
	}

	set<int>::iterator iter_s_hat;
	set<int>::iterator iter_i_hat;
	for (iter_s_hat = s_hat.begin(); iter_s_hat != s_hat.end(); iter_s_hat++) {
		edges_flag[*iter_s_hat] = true;
	}
	for (iter_i_hat = i_hat.begin(); iter_i_hat != i_hat.end(); iter_i_hat++) {
		influEdges_flag[*iter_i_hat] = true;
	}

	list<int>::iterator iter_list;

	set<int> result_seedEdges;
	set<int> result_influEdges;
	set<int> result_influEdgesInSomeBatches = i_hat;
	int result_influScore = 0;
	int h;
	double alpha = 0.9;
	int i = 1;
	while (result_seedEdges.size() < k) {
		int maxEdge = -1;
		int maxInfluScore = 0;
		if (result_seedEdges.size() == 0) {
			for (iter_list = edges.begin(); iter_list != edges.end();
					iter_list++) {
				int tempEdge = *iter_list;
				if (edges_flag[tempEdge] == false) { // check whether the current tempEdge has been included in the previous batches
					double tempInfluScore = 0;
					if (edge_influEdges[tempEdge].size() != 0) {
						set<int>::iterator iter_set;
						for (iter_set = edge_influEdges[tempEdge].begin();
								iter_set != edge_influEdges[tempEdge].end();
								iter_set++) {
							if (new_edge_density_map[*iter_set]
									> lambda * edge_dis_map[*iter_set]
									&& edge_density_map[*iter_set] >= 0
									&& new_edge_density_map[*iter_set]
											> edge_density_map[*iter_set]
									&& influEdges_flag[*iter_set] == false) {
								tempInfluScore++;
							}
						}
						if (tempInfluScore >= maxInfluScore) {
							maxEdge = tempEdge;
							maxInfluScore = tempInfluScore;
						}
					}
				}
			}
			h = maxInfluScore;  // initialize h as the first maximal marginal gain
			cout << "h: " << h << endl;
		} else {
			for (iter_list = edges.begin(); iter_list != edges.end();
					iter_list++) {
				int tempEdge = *iter_list;
				if (edges_flag[tempEdge] == false) { // check whether the current tempEdge has been included in the previous batches
					double tempInfluScore = 0;
					if (edge_influEdges[tempEdge].size() != 0) {
						set<int>::iterator iter_set;
						for (iter_set = edge_influEdges[tempEdge].begin();
								iter_set != edge_influEdges[tempEdge].end();
								iter_set++) {
							if (new_edge_density_map[*iter_set]
									> lambda * edge_dis_map[*iter_set]
									&& edge_density_map[*iter_set] >= 0
									&& new_edge_density_map[*iter_set]
											> edge_density_map[*iter_set]
									&& influEdges_flag[*iter_set] == false) {
								tempInfluScore++;
							}
						}
						if (tempInfluScore >= (double) h * pow(alpha, i)) {
							maxEdge = tempEdge;
							maxInfluScore = tempInfluScore;
							break;
						}
					}
				}
			}
		}
//		can be improved
		if (maxEdge != -1) {
			result_seedEdges.insert(maxEdge);
			edges_flag[maxEdge] = true;
			set<int>::iterator iter;
			for (iter = edge_influEdges[maxEdge].begin();
					iter != edge_influEdges[maxEdge].end(); iter++) {
				if (new_edge_density_map[*iter] > lambda * edge_dis_map[*iter]
						&& edge_density_map[*iter] >= 0
						&& new_edge_density_map[*iter] > edge_density_map[*iter]
						&& influEdges_flag[*iter] == false) {
					result_influEdges.insert(*iter);
					result_influEdgesInSomeBatches.insert(*iter);
					influEdges_flag[*iter] = true;
				}
			}
//			cout << maxEdge << "   " << maxInfluScore << endl;
			result_influScore += maxInfluScore;
			edges.remove(maxEdge);
		} else{
			i++;
		}
	}
	set<int> *result = new set<int> [3];
	result[0] = result_seedEdges;
	result[1] = result_influEdges; // the influenced edges in the i-th batches
	result[2] = result_influEdgesInSomeBatches; // the influenced edges in the first i batches
	cout << result_influScore << " ; actual: " << result_influEdges.size()
			<< endl;
	return result;
}

//int getUnion(int tempEdge, double lambda, map<int, double> edge_dis_map,
//		map<int, set<int>> edge_influEdges, map<int, double> edge_density_map,
//		map<int, double> new_edge_density_map, set<int> result_influEdges,
//		set<int> result_influEdgesInSomeBatches) {
//	int increInfluScore;
//	map<int, int> hashMapInSomeBatches;
//
//	set<int>::iterator iter_set;
//	for (iter_set = result_influEdgesInSomeBatches.begin();
//			iter_set != result_influEdgesInSomeBatches.end(); iter_set++) {
//		hashMapInSomeBatches[*iter_set] = 0;
//	}
//
//	for (iter_set = edge_influEdges[tempEdge].begin();
//			iter_set != edge_influEdges[tempEdge].end(); iter_set++) {
//		if (new_edge_density_map[*iter_set] > lambda * edge_dis_map[*iter_set]
//				&& edge_density_map[*iter_set] >= 0
//				&& new_edge_density_map[*iter_set]
//						> edge_density_map[*iter_set]) {
//			hashMapInSomeBatches[*iter_set] = 0;
//		}
//	}
//
//	increInfluScore = hashMapInSomeBatches.size()
//			- result_influEdgesInSomeBatches.size();
//	return increInfluScore;
//}

//set<int> *greedySel(Graph graph, int k, double lambda,
//		map<int, double> edge_dis_map, map<int, set<int>> edge_influEdges,
//		map<int, double> edge_density_map,
//		map<int, double> new_edge_density_map, set<int> s_hat, set<int> i_hat) {
//	map<int, double> edgeAndInfluNum;
//	int edgesNum = edge_influEdges.size();
//	bool edges_flag[edgesNum];
//
//	list<int> edges;
//	for (int i = 0; i < edgesNum; i++) {
//		edges.push_back(i);
//		edges_flag[i] = false;
//	}
//
//	set<int>::iterator iter_s_hat;
//	set<int>::iterator iter_i_hat;
//	for (iter_s_hat = s_hat.begin(); iter_s_hat != s_hat.end(); iter_s_hat++) {
//		edges_flag[*iter_s_hat] = true;
//	}
//
//	list<int>::iterator iter_list;
//
//	set<int> result_seedEdges;
//	set<int> result_influEdges;
//	set<int> result_influEdgesInSomeBatches = i_hat;
//	int result_influScore = 0;
//	while (result_seedEdges.size() < k) {
//		cout << "enter" << endl;
//		int maxEdge;
//		set<int> maxInfluEdges;
//		set<int> maxInfluEdgesInSomeBatches;
//		double maxInfluScore = 0;
//		for (iter_list = edges.begin(); iter_list != edges.end(); iter_list++) {
//			int tempEdge = *iter_list;
//			if (edges_flag[tempEdge] == false
//					&& edge_influEdges[tempEdge].size() != 0) {
//				int tempInfluScore = getUnion(tempEdge, lambda, edge_dis_map,
//						edge_influEdges, edge_density_map, new_edge_density_map,
//						result_influEdges, result_influEdgesInSomeBatches);
//				if (tempInfluScore > maxInfluScore) {
//					maxEdge = tempEdge;
//					maxInfluScore = tempInfluScore;
//				}
//			}
//		}
//
////		can be improved
//		result_seedEdges.insert(maxEdge);
//		edges_flag[maxEdge] = true;
//		cout << maxEdge << "hhhh" << maxInfluScore << endl;
//		set<int>::iterator iter;
//		for (iter = edge_influEdges[maxEdge].begin();
//				iter != edge_influEdges[maxEdge].end(); iter++) {
//			if (new_edge_density_map[*iter] > lambda * edge_dis_map[maxEdge]
//					&& edge_density_map[*iter] >= 0
//					&& new_edge_density_map[*iter] > edge_density_map[*iter]) {
//				result_influEdges.insert(*iter);
//				result_influEdgesInSomeBatches.insert(*iter);
//			}
//		}
//		result_influScore += maxInfluScore;
//		edges.remove(maxEdge);
//	}
//	set<int> *result = new set<int> [3];
//	result[0] = result_seedEdges;
//	result[1] = result_influEdges; // the influenced edges in the i-th batches
//	result[2] = result_influEdgesInSomeBatches; // the influenced edges in the first i batches
//	cout << result_influScore << " ; actual: " << result_influEdges.size()
//			<< endl;
//	return result;
//}

set<int> *greedySelPruneOne(Graph graph, int k, double lambda,
		map<int, double> edge_dis_map, map<int, set<int>> edge_influEdges,
		map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, set<int> s_hat, set<int> i_hat) {
	clock_t time0, time1, time2, time3, time4;
//	map<int, double> edgeAndInfluScore;
	map<int, double> edgeAndInfluNum;
//	int edgesNum = edge_density_map.size();
	int edgesNum = edge_influEdges.size();

//	for pruning
	for (int i = 0; i < edgesNum; i++) {
		set<int> influEdges = edge_influEdges[i];
		set<int>::iterator iter;
		int influNum = 0;
		for (iter = influEdges.begin(); iter != influEdges.end(); iter++) {
			if (new_edge_density_map[*iter] > lambda * edge_dis_map[*iter]
					&& edge_density_map[*iter] >= 0
					&& new_edge_density_map[*iter] > edge_density_map[*iter]) {
				influNum++;
			}
		}
		edgeAndInfluNum[i] = influNum;
	}

	time0 = clock();
	double duration0, duration1 = 0.0;
	list<int> edges;
	for (int i = 0; i < edgesNum; i++) {
		edges.push_back(i);
	}
	list<int>::iterator iter_list;

	set<int> result_seedEdges;
	set<int> result_seedEdgesInSomeBatches = s_hat;
	set<int> result_influEdges;
	set<int> result_influEdgesInSomeBatches = i_hat;
	double result_influScore = 0;
	while (result_seedEdges.size() < k) {
		int maxEdge;
		double maxInfluScore = 0;
		time1 = clock();
		for (iter_list = edges.begin(); iter_list != edges.end(); iter_list++) {
			int tempEdge = *iter_list;
			set<int>::iterator s_hat_iter = find(
					result_seedEdgesInSomeBatches.begin(),
					result_seedEdgesInSomeBatches.end(), tempEdge);
			if (s_hat_iter == result_seedEdgesInSomeBatches.end()) { // check whether the current tempEdge has been included in the previous batches
				if (edgeAndInfluNum[tempEdge] > maxInfluScore) { // use the upper bound to prune
					double tempInfluScore = 0;
					if (edge_influEdges[tempEdge].size() != 0) {
						set<int>::iterator iter_set;
						time2 = clock();
						for (iter_set = edge_influEdges[tempEdge].begin();
								iter_set != edge_influEdges[tempEdge].end();
								iter_set++) {
							if (new_edge_density_map[*iter_set]
									> lambda * edge_dis_map[*iter_set]
									&& edge_density_map[*iter_set] >= 0
									&& new_edge_density_map[*iter_set]
											> edge_density_map[*iter_set]) {
								set<int>::iterator iter_vec = find(
										result_influEdgesInSomeBatches.begin(),
										result_influEdgesInSomeBatches.end(),
										*iter_set); // iterator is used to save the idx

								if (iter_vec
										== result_influEdgesInSomeBatches.end()) {
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
				}
			}
		}

		time4 = clock();
		duration1 += (double) (time4 - time1) / CLOCKS_PER_SEC;
//		can be improved
		result_seedEdges.insert(maxEdge);
		result_seedEdgesInSomeBatches.insert(maxEdge);
		set<int>::iterator iter;
		for (iter = edge_influEdges[maxEdge].begin();
				iter != edge_influEdges[maxEdge].end(); iter++) {
			if (new_edge_density_map[*iter] > lambda * edge_dis_map[maxEdge]
					&& edge_density_map[*iter] >= 0
					&& new_edge_density_map[*iter] > edge_density_map[*iter]) {
				result_influEdges.insert(*iter);
				result_influEdgesInSomeBatches.insert(*iter);
			}
		}
		result_influScore += maxInfluScore;
//		edges.remove(maxEdge);
	}
	set<int> *result = new set<int> [3];
	result[0] = result_seedEdges;
	result[1] = result_influEdges;
	result[2] = result_influEdgesInSomeBatches;
	cout << result_influScore << " ; actual: " << result_influEdges.size()
			<< endl;
	return result;
}

set<int> *greedySelPruneTwo(Graph graph, int k, double lambda,
		map<int, double> edge_dis_map, map<int, set<int>> edge_influEdges,
		map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, set<int> s_hat, set<int> i_hat) {
	clock_t time0, time1, time2, time3, time4, time5, time6;
	map<int, double> edgeAndInfluNum;
	int edgesNum = edge_influEdges.size();

//	for pruning
	for (int i = 0; i < edgesNum; i++) {
		set<int> influEdges = edge_influEdges[i];
		set<int>::iterator iter;
		int influNum = 0;
		for (iter = influEdges.begin(); iter != influEdges.end(); iter++) {
			if (new_edge_density_map[*iter] > lambda * edge_dis_map[*iter]
					&& edge_density_map[*iter] >= 0
					&& new_edge_density_map[*iter] > edge_density_map[*iter]) {
				influNum++;
			}
		}
		edgeAndInfluNum[i] = influNum;
	}

	time0 = clock();
	double duration0, duration1 = 0.0;
	list<int> edges;
	for (int i = 0; i < edgesNum; i++) {
		edges.push_back(i);
	}
	list<int>::iterator iter_list;

	set<int> result_seedEdges;
	set<int> result_influEdges;
	set<int> result_influEdgesInSomeBatches = i_hat;
	double result_influScore = 0;
	while (result_seedEdges.size() < k) {
		int maxEdge;
		double maxInfluScore = 0;
		time1 = clock();
		for (iter_list = edges.begin(); iter_list != edges.end(); iter_list++) {
			int tempEdge = *iter_list;
			set<int>::iterator s_hat_iter = find(s_hat.begin(), s_hat.end(),
					tempEdge);
			if (s_hat_iter == s_hat.end()) { // check whether the current tempEdge has been included in the previous batches
				if (edgeAndInfluNum[tempEdge] > maxInfluScore) { // use the upper bound to prune
					double tempInfluScore = 0;
					if (edge_influEdges[tempEdge].size() != 0) {
						set<int>::iterator iter_set;
						time2 = clock();
						for (iter_set = edge_influEdges[tempEdge].begin();
								iter_set != edge_influEdges[tempEdge].end();
								iter_set++) {
							if (new_edge_density_map[*iter_set]
									> lambda * edge_dis_map[*iter_set]
									&& edge_density_map[*iter_set] >= 0
									&& new_edge_density_map[*iter_set]
											> edge_density_map[*iter_set]) {
								set<int>::iterator iter_vec = find(
										result_influEdgesInSomeBatches.begin(),
										result_influEdgesInSomeBatches.end(),
										*iter_set); // iterator is used to save the idx

								if (iter_vec
										== result_influEdgesInSomeBatches.end()) {
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
				}
			}
		}

		time4 = clock();
		duration1 += (double) (time4 - time1) / CLOCKS_PER_SEC;
//		can be improved
		result_seedEdges.insert(maxEdge);
		set<int>::iterator iter;
		for (iter = edge_influEdges[maxEdge].begin();
				iter != edge_influEdges[maxEdge].end(); iter++) {
			if (new_edge_density_map[*iter] > lambda * edge_dis_map[maxEdge]
					&& edge_density_map[*iter] >= 0
					&& new_edge_density_map[*iter] > edge_density_map[*iter]) {
				result_influEdges.insert(*iter);
				result_influEdgesInSomeBatches.insert(*iter);
			}
		}
		result_influScore += maxInfluScore;
		edges.remove(maxEdge);

		// update the array edgeAndInfluNum[] to get the new marginal gain for each edge
		time5 = clock();
		for (iter_list = edges.begin(); iter_list != edges.end(); iter_list++) {
			int tempEdge = *iter_list;
			if (edgeAndInfluNum[tempEdge] > maxInfluScore) {
				double tempInfluScore = 0;
				if (edge_influEdges[tempEdge].size() != 0) {
					set<int>::iterator iter_set;
					time2 = clock();
					for (iter_set = edge_influEdges[tempEdge].begin();
							iter_set != edge_influEdges[tempEdge].end();
							iter_set++) {
						if (new_edge_density_map[*iter_set]
								> lambda * edge_dis_map[*iter_set]
								&& edge_density_map[*iter_set] >= 0
								&& new_edge_density_map[*iter_set]
										> edge_density_map[*iter_set]) {
							set<int>::iterator iter_vec = find(
									result_influEdges.begin(),
									result_influEdges.end(), *iter_set); // iterator is used to save the idx

							if (iter_vec == result_influEdges.end()) {
								tempInfluScore++;
//								tempInfluScore +=
//										new_edge_density_map[*iter_set]
//												- edge_density_map[*iter_set];
							}
						}
					}
					edgeAndInfluNum[tempEdge] = tempInfluScore;
				}
			}
		}
		time6 = clock();
		duration0 += (double) (time6 - time5) / CLOCKS_PER_SEC;
	}
	set<int> *result = new set<int> [3];
	result[0] = result_seedEdges;
	result[1] = result_influEdges;
	result[2] = result_influEdgesInSomeBatches;
	cout << result_influScore << " ; actual: " << result_influEdges.size()
			<< endl;
//	string runtime_g = "GreedySel: " + to_string(duration0) + " "
//			+ to_string(duration1) + " seconds" + "\n"
//			+ to_string(result_influEdges.size()) + " "
//			+ to_string(result_influScore);
//	cout << runtime_g << endl;
	return result;
}

set<int> *greedySelPruneThree(Graph graph, int k, double lambda,
		map<int, double> edge_dis_map, map<int, set<int>> edge_influEdges,
		map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, set<int> s_hat, set<int> i_hat) {
	clock_t time0, time1, time2, time3, time4, time5, time6;
	map<int, double> edgeAndInfluNum;
	int edgesNum = edge_influEdges.size();

//	for pruning
	for (int i = 0; i < edgesNum; i++) {
		set<int> influEdges = edge_influEdges[i];
		set<int>::iterator iter;
		int influNum = 0;
		for (iter = influEdges.begin(); iter != influEdges.end(); iter++) {
			if (new_edge_density_map[*iter] > lambda * edge_dis_map[*iter]
					&& edge_density_map[*iter] >= 0
					&& new_edge_density_map[*iter] > edge_density_map[*iter]) {
				influNum++;
			}
		}
		edgeAndInfluNum[i] = influNum;
	}

	time0 = clock();
	double duration0, duration1 = 0.0;
	list<int> edges;
	for (int i = 0; i < edgesNum; i++) {
		edges.push_back(i);
	}
	list<int>::iterator iter_list;

	set<int> result_seedEdges;
	set<int> result_influEdges;
	set<int> result_influEdgesInSomeBatches = i_hat;
	double result_influScore = 0;
	while (result_seedEdges.size() < k) {
		int maxEdge;
		double maxInfluScore = 0;
		time1 = clock();
		for (iter_list = edges.begin(); iter_list != edges.end(); iter_list++) {
			int tempEdge = *iter_list;
			set<int>::iterator s_hat_iter = find(s_hat.begin(), s_hat.end(),
					tempEdge);
			if (s_hat_iter == s_hat.end()) { // check whether the current tempEdge has been included in the previous batches
				if (edgeAndInfluNum[tempEdge] > maxInfluScore) { // use the upper bound to prune
					double tempInfluScore = 0;
					if (edge_influEdges[tempEdge].size() != 0) {
						set<int>::iterator iter_set;
						time2 = clock();
						for (iter_set = edge_influEdges[tempEdge].begin();
								iter_set != edge_influEdges[tempEdge].end();
								iter_set++) {
							if (new_edge_density_map[*iter_set]
									> lambda * edge_dis_map[*iter_set]
									&& edge_density_map[*iter_set] >= 0
									&& new_edge_density_map[*iter_set]
											> edge_density_map[*iter_set]) {
								set<int>::iterator iter_vec = find(
										result_influEdgesInSomeBatches.begin(),
										result_influEdgesInSomeBatches.end(),
										*iter_set); // iterator is used to save the idx

								if (iter_vec
										== result_influEdgesInSomeBatches.end()) {
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
				}
			}
		}

		time4 = clock();
		duration1 += (double) (time4 - time1) / CLOCKS_PER_SEC;
//				can be improved
		result_seedEdges.insert(maxEdge);
		set<int>::iterator iter;
		for (iter = edge_influEdges[maxEdge].begin();
				iter != edge_influEdges[maxEdge].end(); iter++) {
			if (new_edge_density_map[*iter] > lambda * edge_dis_map[maxEdge]
					&& edge_density_map[*iter] >= 0
					&& new_edge_density_map[*iter] > edge_density_map[*iter]) {
				result_influEdges.insert(*iter);
				result_influEdgesInSomeBatches.insert(*iter);
			}
		}
		result_influScore += maxInfluScore;
		edges.remove(maxEdge);

		// update the array edgeAndInfluNum[] to get the new marginal gain for each edge
		time5 = clock();
		for (iter_list = edges.begin(); iter_list != edges.end(); iter_list++) {
			int tempEdge = *iter_list;
			if (edgeAndInfluNum[tempEdge] > maxInfluScore) {
				double tempInfluScore = 0;
				if (edge_influEdges[tempEdge].size() != 0) {
					set<int>::iterator iter_set;
					time2 = clock();
					for (iter_set = edge_influEdges[tempEdge].begin();
							iter_set != edge_influEdges[tempEdge].end();
							iter_set++) {
						if (new_edge_density_map[*iter_set]
								> lambda * edge_dis_map[*iter_set]
								&& edge_density_map[*iter_set] >= 0
								&& new_edge_density_map[*iter_set]
										> edge_density_map[*iter_set]) {
							set<int>::iterator iter_vec = find(
									result_influEdgesInSomeBatches.begin(),
									result_influEdgesInSomeBatches.end(),
									*iter_set); // iterator is used to save the idx

							if (iter_vec == result_influEdges.end()) {
								tempInfluScore++;
//								tempInfluScore +=
//										new_edge_density_map[*iter_set]
//												- edge_density_map[*iter_set];
							}
						}
					}
					edgeAndInfluNum[tempEdge] = tempInfluScore;
				}
			}
		}
		time6 = clock();
		duration0 += (double) (time6 - time5) / CLOCKS_PER_SEC;
	}
	set<int> *result = new set<int> [2];
	result[0] = result_seedEdges;
	result[1] = result_influEdges;
	string runtime_g = "GreedySel: " + to_string(duration0) + " "
			+ to_string(duration1) + " seconds" + "\n"
			+ to_string(result_influEdges.size()) + " "
			+ to_string(result_influScore);
	cout << runtime_g << endl;
	return result;
}

vector<int> *greedySelUpperBoundPerBatch(Graph graph, int k, double lambda,
		map<int, double> edge_dis_map, map<int, set<int>> edge_influEdges,
		map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, set<int> s_hat, set<int> i_hat) {
	clock_t time0, time1, time2, time3, time4;
	map<int, double> edgeAndInfluNum;
	int edgesNum = edge_influEdges.size();
	time0 = clock();
	bool edges_flag[edgesNum];
	bool influEdges_flag[edgesNum];
	double duration0, duration1 = 0.0;
	list<int> edges;
	for (int i = 0; i < edgesNum; i++) {
		edges.push_back(i);
		edges_flag[i] = false;
		influEdges_flag[i] = false;
	}

	set<int>::iterator iter_s_hat;
	set<int>::iterator iter_i_hat;
	for (iter_s_hat = s_hat.begin(); iter_s_hat != s_hat.end(); iter_s_hat++) {
		edges_flag[*iter_s_hat] = true;
	}
	for (iter_i_hat = i_hat.begin(); iter_i_hat != i_hat.end(); iter_i_hat++) {
		influEdges_flag[*iter_i_hat] = true;
	}

	list<int>::iterator iter_list;

	vector<int> result_seedEdges;
	vector<int> result_influEdges;
//	set<int> result_influEdgesInSomeBatches = i_hat;

	map<int, int> edgeAndNum;
	map<int, set<int>> edgeAndInfluEdges;

	for (iter_list = edges.begin(); iter_list != edges.end(); iter_list++) {
		int tempEdge = *iter_list;
		if (edges_flag[tempEdge] == false) {
			double tempInfluScore = 0;
			if (edge_influEdges[tempEdge].size() != 0) {
				set<int>::iterator iter_set;
				time2 = clock();
				set<int> tempInfluEdges;
				for (iter_set = edge_influEdges[tempEdge].begin();
						iter_set != edge_influEdges[tempEdge].end();
						iter_set++) {
					if (new_edge_density_map[*iter_set]
							> lambda * edge_dis_map[*iter_set]
							&& edge_density_map[*iter_set] >= 0
							&& new_edge_density_map[*iter_set]
									> edge_density_map[*iter_set]) {
						if (influEdges_flag[*iter_set] == false) {
							tempInfluScore++;
							tempInfluEdges.insert(*iter_set);
						}
					}
				}
				edgeAndNum[tempEdge] = tempInfluScore;
				edgeAndInfluEdges[tempEdge] = tempInfluEdges;
			}
		}
	}

	vector<pair<int, double>> vtMap;
	for (auto it = edgeAndNum.begin(); it != edgeAndNum.end(); it++) {
		vtMap.push_back(make_pair(it->first, it->second));
	}
	sort(vtMap.begin(), vtMap.end(),
			[](const pair<int, int> &x, const pair<int, int> &y) -> int {
				return x.second > y.second;
			});

	int j = 0;
	auto it = vtMap.begin();
	while (j < k) {
		int edge = it->first;
		set<int> influEdges = edgeAndInfluEdges[edge];
		result_seedEdges.push_back(edge);
		set<int>::iterator iter;
		for (iter = influEdges.begin(); iter != influEdges.end(); iter++) {
			if (new_edge_density_map[*iter] > lambda * edge_dis_map[*iter]
					&& edge_density_map[*iter] >= 0
					&& new_edge_density_map[*iter] > edge_density_map[*iter]) {
				result_influEdges.push_back(*iter);
			}
		}
		it++;
		j++;
	}

	vector<int> *result = new vector<int> [3];
	result[0] = result_seedEdges;
	result[1] = result_influEdges; // the influenced edges in the i-th batches
//	result[2] = result_influEdgesInSomeBatches; // the influenced edges in the first i batches
	return result;
}

set<int> *greedySelSampling(Graph graph, int k, double lambda,
		map<int, double> edge_dis_map, map<int, set<int>> edge_influEdges,
		map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, set<int> s_hat, set<int> i_hat) {
	map<int, double> edgeAndInfluNum;
	int edgesNum = edge_influEdges.size();
	map<int, set<int>> r_edge_influEdges;

	bool edges_flag[edgesNum];
	bool influEdges_flag[edgesNum];
//	initialize
	for (int i = 0; i < edgesNum; i++) {
		set<int> r_influEdges;
		r_edge_influEdges[i] = r_influEdges;
	}
	map<int, set<int>> numArr;
	for (int i = 0; i < edgesNum; i++) {
		set<int> value;
		numArr[i] = value;
	}

	for (int i = 0; i < edgesNum; i++) {
		influEdges_flag[i] = false;
		edges_flag[i] = false;
		set<int> influEdges = edge_influEdges[i];
		set<int>::iterator iter;
		for (iter = influEdges.begin(); iter != influEdges.end(); iter++) {
			if (*iter < edgesNum
					&& new_edge_density_map[*iter]
							> lambda * edge_dis_map[*iter]
					&& edge_density_map[*iter] >= 0
					&& new_edge_density_map[*iter] > edge_density_map[*iter]) {
				r_edge_influEdges[*iter].insert(i);
				numArr[i].insert(*iter);
			}
		}
	}

	int maxNum = 0;
	int zeroNum = 0;
	int oneNum = 0;
	int twoNum = 0;
	int threeNum = 0;
	int fourNum = 0;

	int numbers[9];
	for (int i = 0; i < 9; i++) {
		numbers[i] = 0;
	}
	for (int i = 0; i < edgesNum; i++) {
		numbers[numArr[i].size() / 5]++;
		if (numArr[i].size() > maxNum) {
			maxNum = numArr[i].size();
		}
		if (numArr[i].size() == 0) {
			zeroNum++;
		}
		if (numArr[i].size() == 1) {
			oneNum++;
		}
		if (numArr[i].size() == 2) {
			twoNum++;
		}
		if (numArr[i].size() == 3) {
			threeNum++;
		}
		if (numArr[i].size() == 4) {
			fourNum++;
		}
	}
////	2699 47 153 213 202 (3314) 971 586 455 317 217 134 46 7 2 0

	cout << "The maximal number is: " << maxNum << " " << zeroNum << " "
			<< oneNum << " " << twoNum << " " << threeNum << " " << fourNum
			<< endl;
	for (int i = 0; i < 9; i++) {
		cout << numbers[i] << " ";
	}
//	cout << endl;

	set<int>::iterator iter_s_hat;
	for (iter_s_hat = s_hat.begin(); iter_s_hat != s_hat.end(); iter_s_hat++) {
		edges_flag[*iter_s_hat] = true;
	}

	set<int>::iterator iter_i_hat;
	for (iter_i_hat = i_hat.begin(); iter_i_hat != i_hat.end(); iter_i_hat++) {
		influEdges_flag[*iter_i_hat] = true;
	}

	int size = edgesNum;
//	int size = 500;
//	int size = 1000;

//	double eplison = 0.02;
//	double delta = 0.05;
//	double part1 = 2 * (double) maxNum / (double) edgesNum + eplison;
//	double part2 = eplison * eplison;
//	double part3 = log10(2.0 / delta);
//	int size = part1 / part2 * part3;

//	cout << "hha" << " " << edgesNum << " " << part1 << " " << part2 << " "
//			<< part3 << " " << size << endl;

	set<int>::iterator iter_set;

	set<int> result_seedEdges;
	set<int> result_influEdges;
	set<int> result_influEdgesInSomeBatches = i_hat;
	double result_influScore = 0;
//	cout << "start" << endl;
	while (result_seedEdges.size() < k) {
		int maxEdge;
		int maxInfluNum = 0;
		for (int j = 0; j < size; j++) {
			set<int> temp_edges = r_edge_influEdges[j]; // j: influEdge
			if (influEdges_flag[j] == false) {
//				cout << "1" << endl;
				for (iter_set = temp_edges.begin();
						iter_set != temp_edges.end(); iter_set++) {
					int tempEdge = *iter_set;  // tempEdge: seedEdge
					int tempInfluNum = 0;
//					cout << "2" << endl;
					if (tempEdge < edgesNum) {
						if (edges_flag[tempEdge] == false) { // check whether the current tempEdge has been included in the previous batches
							set<int> influEdges = edge_influEdges[tempEdge];
							set<int>::iterator iter_t;
							for (iter_t = influEdges.begin();
									iter_t != influEdges.end(); iter_t++) {
								int influEdge = *iter_t;
								if (influEdges_flag[influEdge] == false) {
									tempInfluNum++;
								}
							}
							if (tempInfluNum > maxInfluNum) {
								maxEdge = tempEdge;
								maxInfluNum = tempInfluNum;
							}
						}
					}
				}
			}
		}

//		can be improved
		result_seedEdges.insert(maxEdge);
		edges_flag[maxEdge] = true;
//		cout << maxEdge << endl;
		set<int>::iterator iter;
		for (iter = edge_influEdges[maxEdge].begin();
				iter != edge_influEdges[maxEdge].end(); iter++) {
			if (new_edge_density_map[*iter] > lambda * edge_dis_map[maxEdge]
					&& edge_density_map[*iter] >= 0
					&& new_edge_density_map[*iter] > edge_density_map[*iter]) {
				result_influEdges.insert(*iter);
				result_influEdgesInSomeBatches.insert(*iter);
				if (*iter < edgesNum) {
					influEdges_flag[*iter] = true;
				} else {
				}
			}
		}
		result_influScore += maxInfluNum;

	}
	set<int> *result = new set<int> [3];
	result[0] = result_seedEdges;
	result[1] = result_influEdges;
	result[2] = result_influEdgesInSomeBatches;
//	cout << result_influScore << " ; actual: " << result_influEdges.size()
//			<< endl;
	return result;
}

set<int> *greedySelAdaptSamplingEachStep(Graph graph, int k, double lambda,
		map<int, double> edge_dis_map, map<int, set<int>> edge_influEdges,
		map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, set<int> s_hat, set<int> i_hat,
		int size) {
//	cout << "enter" << endl;
	map<int, double> edgeAndInfluNum;
	int edgesNum = edge_influEdges.size();
	map<int, set<int>> r_edge_influEdges;

	bool edges_flag[edgesNum];
	bool influEdges_flag[edgesNum];
//	initialize
	for (int i = 0; i < edgesNum; i++) {
		set<int> r_influEdges;
		r_edge_influEdges[i] = r_influEdges;
	}
	for (int i = 0; i < edgesNum; i++) {
		influEdges_flag[i] = false;
		edges_flag[i] = false;
		set<int> influEdges = edge_influEdges[i];
		set<int>::iterator iter;
		for (iter = influEdges.begin(); iter != influEdges.end(); iter++) {
			if (*iter < edgesNum
					&& new_edge_density_map[*iter]
							> lambda * edge_dis_map[*iter]
					&& edge_density_map[*iter] >= 0
					&& new_edge_density_map[*iter] > edge_density_map[*iter]) {
				r_edge_influEdges[*iter].insert(i);
			}
		}
	}

	set<int>::iterator iter_s_hat;
	for (iter_s_hat = s_hat.begin(); iter_s_hat != s_hat.end(); iter_s_hat++) {
		edges_flag[*iter_s_hat] = true;
	}

	set<int>::iterator iter_i_hat;
	for (iter_i_hat = i_hat.begin(); iter_i_hat != i_hat.end(); iter_i_hat++) {
		influEdges_flag[*iter_i_hat] = true;
	}

//	int size = edgesNum;

//	cout << "hha" << " " << edgesNum << " " << part1 << " " << part2 << " "
//			<< part3 << " " << size << endl;

	set<int>::iterator iter_set;

	set<int> result_seedEdges;
	set<int> result_influEdges;
	set<int> result_influEdgesInSomeBatches = i_hat;
	double result_influScore = 0;
//	cout << "start" << endl;
	while (result_seedEdges.size() < k) {
		int maxEdge;
		int maxInfluNum = 0;
		for (int j = 0; j < size; j++) {
			set<int> temp_edges = r_edge_influEdges[j]; // j: influEdge
			if (influEdges_flag[j] == false) {
				//				cout << "1" << endl;
				for (iter_set = temp_edges.begin();
						iter_set != temp_edges.end(); iter_set++) {
					int tempEdge = *iter_set;  // tempEdge: seedEdge
					int tempInfluNum = 0;
					//					cout << "2" << endl;
					if (tempEdge < edgesNum) {
						if (edges_flag[tempEdge] == false) { // check whether the current tempEdge has been included in the previous batches
							set<int> influEdges = edge_influEdges[tempEdge];
							set<int>::iterator iter_t;
							for (iter_t = influEdges.begin();
									iter_t != influEdges.end(); iter_t++) {
								int influEdge = *iter_t;
								if (influEdges_flag[influEdge] == false) {
									tempInfluNum++;
								}
							}
							if (tempInfluNum > maxInfluNum) {
								maxEdge = tempEdge;
								maxInfluNum = tempInfluNum;
							}
						}
					}
				}
			}
		}

		//		can be improved
		result_seedEdges.insert(maxEdge);
		edges_flag[maxEdge] = true;
//		cout << maxEdge << endl;
		set<int>::iterator iter;
		for (iter = edge_influEdges[maxEdge].begin();
				iter != edge_influEdges[maxEdge].end(); iter++) {
			if (new_edge_density_map[*iter] > lambda * edge_dis_map[maxEdge]
					&& edge_density_map[*iter] >= 0
					&& new_edge_density_map[*iter] > edge_density_map[*iter]) {
				result_influEdges.insert(*iter);
				result_influEdgesInSomeBatches.insert(*iter);
				if (*iter < edgesNum) {
					influEdges_flag[*iter] = true;
				} else {
				}
			}
		}
		result_influScore += maxInfluNum;

	}
	set<int> *result = new set<int> [3];
	result[0] = result_seedEdges;
	result[1] = result_influEdges;
	result[2] = result_influEdgesInSomeBatches;
//	cout << result_influScore << " ; actual: " << result_influEdges.size()
//			<< endl;
	return result;
}

set<int> *greedySelAdaptSampling(Graph graph, int k, double lambda,
		map<int, double> edge_dis_map, map<int, set<int>> edge_influEdges,
		map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, set<int> s_hat, set<int> i_hat) {
	map<int, double> edgeAndInfluNum;
	int edgesNum = edge_influEdges.size();

	map<int, set<int>> r_edge_influEdges;

	bool edges_flag[edgesNum];
	bool influEdges_flag[edgesNum];
//	initialize
	for (int i = 0; i < edgesNum; i++) {
		set<int> r_influEdges;
		r_edge_influEdges[i] = r_influEdges;
	}
	for (int i = 0; i < edgesNum; i++) {
		influEdges_flag[i] = false;
		edges_flag[i] = false;
		set<int> influEdges = edge_influEdges[i];
		set<int>::iterator iter;
		for (iter = influEdges.begin(); iter != influEdges.end(); iter++) {
			if (*iter < edgesNum
					&& new_edge_density_map[*iter]
							> lambda * edge_dis_map[*iter]
					&& edge_density_map[*iter] >= 0
					&& new_edge_density_map[*iter] > edge_density_map[*iter]) {
				r_edge_influEdges[*iter].insert(i);
			}
		}
	}

	int effectEdgesNum = 0;
	for (int i = 0; i < edgesNum; i++) {
		if (r_edge_influEdges[i].size() > 0) {
			effectEdgesNum++;
		}
	}

	set<int>::iterator iter_s_hat;
	for (iter_s_hat = s_hat.begin(); iter_s_hat != s_hat.end(); iter_s_hat++) {
		edges_flag[*iter_s_hat] = true;
	}

	set<int>::iterator iter_i_hat;
	for (iter_i_hat = i_hat.begin(); iter_i_hat != i_hat.end(); iter_i_hat++) {
		influEdges_flag[*iter_i_hat] = true;
	}

	set<int>::iterator iter_set;
	set<int> result_seedEdges;
	set<int> result_influEdges;
	set<int> result_influEdgesInSomeBatches = i_hat;

	set<int> *result = new set<int> [3];
	int size = 50;
	int iteration = 0;
	int deltaSize = pow(2, iteration);
	result = greedySelAdaptSamplingEachStep(graph, k, lambda, edge_dis_map,
			edge_influEdges, edge_density_map, new_edge_density_map, s_hat,
			i_hat, size);
//	cout << "out" << endl;

	result_influEdges = result[1];
	set<int>::iterator influEdges_iter;
	double prop_upper = 0.0;
	for (int i = 0; i < edgesNum; i++) {
		set<int>::iterator result_influEdges_iter = find(
				result_influEdges.begin(), result_influEdges.end(), i);
		double pi = 1
				- pow(
						(1
								- (double) r_edge_influEdges[i].size()
										/ (double) edgesNum), k);
		if (result_influEdges_iter == result_influEdges.end()) {
			prop_upper += 1 - pi;
		} else {
			prop_upper += pi;
		}
	}

//	double probability = prop_upper / (double) edgesNum;
	double probability = prop_upper / (double) effectEdgesNum;
	double eplison = 0.02;
	double delta = 0.05;
	double part1 = 2 + eplison;
	double part2 = eplison * eplison * probability;
	double part3 = log10(2.0 / delta);

	int sizeThreshold = part1 / part2 * part3;
	cout << "sizeThreshod: " << sizeThreshold << " " << probability << " "
			<< part1 << " " << part2 << " " << part3 << endl;
	while (size < sizeThreshold) {
		size = sizeThreshold;
		size += deltaSize;
		iteration++;
		deltaSize = pow(2, iteration);
		result = greedySelAdaptSamplingEachStep(graph, k, lambda, edge_dis_map,
				edge_influEdges, edge_density_map, new_edge_density_map, s_hat,
				i_hat, size);

		result_influEdges = result[1];
		set<int>::iterator influEdges_iter;
		double prop_upper = 0.0;
		for (int i = 0; i < edgesNum; i++) {
			set<int>::iterator result_influEdges_iter = find(
					result_influEdges.begin(), result_influEdges.end(), i);
			double pi = 1
					- pow(
							(1
									- (double) r_edge_influEdges[i].size()
											/ (double) edgesNum), k);
			if (result_influEdges_iter == result_influEdges.end()) {
				prop_upper += 1 - pi;
			} else {
				prop_upper += pi;
			}
		}
		probability = prop_upper / (double) effectEdgesNum;

		part2 = eplison * eplison * probability;
		sizeThreshold = part1 / part2 * part3;
		cout << deltaSize << " " << size << " " << sizeThreshold << endl;
	}
	cout << "actual: " << result_influEdges.size() << endl;
	return result;
//	result[0] = result_seedEdges;
//	result[1] = result_influEdges;
//	result[2] = result_influEdgesInSomeBatches;
//	return result;
}

void computeBound() {
	double alpha = 0.8;
	int k = 20;
	double garma = (1 - pow(alpha, k)) / (1 - alpha) / pow(alpha, k - 1);
	double result = 1 - exp(-k / garma);
	cout << result << endl;

	double result1 = 1 - exp(- alpha);
	cout << result1 << endl;
//	alpha: 0.9; k = 5  -> 0.55
//	alpha: 0.9; k = 20 -> 0.26
//	alpha: 0.8; k = 20 -> 0.056
//	alpha: 0.9; k = 50 -> 0.028

}
