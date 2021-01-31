#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <algorithm>
#include <math.h>
#include <queue>

#include "time.h"
#include "edge.h"
#include "graph.h"
#include "fileProcess.h"
#include "branchAndbound.h"
using namespace std;

set<int> *branchAndBoundSel(Graph graph, int k, double lambda,
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
	set<int> edges;
	for (int i = 0; i < edgesNum; i++) {
		edges.insert(i);
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

	set<int> result_seedEdges;
	set<int> result_influEdges;
	set<int> result_influEdgesInSomeBatches = i_hat;
	priority_queue<global_max, vector<global_max>, cmp_max> H;

	int global_lower = 0;
	int global_upper = 9999;

	set<int> selectedSeedEdges;		// I(S)
	set<int> selectedInfluEdges; 	// I(S_i)
	set<int> unscannedEdges = edges; 	// \bar{S_i}
	global_max g;
	g.influScore = edges.size();
	g.seedEdges = selectedSeedEdges;
	g.influEdges = selectedInfluEdges;
	g.unscannedEdges = unscannedEdges;
	H.push(g);
	while (global_lower < global_upper) {
		cout << global_lower << " " << global_upper << " " << H.size() << endl;
		g = H.top();
		H.pop();
		global_upper = g.influScore;
		selectedSeedEdges = g.seedEdges;
		selectedInfluEdges = g.influEdges;
		unscannedEdges = g.unscannedEdges;
		int edge = *unscannedEdges.begin();
		cout << unscannedEdges.size() << endl;
		unscannedEdges.erase(edge);

		if (selectedSeedEdges.size() < k) {
			cout << "enter: old; " << selectedSeedEdges.size() << " "
					<< unscannedEdges.size() << endl;
//			for the raw branch
			lower_upper_bound old_lower_bound = computeLowerBound(
					selectedSeedEdges, unscannedEdges, lambda, edge_dis_map,
					edge_influEdges, edge_density_map, new_edge_density_map,
					selectedInfluEdges, k);
			cout << "hhah1" << endl;

			int old_upper_bound = computeUpperBound1(selectedSeedEdges,
					unscannedEdges, lambda, edge_dis_map, edge_influEdges,
					edge_density_map, new_edge_density_map, selectedInfluEdges,
					k);
			cout << "hhah1 out" << endl;

			cout << "old_lower_bound.influScore_k: "
					<< old_lower_bound.influScore_k << " "
					<< old_lower_bound.influEdges.size() << endl;
			if (old_lower_bound.influScore_k > global_lower) {
				global_lower = old_lower_bound.influScore_k;
				result_seedEdges = old_lower_bound.seedEdges_k;
				result_influEdges = old_lower_bound.influEdges;
			}
			cout << "old_upper_bound: " << old_upper_bound << endl;
			if (old_upper_bound > global_lower) {
				global_max old_g;
				old_g.influScore = old_upper_bound;
				old_g.influEdges = selectedInfluEdges;
				old_g.seedEdges = selectedSeedEdges;
				old_g.unscannedEdges = unscannedEdges;
				H.push(old_g);
			}

//			for the branch adding a new element
			cout << "enter: new; " << selectedSeedEdges.size() << " "
					<< unscannedEdges.size() << endl;
			selectedSeedEdges.insert(edge);
			set<int>::iterator iter;
			for (iter = edge_influEdges[edge].begin();
					iter != edge_influEdges[edge].end(); iter++) {
				if (new_edge_density_map[*iter] > lambda * edge_dis_map[edge]
						&& edge_density_map[*iter] >= 0
						&& new_edge_density_map[*iter]
								> edge_density_map[*iter]) {
					selectedInfluEdges.insert(*iter);
				}
			}
			lower_upper_bound new_lower_bound = computeLowerBound(
					selectedSeedEdges, unscannedEdges, lambda, edge_dis_map,
					edge_influEdges, edge_density_map, new_edge_density_map,
					selectedInfluEdges, k);
			int new_upper_bound = computeUpperBound1(selectedSeedEdges,
					unscannedEdges, lambda, edge_dis_map, edge_influEdges,
					edge_density_map, new_edge_density_map, selectedInfluEdges,
					k);

			if (new_lower_bound.influScore_k > global_lower) {
				global_lower = new_lower_bound.influScore_k;
				result_seedEdges = new_lower_bound.seedEdges_k;
				result_influEdges = new_lower_bound.influEdges;
			}
			cout << "new_upper_bound: " << new_upper_bound << endl;
			if (new_upper_bound > global_lower) {
				global_max new_g;
				new_g.influScore = new_upper_bound;
				new_g.influEdges = selectedInfluEdges;
				new_g.seedEdges = selectedSeedEdges;
				new_g.unscannedEdges = unscannedEdges;
				H.push(new_g);
			}
		}
	}
	set<int> *result = new set<int> [3];
	result[0] = result_seedEdges;
	result[1] = result_influEdges; // the influenced edges in the i-th batches
	result[2] = result_influEdgesInSomeBatches; // the influenced edges in the first i batches
	cout << "actual: " << result_influEdges.size() << endl;
	return result;
}

lower_upper_bound computeLowerBound(set<int> seedEdges, set<int> unscannedEdges,
		double lambda, map<int, double> edge_dis_map,
		map<int, set<int>> edge_influEdges, map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, set<int> influEdges, int k) {
	lower_upper_bound bound;
	set<int> result_seedEdges = seedEdges;
	set<int> result_influEdges = influEdges;
	int result_influScore = 0;
	set<int>::iterator iter_set;

	while (result_seedEdges.size() < k) {
		int maxEdge;
		int maxInfluScore = 0;
		for (iter_set = unscannedEdges.begin();
				iter_set != unscannedEdges.end(); iter_set++) {
			int tempEdge = *iter_set;
			int tempInfluScore = 0;
			if (edge_influEdges[tempEdge].size() != 0) {
				set<int>::iterator iter_set;
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
//						if (influEdges_flag[*iter_set] == false) {
							tempInfluScore++;
						}
					}
				}
				if (tempInfluScore > maxInfluScore) {
					maxEdge = tempEdge;
					maxInfluScore = tempInfluScore;
				}
			}
		}
//		can be improved
		result_seedEdges.insert(maxEdge);
		unscannedEdges.erase(maxEdge);
		result_influScore += maxInfluScore;
//		edges_flag[maxEdge] = true;
		set<int>::iterator iter;
		for (iter = edge_influEdges[maxEdge].begin();
				iter != edge_influEdges[maxEdge].end(); iter++) {
			if (new_edge_density_map[*iter] > lambda * edge_dis_map[*iter]
					&& edge_density_map[*iter] >= 0
					&& new_edge_density_map[*iter] > edge_density_map[*iter]) {
				result_influEdges.insert(*iter);
			}
		}
	}
	bound.influScore_k = result_influScore;
	bound.seedEdges_k = result_seedEdges;
	bound.influEdges = result_influEdges;
//	cout << result_influEdges.size() << " " << result_influScore << endl;
	return bound;
}

int computeUpperBound0(set<int> seedEdges, set<int> unscannedEdges,
		double lambda, map<int, double> edge_dis_map,
		map<int, set<int>> edge_influEdges, map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, set<int> influEdges, int k) {
	set<int> result_seedEdges = seedEdges;
	set<int> result_influEdges = influEdges;
	int result_influScore = 0;
	set<int>::iterator iter_set;

	int maxEdge;
	int maxInfluScore = 0;
	for (iter_set = unscannedEdges.begin(); iter_set != unscannedEdges.end();
			iter_set++) {
		int tempEdge = *iter_set;
		int tempInfluScore = 0;
		if (edge_influEdges[tempEdge].size() != 0) {
			set<int>::iterator iter_set;
			for (iter_set = edge_influEdges[tempEdge].begin();
					iter_set != edge_influEdges[tempEdge].end(); iter_set++) {
				if (new_edge_density_map[*iter_set]
						> lambda * edge_dis_map[*iter_set]
						&& edge_density_map[*iter_set] >= 0
						&& new_edge_density_map[*iter_set]
								> edge_density_map[*iter_set]) {
					set<int>::iterator iter_vec = find(
							result_influEdges.begin(), result_influEdges.end(),
							*iter_set); // iterator is used to save the idx
					if (iter_vec == result_influEdges.end()) {
						tempInfluScore++;
					}
				}
			}
			if (tempInfluScore > maxInfluScore) {
				maxEdge = tempEdge;
				maxInfluScore = tempInfluScore;
			}
		}
	}
	result_influScore = result_influEdges.size()
			+ maxInfluScore * (k - result_seedEdges.size());
	return result_influScore;
}

int computeUpperBound1(set<int> seedEdges, set<int> unscannedEdges,
		double lambda, map<int, double> edge_dis_map,
		map<int, set<int>> edge_influEdges, map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, set<int> influEdges, int k) {
	set<int> result_seedEdges = seedEdges;
	set<int> result_influEdges = influEdges;
	int result_influScore = 0;
	set<int>::iterator iter_set;

	priority_queue<global_max, vector<global_max>, cmp_max> H;

	for (iter_set = unscannedEdges.begin(); iter_set != unscannedEdges.end();
			iter_set++) {
		int tempEdge = *iter_set;
		int tempInfluScore = 0;
		global_max new_g;
		if (edge_influEdges[tempEdge].size() != 0) {
			set<int>::iterator iter_set;
			for (iter_set = edge_influEdges[tempEdge].begin();
					iter_set != edge_influEdges[tempEdge].end(); iter_set++) {
				if (new_edge_density_map[*iter_set]
						> lambda * edge_dis_map[*iter_set]
						&& edge_density_map[*iter_set] >= 0
						&& new_edge_density_map[*iter_set]
								> edge_density_map[*iter_set]) {
					set<int>::iterator iter_vec = find(
							result_influEdges.begin(), result_influEdges.end(),
							*iter_set); // iterator is used to save the idx
					if (iter_vec == result_influEdges.end()) {
						tempInfluScore++;
					}
				}
			}
			new_g.influScore = tempInfluScore;
		} else {
			new_g.influScore = 0;
		}
		new_g.seedEdges.insert(tempEdge);
		H.push(new_g);
	}
	int i = 0;
	while (i < k - result_seedEdges.size()) {
//		cout << "here 2" << endl;
		global_max g = H.top();
//		cout << "here 1" << endl;
		H.pop();
//		cout << "here" << endl;
		result_influScore += result_influEdges.size() + g.influScore;
		i++;
	}
	return result_influScore;
}

vector<vector<int>> ha;
vector<int> v;
void hehe(int p, int dep, int k, int n) {
	if (dep == k) {
		ha.push_back(v);
		return;
	}

	for (int i = p; i <= n; ++i) {
		v[dep] = i;
		hehe(i + 1, dep + 1, k, n);
	}

}
vector<vector<int>> combine(int n, int k) {
	ha.clear();
	v.clear();
	v.resize(k);
	hehe(1, 0, k, n);
	return ha;
}

void exact(Graph graph, int k, double lambda, map<int, double> edge_dis_map,
		map<int, set<int>> edge_influEdges, map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, set<int> s_hat, set<int> i_hat) {
	int edgesNum = edge_influEdges.size();
	combine(edgesNum, k);

	set<int> result_seedEdges;
	set<int> result_influEdges;
	int result_influScore = 0;
	int len1 = ha.size();
	cout << "" << len1 << endl;
	for (int i = 0; i < len1; i++) {
		int len2 = ha[i].size();
		set<int> temp_result_seedEdges;
		set<int> temp_result_influEdges;
		int temp_result_influScore = 0;
		for (int j = 0; j < len2; j++) {
			int tempEdge = ha[i][j];
			temp_result_seedEdges.insert(tempEdge);
			if (edge_influEdges[tempEdge].size() != 0) {
				set<int>::iterator iter_set;
				for (iter_set = edge_influEdges[tempEdge].begin();
						iter_set != edge_influEdges[tempEdge].end();
						iter_set++) {
					if (new_edge_density_map[*iter_set]
							> lambda * edge_dis_map[*iter_set]
							&& edge_density_map[*iter_set] >= 0
							&& new_edge_density_map[*iter_set]
									> edge_density_map[*iter_set]) {
						set<int>::iterator iter_vec = find(
								temp_result_influEdges.begin(),
								temp_result_influEdges.end(), *iter_set); // iterator is used to save the idx
						if (iter_vec == temp_result_influEdges.end()) {
							temp_result_influScore++;
							temp_result_influEdges.insert(*iter_vec);
						}
					}
				}
			}
		}
		cout << i << " " << temp_result_influScore << " " << result_influScore
				<< endl;
		if (temp_result_influScore > result_influScore) {
			result_seedEdges = temp_result_seedEdges;
			result_influEdges = temp_result_influEdges;
			result_influScore = temp_result_influScore;
		}
	}
	cout << result_influScore << endl;
	set<int>::iterator iter_seed;
	for (iter_seed = result_seedEdges.begin();
			iter_seed != result_seedEdges.end(); iter_seed++) {
		cout << *iter_seed << " ";
	}
}
