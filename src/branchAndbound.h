/*
 * branchAndbound.h
 *
 *  Created on: 24 Jul 2020
 *      Author: hui
 */

#ifndef BRANCHANDBOUND_H_
#define BRANCHANDBOUND_H_

struct global_max {
	set<int> seedEdges;
	set<int> influEdges;
	set<int> unscannedEdges;
	int influScore;
};
struct cmp_max {
	bool operator()(const global_max &a, const global_max &b) {
		return a.influScore <= b.influScore;
	}
};

struct lower_upper_bound {
	set<int> seedEdges_k;
	int influScore_k;
	set<int> influEdges;
};

set<int> *branchAndBoundSel(Graph graph, int k, double lambda,
		map<int, double> edge_dis_map, map<int, set<int>> edge_influEdges,
		map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, set<int> s_hat, set<int> i_hat);

lower_upper_bound computeLowerBound(set<int> seedEdges, set<int> unscannedEdges,
		double lambda, map<int, double> edge_dis_map,
		map<int, set<int>> edge_influEdges, map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, set<int> influEdges, int k);

// to do
int computeUpperBound0(set<int> seedEdges, set<int> unscannedEdges,
		double lambda, map<int, double> edge_dis_map,
		map<int, set<int>> edge_influEdges, map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, set<int> influEdges, int k);

// the upper bound is defined as the sum of the largest k-|S| marginal gain
int computeUpperBound1(set<int> seedEdges, set<int> unscannedEdges,
		double lambda, map<int, double> edge_dis_map,
		map<int, set<int>> edge_influEdges, map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, set<int> influEdges, int k);

void exact(Graph graph, int k, double lambda, map<int, double> edge_dis_map,
		map<int, set<int>> edge_influEdges, map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, set<int> s_hat, set<int> i_hat);

#endif /* BRANCHANDBOUND_H_ */
