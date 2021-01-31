/*
 * greedySel.h
 *
 *  Created on: 20 May 2020
 *      Author: hui
 */

#ifndef GREEDYSEL_H_
#define GREEDYSEL_H_

#include <iostream>
#include <map>
#include <vector>

#include "graph.h"

using namespace std;

struct union_result {
	set<int> result_influEdges;
	set<int> result_influEdgesInSomeBatches;
	int increInfluScore;
};

int getUnion(int tempEdge, double lambda,
		map<int, double> edge_dis_map, map<int, set<int>> edge_influEdges,
		map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, set<int> result_influEdges,
		set<int> result_influEdgesInSomeBatches);

set<int> *greedySel(Graph graph, int k, double lambda,
		map<int, double> edge_dis_map, map<int, set<int>> edge_influEdges,
		map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, set<int> s_hat, set<int> i_hat);

set<int> *greedySelProgressive(Graph graph, int k, double lambda,
		map<int, double> edge_dis_map, map<int, set<int>> edge_influEdges,
		map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, set<int> s_hat, set<int> i_hat);

//set<int> *greedySel1(Graph graph, int k, double lambda,
//		map<int, double> edge_dis_map, map<int, set<int>> edge_influEdges,
//		map<int, double> edge_density_map,
//		map<int, double> new_edge_density_map, set<int> s_hat, set<int> i_hat);

set<int> *greedySelPruneOne(Graph graph, int k, double lambda,
		map<int, double> edge_dis_map, map<int, set<int>> edge_influEdges,
		map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, set<int> s_hat, set<int> i_hat);

set<int> *greedySelPruneTwo(Graph graph, int k, double lambda,
		map<int, double> edge_dis_map, map<int, set<int>> edge_influEdges,
		map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, set<int> s_hat, set<int> i_hat);

set<int> *greedySelPruneThree(Graph graph, int k, double lambda,
		map<int, double> edge_dis_map, map<int, set<int>> edge_influEdges,
		map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, set<int> s_hat, set<int> i_hat);

vector<int> *greedySelUpperBoundPerBatch(Graph graph, int k, double lambda,
		map<int, double> edge_dis_map, map<int, set<int>> edge_influEdges,
		map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, set<int> s_hat, set<int> i_hat);

set<int> *greedySelSampling(Graph graph, int k, double lambda,
		map<int, double> edge_dis_map, map<int, set<int>> edge_influEdges,
		map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, set<int> s_hat, set<int> i_hat);

set<int> *greedySelAdaptSampling(Graph graph, int k, double lambda,
		map<int, double> edge_dis_map, map<int, set<int>> edge_influEdges,
		map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, set<int> s_hat, set<int> i_hat);

set<int> *greedySelAdaptSamplingEachStep(Graph graph, int k, double lambda,
		map<int, double> edge_dis_map, map<int, set<int>> edge_influEdges,
		map<int, double> edge_density_map,
		map<int, double> new_edge_density_map, set<int> s_hat, set<int> i_hat,
		int size);

void computeBound();

#endif /* GREEDYSEL_H_ */
