/*
 * globalSel.h
 *
 *  Created on: 20 May 2020
 *      Author: hui
 */

#ifndef GLOBALSEL_H_
#define GLOBALSEL_H_

#include <iostream>
#include <map>
#include <vector>

#include "graph.h"

using namespace std;

set<int> globalSel(int algorithm, Graph graph, int K, int batchSize,
		double lambda, map<int, double> edge_dis_map,
		map<int, set<int>> edge_influEdges, map<int, double> edge_density_map,
		map<string, double> edge_pro_map, map<int, double> edge_residual_map,
		map<int, set<Path>> edge_influPaths, string output_file);

set<int> globalSelEven(int algorithm, Graph graph, int K, int batchSize,
		double lambda, map<int, double> edge_dis_map,
		map<int, set<int>> edge_influEdges, map<int, double> edge_density_map,
		map<string, double> edge_pro_map, map<int, double> edge_residual_map, map<int, set<Path>> edge_influPaths,
		string output_file);

set<int> globalSelPruneOne(int algorithm, Graph graph, int K, int batchSize,
		double lambda, map<int, double> edge_dis_map,
		map<int, set<int>> edge_influEdges, map<int, double> edge_density_map,
		map<string, double> edge_pro_map, map<int, double> edge_residual_map, map<int, set<Path>> edge_influPaths,
		string output_file);

#endif /* GLOBALSEL_H_ */
