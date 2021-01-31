#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <algorithm>

#include "time.h"
#include "edge.h"
#include "graph.h"
#include "fileProcess.h"
#include "globalSel.h"
#include "greedySel.h"
using namespace std;

set<int> globalSel(int algorithm, Graph graph, int K, int batchSize,
		double lambda, map<int, double> edge_dis_map,
		map<int, set<int>> edge_influEdges, map<int, double> edge_density_map,
		map<string, double> edge_pro_map, map<int, double> edge_residual_map,
		map<int, set<Path>> edge_influPaths, string output_file) {
	clock_t time0, time1, time2, time3, time4;
	set<int> result_seedEdges;
	set<int> result_influEdges;
//	map<int, double> edgeAndInfluScore;
	map<int, int> edgeAndInfluNum;
	int edgesNum = edge_influEdges.size();
	map<int, double> new_edge_density_map; // which will be updated after each batch

//	S \xi
	set<int> seedEdgesInOneBatch[batchSize + 1][K + 1];
//	bool seedEdgesInOneBatch_flag[batchSize + 1][K + 1][edgesNum];
//	\Phi(S) \mathbb{I}
	set<int> seedEdgesInSomeBatch[batchSize + 1][K + 1];
	bool seedEdgesInSomeBatch_flag[batchSize + 1][K + 1][edgesNum];

	set<int> influEdgesInOneBatch[batchSize + 1][K + 1];
//	bool influEdgesInOneBatch_flag[batchSize + 1][K + 1][edgesNum];
	set<int> influEdgesInSomeBatch[batchSize + 1][K + 1];
	bool influEdgesInSomeBatch_flag[batchSize + 1][K + 1][edgesNum];

	for (int i = 1; i <= batchSize; i++) {
		for (int k = 1; k <= K; k++) {
//			wrong
			map<int, set<Path>> edge_influEdges_replace;
			new_edge_density_map = graph.getEdgeDensity_afterDiffusion(
					edge_density_map, edge_pro_map, edge_residual_map,
					edge_influEdges_replace);
			set<int> s_hat = seedEdgesInSomeBatch[i - 1][K - k];
			set<int> i_hat = influEdgesInSomeBatch[i - 1][K - k];
//			bool *s_hat_flag = new bool[edgesNum];
//			bool *i_hat_flag = new bool[edgesNum];
//			s_hat_flag = seedEdgesInSomeBatch_flag[i - 1][K - k];
//			i_hat_flag = influEdgesInSomeBatch_flag[i - 1][K - k];
			set<int> *resultByBatchSel;
			algorithm = 5;
			if (algorithm == 0) {
				resultByBatchSel = greedySel(graph, k, lambda, edge_dis_map,
						edge_influEdges, edge_density_map, new_edge_density_map,
						s_hat, i_hat);
			} else if (algorithm == 1) {
				resultByBatchSel = greedySelPruneOne(graph, k, lambda,
						edge_dis_map, edge_influEdges, edge_density_map,
						new_edge_density_map, s_hat, i_hat);
			} else if (algorithm == 2) {
				resultByBatchSel = greedySelPruneTwo(graph, k, lambda,
						edge_dis_map, edge_influEdges, edge_density_map,
						new_edge_density_map, s_hat, i_hat);
			} else if (algorithm == 3) {
				resultByBatchSel = greedySelSampling(graph, k, lambda,
						edge_dis_map, edge_influEdges, edge_density_map,
						new_edge_density_map, s_hat, i_hat);
			} else if (algorithm == 4) {
				resultByBatchSel = greedySelAdaptSampling(graph, k, lambda,
						edge_dis_map, edge_influEdges, edge_density_map,
						new_edge_density_map, s_hat, i_hat);
			} else if (algorithm == 5) {
				resultByBatchSel = greedySelProgressive(graph, k, lambda,
						edge_dis_map, edge_influEdges, edge_density_map,
						new_edge_density_map, s_hat, i_hat);
			}
			seedEdgesInOneBatch[i][k] = resultByBatchSel[0];
			influEdgesInOneBatch[i][k] = resultByBatchSel[1];
			cout << "i: " << i << ";k: " << k << "; seedNum: "
					<< seedEdgesInOneBatch[i][k].size() << "; influNum: "
					<< influEdgesInOneBatch[i][k].size() << endl;
			int max_influNum = 0;
			int max_j = 0;
			for (int j = 0; j <= k; j++) {
				int temp_influNum = influEdgesInSomeBatch[i - 1][k - j].size()
						+ influEdgesInOneBatch[i][j].size();
//				cout << temp_influNum << " " << max_influNum << " "
//						<< influEdgesInSomeBatch[i - 1][k - j].size() << " "
//						<< influEdgesInOneBatch[i][j].size() << endl;
				if (temp_influNum > max_influNum) {
					max_influNum = temp_influNum;
					max_j = j;
				}
			}
			set<int> union_seedEdges;
			set<int> sEdges1 = seedEdgesInSomeBatch[i - 1][k - max_j];
			set<int> sEdges2 = seedEdgesInOneBatch[i][max_j];
			set_union(sEdges1.begin(), sEdges1.end(), sEdges2.begin(),
					sEdges2.end(),
					inserter(union_seedEdges, union_seedEdges.begin()));
			seedEdgesInSomeBatch[i][k] = union_seedEdges;

			set<int> union_influEdges;
			set<int> iEdges1 = influEdgesInSomeBatch[i - 1][k - max_j];
			set<int> iEdges2 = influEdgesInOneBatch[i][max_j];
			set_union(iEdges1.begin(), iEdges1.end(), iEdges2.begin(),
					iEdges2.end(),
					inserter(union_influEdges, union_influEdges.begin()));
			influEdgesInSomeBatch[i][k] = union_influEdges;
			edge_density_map = new_edge_density_map;
		}
	}
	result_seedEdges = seedEdgesInSomeBatch[batchSize][K];
	result_influEdges = influEdgesInSomeBatch[batchSize][K];
	cout << "The result is: " << result_seedEdges.size() << " "
			<< result_influEdges.size() << endl;
	cout << "The result is: " << result_seedEdges.size() << " "
			<< result_influEdges.size() << endl;
	string runtime_g = "GlobalSel: " + to_string(result_seedEdges.size()) + " "
			+ to_string(result_influEdges.size());
	writeOutputFile(output_file, runtime_g);
	return result_seedEdges;
}

set<int> globalSelEven(int algorithm, Graph graph, int K, int batchSize,
		double lambda, map<int, double> edge_dis_map,
		map<int, set<int>> edge_influEdges, map<int, double> edge_density_map,
		map<string, double> edge_pro_map, map<int, double> edge_residual_map,
		map<int, set<Path>> edge_influPaths, string output_file) {
	clock_t time0, time1, time2, time3, time4;
	set<int> result_seedEdges;
	set<int> result_influEdges;
	//	map<int, double> edgeAndInfluScore;
	map<int, int> edgeAndInfluNum;
	int edgesNum = edge_influEdges.size();
	map<int, double> new_edge_density_map; // which will be updated after each batch

	int k = 0;
	if (K % batchSize == 0) {
		k = K / batchSize;
	} else {
		k = K / batchSize + 1;
		batchSize = K / k + 1;
	}
	//	S \xi
	set<int> seedEdgesInOneBatch[batchSize + 1];
	//	\Phi(S) \mathbb{I}
	set<int> seedEdgesInSomeBatch[batchSize + 1];

	set<int> influEdgesInOneBatch[batchSize + 1];
	set<int> influEdgesInSomeBatch[batchSize + 1];

	for (int i = 1; i < batchSize + 1; i++) {
		//			wrong
		map<int, set<Path>> edge_influEdges_replace;
		new_edge_density_map = graph.getEdgeDensity_afterDiffusion(
				edge_density_map, edge_pro_map, edge_residual_map,
				edge_influEdges_replace);
		set<int> s_hat = seedEdgesInSomeBatch[i - 1];
		set<int> i_hat = influEdgesInSomeBatch[i - 1];
		set<int> *resultByBatchSel;
		if (i == batchSize) {
			k = K - k * (i - 1); // the number of remaining edges in the final batch
		}
		if (algorithm == 0) {
			resultByBatchSel = greedySel(graph, k, lambda, edge_dis_map,
					edge_influEdges, edge_density_map, new_edge_density_map,
					s_hat, i_hat);
		} else if (algorithm == 1) {
			resultByBatchSel = greedySelPruneOne(graph, k, lambda, edge_dis_map,
					edge_influEdges, edge_density_map, new_edge_density_map,
					s_hat, i_hat);
		} else if (algorithm == 2) {
			resultByBatchSel = greedySelPruneTwo(graph, k, lambda, edge_dis_map,
					edge_influEdges, edge_density_map, new_edge_density_map,
					s_hat, i_hat);
		}
		seedEdgesInOneBatch[i] = resultByBatchSel[0];
		influEdgesInOneBatch[i] = resultByBatchSel[1];
		influEdgesInSomeBatch[i] = resultByBatchSel[2];
		cout << "i: " << i << ";k: " << k << "; seedNum: "
				<< seedEdgesInOneBatch[i].size() << "; influNum: "
				<< influEdgesInOneBatch[i].size() << endl;

		set<int> union_seedEdges;
		set<int> sEdges1 = seedEdgesInSomeBatch[i - 1];
		set<int> sEdges2 = seedEdgesInOneBatch[i];
		set_union(sEdges1.begin(), sEdges1.end(), sEdges2.begin(),
				sEdges2.end(),
				inserter(union_seedEdges, union_seedEdges.begin()));
		seedEdgesInSomeBatch[i] = union_seedEdges;

//		set<int> union_influEdges;
//		set<int> iEdges1 = influEdgesInSomeBatch[i - 1];
//		set<int> iEdges2 = influEdgesInOneBatch[i];
//		set_union(iEdges1.begin(), iEdges1.end(), iEdges2.begin(),
//				iEdges2.end(),
//				inserter(union_influEdges, union_influEdges.begin()));
//		influEdgesInSomeBatch[i] = union_influEdges;
		edge_density_map = new_edge_density_map;
	}
	result_seedEdges = seedEdgesInSomeBatch[batchSize];
	result_influEdges = influEdgesInSomeBatch[batchSize];
	cout << "The result is: " << result_seedEdges.size() << " "
			<< result_influEdges.size() << endl;
	string runtime_g = "GlobalSelPruneOne: "
			+ to_string(result_seedEdges.size()) + " "
			+ to_string(result_influEdges.size());
	writeOutputFile(output_file, runtime_g);
	return result_seedEdges;
}

set<int> globalSelPruneOne(int algorithm, Graph graph, int K, int batchSize,
		double lambda, map<int, double> edge_dis_map,
		map<int, set<int>> edge_influEdges, map<int, double> edge_density_map,
		map<string, double> edge_pro_map, map<int, double> edge_residual_map,
		map<int, set<Path>> edge_influPaths, string output_file) {
	clock_t time0, time1, time2, time3, time4;
	set<int> result_seedEdges;
	set<int> result_influEdges;
	map<int, int> edgeAndInfluNum;
	int edgesNum = edge_influEdges.size();
	map<int, double> new_edge_density_map; // which will be updated after each batch

//	S \xi
	set<int> seedEdgesInOneBatch[batchSize + 1][K + 1];
	set<int> seedEdgesInSomeBatch[batchSize + 1][K + 1];

//	\Phi(S) \mathbb{I}
	set<int> influEdgesInOneBatch[batchSize + 1][K + 1];
	vector<int> influEdgesInOneBatch_upper[batchSize + 1][K + 1];
	set<int> influEdgesInSomeBatch[batchSize + 1][K + 1];
	set<int> influEdgesInSomeBatch_lower[batchSize + 1][K + 1];

	for (int i = 1; i <= batchSize; i++) {
		for (int k = 1; k <= K; k++) {
			//			wrong
			map<int, set<Path>> edge_influEdges_replace;
			new_edge_density_map = graph.getEdgeDensity_afterDiffusion(
					edge_density_map, edge_pro_map, edge_residual_map,
					edge_influEdges_replace);
			influEdgesInSomeBatch_lower[i][k] = influEdgesInSomeBatch[i - 1][k];

			for (int j = 0; j <= k; j++) {
				vector<int> *resultByBatchSel_upper;
				set<int> s_hat = seedEdgesInSomeBatch[i - 1][K - k];
				set<int> i_hat = influEdgesInSomeBatch[i - 1][K - k];
				resultByBatchSel_upper = greedySelUpperBoundPerBatch(graph, k,
						lambda, edge_dis_map, edge_influEdges, edge_density_map,
						new_edge_density_map, s_hat, i_hat);
//				EstimateBound
				influEdgesInOneBatch_upper[i][j] = resultByBatchSel_upper[1];
				int temp_influNum_upper =
						influEdgesInSomeBatch[i - 1][k - j].size()
								+ influEdgesInOneBatch_upper[i][j].size();
				if (influEdgesInSomeBatch_lower[i][j].size()
						<= temp_influNum_upper) {
					set<int> *resultByBatchSel;
					if (algorithm == 0) {
						resultByBatchSel = greedySel(graph, k, lambda,
								edge_dis_map, edge_influEdges, edge_density_map,
								new_edge_density_map, s_hat, i_hat);
					} else if (algorithm == 1) {
						resultByBatchSel = greedySelPruneOne(graph, k, lambda,
								edge_dis_map, edge_influEdges, edge_density_map,
								new_edge_density_map, s_hat, i_hat);
					} else if (algorithm == 2) {
						resultByBatchSel = greedySelPruneTwo(graph, k, lambda,
								edge_dis_map, edge_influEdges, edge_density_map,
								new_edge_density_map, s_hat, i_hat);
					}
					seedEdgesInOneBatch[i][k] = resultByBatchSel[0];
					influEdgesInOneBatch[i][k] = resultByBatchSel[1];
					cout << "i: " << i << ";k: " << k << "; seedNum: "
							<< seedEdgesInOneBatch[i][k].size()
							<< "; influNum: "
							<< influEdgesInOneBatch[i][k].size() << endl;

					set<int> union_seedEdges;
					set<int> sEdges1 = seedEdgesInSomeBatch[i - 1][k - j];
					set<int> sEdges2 = seedEdgesInOneBatch[i][j];
					set_union(sEdges1.begin(), sEdges1.end(), sEdges2.begin(),
							sEdges2.end(),
							inserter(union_seedEdges, union_seedEdges.begin()));
					seedEdgesInSomeBatch[i][k] = union_seedEdges;

					set<int> union_influEdges;
					set<int> iEdges1 = influEdgesInSomeBatch[i - 1][k - j];
					set<int> iEdges2 = influEdgesInOneBatch[i][j];
					set_union(iEdges1.begin(), iEdges1.end(), iEdges2.begin(),
							iEdges2.end(),
							inserter(union_influEdges,
									union_influEdges.begin()));
					influEdgesInSomeBatch_lower[i][k] = union_influEdges;
				}
			}
			influEdgesInSomeBatch[i][k] = influEdgesInSomeBatch_lower[i][k];
		}
		edge_density_map = new_edge_density_map;
	}
	result_seedEdges = seedEdgesInSomeBatch[batchSize][K];
	result_influEdges = influEdgesInSomeBatch[batchSize][K];
	cout << "The result is: " << result_seedEdges.size() << " "
			<< result_influEdges.size() << endl;
	return result_seedEdges;
}
