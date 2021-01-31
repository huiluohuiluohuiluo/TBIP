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
#include "duration.h"
#include "partition.h"

using namespace std;
map<int, double> edgeAndEstScore_d;
class Cmpare {
public:
	bool operator()(int a, int b) {
		return edgeAndEstScore_d[a] > edgeAndEstScore_d[b];
	}
};

// edge_influEdges: the influence to the in-neighbor edges (to compute the influence)
// r_edge_influEdges: the influence to the out-neighbor edges (to compute the traffic spread)
set<int> *durationSel(Graph graph, int K, int M, int batchSize,
		int durationSize, double lambda, map<int, double> edge_dis_map,
		map<int, set<int>> edge_influEdges,
		map<int, set<int>> r_edge_influEdges, map<int, double> edge_density_map,
		map<int, double> new_edge_density_map,
		map<int, map<int, double>> edge_pro_map,
		map<int, double> edge_residual_map, string input_partition_file,
		string output_file, string parType, int pruneFlag, int algFlag,
		string city, string input_partition_minCut_file,
		map<int, Edge> edges_map, string input_influEdges_file) {
	clock_t time0, time1, time2, time3, time4, time01, time02;
	ifstream file(input_influEdges_file);
	int edgesNum = edge_influEdges.size();
	map<int, set<int>> edgesAndInfluEdges;
	double duration01 = 0.0;
	double duration0 = 0.0;
	double duration1 = 0.0;
	cout << "he" << endl;
	if (!file) {
		cout << "h" << endl;
		ofstream outfile_edges(input_influEdges_file, ios::app);
		int slides = batchSize - durationSize + 1;
		cout << "The number of slides is: " << slides << endl;

		time0 = clock();
		//	get the influence edges in all batches
		map<int, map<int, set<int>>> edgesWithInfluInAllBatches;
		//	i: batchSize(1); j: edgeId
		//	int sum = 0;
		for (int i = 1; i <= batchSize; i++) { // batches
			map<int, set<int>> a_edge_influEdges;
			for (int j = 0; j < edgesNum; j++) {	// edges
				set<int> a_influEdges;
				set<int>::iterator iter_set;
				for (iter_set = edge_influEdges[j].begin();
						iter_set != edge_influEdges[j].end(); iter_set++) {	// influEdges
					if (new_edge_density_map[*iter_set]
							> lambda * edge_dis_map[*iter_set]
							&& edge_density_map[*iter_set] >= 0) {
						a_influEdges.insert(*iter_set);
					}
				}
				a_edge_influEdges[j] = a_influEdges;
			}
			edgesWithInfluInAllBatches[i] = a_edge_influEdges;

			//		update the traffic volume
			time01 = clock();
			edge_density_map = new_edge_density_map;
			double sum_density = 0.0;
			for (int ii = 0; ii < edgesNum; ii++) {
				sum_density += edge_density_map[ii];
			}
			cout << "The batch size is: " << i << "; the sum of density: "
					<< sum_density << endl;
			new_edge_density_map = graph.getEdgeDensity_afterDiffusion_Pro(
					edge_density_map, edge_pro_map, edge_residual_map,
					r_edge_influEdges);

		}
		time1 = clock();
		//	int congestedNum = sum / batchSize;
		//	cout << "finish 1" << endl;
		//	cout << "The number of congested roads is: " << congestedNum << endl;

		//	consider duration factor and filter some in edgesWithInfluInAllBatches
		//	int pruneFlag = 0;
		//	0: without prune; 1: prune
		if (pruneFlag == 0) {
			for (int i = 1; i <= slides; i++) { // batches
				map<int, set<int>> edge_influEdges =
						edgesWithInfluInAllBatches[i];
				for (int j = 0; j < edgesNum; j++) {
					set<int> intersection_Edges =
							edgesWithInfluInAllBatches[i][j];
					for (int k = 1; k < durationSize; k++) {
						set<int> sEdges1 = intersection_Edges;
						set<int> sEdges2 = edgesWithInfluInAllBatches[i + k][j];
						set_intersection(sEdges1.begin(), sEdges1.end(),
								sEdges2.begin(), sEdges2.end(),
								inserter(intersection_Edges,
										intersection_Edges.begin()));
					}
					//					edgesAndInfluEdges[j] = intersection_Edges;
					set<int>::iterator iter_set;
					for (iter_set = intersection_Edges.begin();
							iter_set != intersection_Edges.end(); iter_set++) {
						edgesAndInfluEdges[j].insert(*iter_set);
					}
				}
			}

			for (int j = 0; j < edgesNum; j++) {
				string line = to_string(j) + ";";
				set<int>::iterator iter_set;
				for (iter_set = edgesAndInfluEdges[j].begin();
						iter_set != edgesAndInfluEdges[j].end(); iter_set++) {
					line += to_string(*iter_set) + ",";
				}
				line.pop_back();
				outfile_edges << line << endl;
			}
		} else if (pruneFlag == 1) {
			cout << "enter Prunning ..." << endl;
			//		prunning
			//		i: batchSize(1); j: edgeId; k: duration(0)
			//		map<int, map<int, set<int>>> influEdgesIntersection;
			//		influEdgesIntersection: \ksi; edgesWithInfluInAllBatches: A
			map<int, map<int, map<int, set<int>>> > startAndEndInfluEdges;
			//		startAndEndInfluEdges: B; i: edgeId; j: startSlide (include, starts from startSlide-1); k: endSlide (include);
			for (int j = 0; j < edgesNum; j++) {
				set<int> intersection_Edges;
				for (int k = durationSize; k > 0; k--) {
					set<int> sEdges1 = intersection_Edges;
					set<int> sEdges2 = edgesWithInfluInAllBatches[k][j];
					set_intersection(sEdges1.begin(), sEdges1.end(),
							sEdges2.begin(), sEdges2.end(),
							inserter(intersection_Edges,
									intersection_Edges.begin()));
					startAndEndInfluEdges[j][k][durationSize] =
							intersection_Edges;
				}
				//			startAndEndInfluEdges[2][j] = intersection_Edges;
				////			connect with the first slide
				//			set<int> sEdges1 = edgesWithInfluInAllBatches[1][j];
				//			set<int> sEdges2 = intersection_Edges;
				//			set_intersection(sEdges1.begin(), sEdges1.end(),
				//					sEdges2.begin(), sEdges2.end(),
				//					inserter(intersection_Edges,
				//							intersection_Edges.begin()));
				//			influEdgesIntersection[1][j] = intersection_Edges;
				//			edgesAndInfluEdges[j] = intersection_Edges;
			}

			// i: starting batch; j: edgeId; k: duration size;
			for (int i = 2; i <= slides; i++) { // batches
				for (int j = 0; j < edgesNum; j++) {
					set<int> intersection_Edges;
					set<int> iEdges1 = startAndEndInfluEdges[j][i][i
							+ durationSize - 2];
					set<int> iEdges2 = edgesWithInfluInAllBatches[i
							+ durationSize - 1][j];
					set_intersection(iEdges1.begin(), iEdges1.end(),
							iEdges2.begin(), iEdges2.end(),
							inserter(intersection_Edges,
									intersection_Edges.begin()));

					//				influEdgesIntersection[i][j] = intersection_Edges;
					set<int>::iterator iter_set;
					for (iter_set = intersection_Edges.begin();
							iter_set != intersection_Edges.end(); iter_set++) {
						edgesAndInfluEdges[j].insert(*iter_set);
					}
				}
			}
		}

		cout << "finish 2" << endl;
	} else {
		ifstream file(input_influEdges_file);
		string line;
		while (getline(file, line)) {
			vector<string> line_vec;
			split(line, line_vec, ';');
			int edgeId = stoi(line_vec[0]);
			vector<string> line_vec_1;
			set<int> influEdges;
			if (line_vec.size() > 1) {
				split(line_vec[1], line_vec_1, ',');
				int influEdgesNum = line_vec_1.size();
				for (int i = 0; i < influEdgesNum; i++) {
					influEdges.insert(stoi(line_vec_1[i]));
				}
			}
			edgesAndInfluEdges[edgeId] = influEdges;
		}
	}
	time2 = clock();
//	find the top-k coverage based on edgesAndInfluEdges
	set<int> *result = new set<int> [3];

	if (algFlag == 0) {
		result = greedy(K, edgesNum, edgesAndInfluEdges);
	} else if (algFlag == 1) {
		result = dp(K, edgesNum, M, edgesAndInfluEdges, input_partition_file,
				parType);
	} else if (algFlag == 2) {
		result = selection(K, edgesNum, M, edgesAndInfluEdges,
				input_partition_file, parType);
	} else if (algFlag == 3) {
		result = partition(K, edgesNum, M, edgesAndInfluEdges,
				input_partition_file, parType);
	} else if (algFlag == 4) {
		result = greedySelSampling(K, edgesNum, edgesAndInfluEdges, city);
	} else if (algFlag == 5) {
		result = topK_minCut(K, edgesNum, city, edgesAndInfluEdges,
				input_partition_minCut_file, edges_map, edge_density_map);
	}

	time3 = clock();
	duration0 = (double) (time1 - time0) / CLOCKS_PER_SEC;
	duration1 = (double) (time2 - time1) / CLOCKS_PER_SEC;
	double duration2 = (double) (time3 - time2) / CLOCKS_PER_SEC;

	set<int> result_seedEdges = result[0];
	set<int> result_influEdges = result[1];
	cout << "The number of covered edges is: " << result_influEdges.size()
			<< endl;
	cout << duration01 << " " << duration0 << " " << duration1 << " "
			<< duration2 << endl;
	string resultStr = to_string(algFlag) + " " + to_string(pruneFlag) + " "
			+ to_string(result_influEdges.size());
	cout << "The seed edges are as follows: " << result_seedEdges.size()
			<< endl;
	set<int>::iterator iter_seed;
	for (iter_seed = result_seedEdges.begin();
			iter_seed != result_seedEdges.end(); iter_seed++) {
		cout << *iter_seed << " ";
	}
	cout << endl;
	cout << "The influenced edges by the seed edges are as follows: "
			<< result_influEdges.size() << endl;
	set<int>::iterator iter_influ;
	for (iter_influ = result_influEdges.begin();
			iter_influ != result_influEdges.end(); iter_influ++) {
		cout << *iter_influ << " ";
	}

	cout << endl;
	writeOutputFile(output_file, resultStr);
	string timeStr = to_string(duration01) + " " + to_string(duration0) + " "
			+ to_string(duration1) + " " + to_string(duration2);
	writeOutputFile(output_file, timeStr);
	return result;
}

set<int> *greedy(int K, int edgesNum, map<int, set<int>> edgesAndInfluEdges) {
	list<int>::iterator iter_list;
	set<int> result_seedEdges;
	set<int> result_influEdges;
	int result_influScore = 0;
	bool edges_flag[edgesNum];
	bool influEdges_flag[edgesNum];

	for (int i = 0; i < edgesNum; i++) {
		edges_flag[i] = false;
		influEdges_flag[i] = false;
	}
	cout << "The number of K is: " << K << endl;
	while (result_seedEdges.size() < K) {
		int maxEdge;
		double maxInfluScore = -1;
		for (int i = 0; i < edgesNum; i++) {
			if (edges_flag[i] == false) { // check whether the current tempEdge has been included in the previous batches
//				cout << "enter" << endl;
				double tempInfluScore = 0;
				if (edgesAndInfluEdges[i].size() != 0) {
					set<int>::iterator iter_set;
					for (iter_set = edgesAndInfluEdges[i].begin();
							iter_set != edgesAndInfluEdges[i].end();
							iter_set++) {
						if (influEdges_flag[*iter_set] == false) {
							tempInfluScore++;
						}
					}
					if (tempInfluScore > maxInfluScore) {
						maxEdge = i;
						maxInfluScore = tempInfluScore;
					}
				}
			}
		}

		result_seedEdges.insert(maxEdge);
		edges_flag[maxEdge] = true;
		set<int>::iterator iter;
		cout << maxEdge << ";";
		for (iter = edgesAndInfluEdges[maxEdge].begin();
				iter != edgesAndInfluEdges[maxEdge].end(); iter++) {
			if (influEdges_flag[*iter] == false) {
				result_influEdges.insert(*iter);
				cout << *iter << ",";
				influEdges_flag[*iter] = true;
			}
		}
		cout << endl;
		result_influScore += maxInfluScore;
	}
	set<int> *result = new set<int> [3];
	result[0] = result_seedEdges;
	result[1] = result_influEdges;
	return result;
}

set<int> *greedySelSampling(int K, int edgesNum,
		map<int, set<int>> edgesAndInfluEdges, string city) {
	map<int, double> edgeAndInfluNum;

	bool edges_flag[edgesNum];
	bool influEdges_flag[edgesNum];
	map<int, set<int>> r_edge_influEdges;
//	int size = edgesNum;
//	int size = 500;
//	int size = 1000;

	map<int, set<int>> numArr;
	for (int i = 0; i < edgesNum; i++) {
		set<int> value;
		numArr[i] = value;
	}

	for (int i = 0; i < edgesNum; i++) {
		influEdges_flag[i] = false;
		edges_flag[i] = false;
		set<int> influEdges = edgesAndInfluEdges[i];
		set<int>::iterator iter;
		for (iter = influEdges.begin(); iter != influEdges.end(); iter++) {
			r_edge_influEdges[*iter].insert(i);
			numArr[i].insert(*iter);
		}
	}

//	int maxNum = 0;
//
//	for (int i = 0; i < edgesNum; i++) {
//		if (numArr[i].size() > maxNum) {
//			maxNum = numArr[i].size();
//		}
//	}
//
//	double eplison = 0.02;
//	double delta = 0.05;
//	double part1 = 2 * (double) maxNum / (double) edgesNum + eplison;
//	double part2 = eplison * eplison;
//	double part3 = log10(2.0 / delta);
//	int size = part1 / part2 * part3;
//
//	cout << "hha" << " " << edgesNum << " " << part1 << " " << part2 << " "
//			<< part3 << " " << " " << K << " " << size << endl;
	int size;

//	if (city == "chengdu") {
//		size = 2000;
//	} else if (city == "xian") {
//		size = 1500;
//	} else if (city == "porto") {
//		size = 70000;
//	}

	if (city == "chengdu") {
		size = 2000;	// 6135
	} else if (city == "xian") {
		size = 1500;	// 5045
	} else if (city == "porto") {
		size = 30000;	//108571
	}
	set<int>::iterator iter_set;
	set<int> result_seedEdges;
	set<int> result_influEdges;
	double result_influScore = 0;
	while (result_seedEdges.size() < K) {
		int maxEdge;
		int maxInfluNum = 0;
		for (int j = 0; j < size; j++) {
			set<int> temp_edges = r_edge_influEdges[j]; // j: influEdge
			if (influEdges_flag[j] == false) {
				for (iter_set = temp_edges.begin();
						iter_set != temp_edges.end(); iter_set++) {
					int tempEdge = *iter_set;  // tempEdge: seedEdge
					int tempInfluNum = 0;
//					if (tempEdge < edgesNum) {
					if (edges_flag[tempEdge] == false) { // check whether the current tempEdge has been included in the previous batches
						set<int> influEdges = edgesAndInfluEdges[tempEdge];
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
//					}
				}
			}
//			cout << maxEdge << " " << maxInfluNum << endl;
		}

//		can be improved
		result_seedEdges.insert(maxEdge);
		edges_flag[maxEdge] = true;
		cout << maxEdge << ";";
		set<int>::iterator iter;
		for (iter = edgesAndInfluEdges[maxEdge].begin();
				iter != edgesAndInfluEdges[maxEdge].end(); iter++) {
			if (influEdges_flag[*iter] == false) {
				result_influEdges.insert(*iter);
				cout << *iter << ",";
				influEdges_flag[*iter] = true;
			}
		}
		cout << endl;
		result_influScore += maxInfluNum;
//		cout << "result_influScore: " << result_influScore << endl;
	}
	set<int> *result = new set<int> [3];
	result[0] = result_seedEdges;
	result[1] = result_influEdges;
//	cout << result_influScore << " ; actual: " << result_influEdges.size()
//			<< endl;
	return result;
}

set<int> *dp(int K, int edgesNum, int M, map<int, set<int>> edgesAndInfluEdges,
		string input_partition_file, string parType) {
	map<int, double> edgeAndInfluScore;
	bool edges_flag[edgesNum];
	bool influEdges_flag[edgesNum];

	for (int i = 0; i < edgesNum; i++) {
		edges_flag[i] = false;
		influEdges_flag[i] = false;
		set<int> influEdges = edgesAndInfluEdges[i];
		edgeAndInfluScore[i] = influEdges.size();
	}

//	string input_densityPartition = baseFile + city + "/partition/"
//			+ to_string(hour) + ".txt.part." + to_string(M);
//	cout << input_densityPartition << endl;

//	Method 1: follow the Ping's work
//	list<set<int>> mClusters = getThetaPartition(theta, lambda, input_file,
//			edge_influEdges, edge_density_map, new_edge_density_map);
//	M = mClusters.size();

//	Method 2: partition the graph by the traffic density/volume
//	cluster_info clusterInfo = getDensityPartition(M, edgesNum,
//			input_densityPartition);
//	map<int, set<int>> mClusters = clusterInfo.mClusters;
//	map<int, int> edgeAndPartition = clusterInfo.edgeAndPartition;
//	cout << "hhhhh: " << mClusters.size() << " " << edgeAndPartition.size()
//			<< endl;

//	Method 3: partition the graph by the k-means method
	cluster_info clusterInfo;
	clock_t time0, time1, time2, time3;
	if (parType == "kmeans") {
		clusterInfo = getKMeansPartition(M, input_partition_file,
				edgesAndInfluEdges);
	} else if (parType == "grid") {
		clusterInfo = getKMeansPartition(M, input_partition_file,
				edgesAndInfluEdges);
	} else {
		clusterInfo = getPartitionMinCut(M, input_partition_file);
	}
	map<int, set<int>> mClusters = clusterInfo.mClusters; // clusterId(0): edges in the cluster
	map<int, int> edgeAndPartition = clusterInfo.edgeAndPartition;

	double duration1 = 0.0;
	time0 = clock();

	set<int> result_seedEdges;
	set<int> result_influEdges;
	double result_influScore = 0;
//	score_arr[m][k] denotes the influence score of current found k seed edges from the first m partitions
	double score_arr[M + 1][K + 1];
//	bit_arr[m][k] means the partition id where the k-th seed edge resides in which cluster of the first m clusters
	int bit_arr[M + 1][K + 1];  // the cluster id starts from 1
	for (int k = 1; k < K + 1; k++) {
		score_arr[0][k] = 0;
		bit_arr[0][k] = 0;
	}
	for (int m = 1; m < M + 1; m++) {
		score_arr[m][0] = 0;
	}
	for (int k = 1; k < K + 1; k++) {
		list<set<int>>::iterator iter_list;
		int edgePos[M + 1];
		for (int m = 1; m < M + 1; m++) {
			edgePos[m] = 0;
		}
		for (int m = 1; m < M + 1; m++) {
			set<int> edgesInCluster = mClusters[m - 1]; // edges in a cluster
			set<int>::iterator key_iter_edge;
			int maxEdge_inner;
			double maxInfluScore_inner = -1;
			time2 = clock();
			for (key_iter_edge = edgesInCluster.begin();
					key_iter_edge != edgesInCluster.end(); key_iter_edge++) {
				int edgeId = *key_iter_edge;
				double tempInfluScore_inner = 0;
				// judge whether edgeId has been considered as a seed edge
				if (edges_flag[edgeId] == false) {
					set<int> influEdges = edgesAndInfluEdges[edgeId];
					set<int>::iterator iter_influEdges;
					for (iter_influEdges = influEdges.begin();
							iter_influEdges != influEdges.end();
							iter_influEdges++) {
						int influEdgeId = *iter_influEdges;
						if (edgeAndPartition[edgeId]
								== edgeAndPartition[influEdgeId]) {
							if (influEdges_flag[influEdgeId] == false) {
								tempInfluScore_inner++;
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
//			double actual_maxInfluScore_inner = maxInfluScore_inner;
			// judge whether edgeId has been considered as a seed edge
			double actual_maxInfluScore_inner = 0;
			set<int> influEdges = edgesAndInfluEdges[maxEdge_inner];
			set<int>::iterator iter_influEdges;
			for (iter_influEdges = influEdges.begin();
					iter_influEdges != influEdges.end(); iter_influEdges++) {
				int influEdgeId = *iter_influEdges;
				if (influEdges_flag[influEdgeId] == false) {
					actual_maxInfluScore_inner++;
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
				result_seedEdges.insert(edge_seed);
				edges_flag[edge_seed] = true;
				set<int> influEdges = edgesAndInfluEdges[edge_seed];
				set<int>::iterator iter_influEdges;
				for (iter_influEdges = influEdges.begin();
						iter_influEdges != influEdges.end();
						iter_influEdges++) {
					int influEdgeId = *iter_influEdges;
					result_influEdges.insert(influEdgeId);
					influEdges_flag[influEdgeId] = true;
				}
				cout << "We are coming: " << edge_seed << " "
						<< result_influEdges.size() << endl;
			}
		} // end the second loop
	} // end the first loop

	double score = 0.0;
	set<int> result_influEdges1;

	set<int>::iterator iter_result;
	for (iter_result = result_seedEdges.begin();
			iter_result != result_seedEdges.end(); iter_result++) {
		int edgeId = *iter_result;
		cout << edgeId << ";";
		set<int> influEdges = edgesAndInfluEdges[edgeId];
		set<int>::iterator iter_influEdges;
		for (iter_influEdges = influEdges.begin();
				iter_influEdges != influEdges.end(); iter_influEdges++) {
			int influEdgeId = *iter_influEdges;
			if (influEdges_flag[influEdgeId] == true) {
				cout << influEdgeId << ",";
				result_influEdges1.insert(influEdgeId);
				influEdges_flag[influEdgeId] = true;
			}
		}
		cout << endl;
	}
	cout << "actual score result is: " << score << endl;

	result_influScore = score_arr[M][K];
	time1 = clock();
	double duration = (double) (time1 - time0) / CLOCKS_PER_SEC;
	string runtime_d = "DP: " + to_string(duration) + " " + to_string(duration1)
			+ " seconds" + "\n" + to_string(result_influEdges.size()) + " "
			+ to_string(result_influScore) + " "
			+ to_string(result_seedEdges.size());
	cout << runtime_d << endl;

	set<int> *result = new set<int> [3];
	result[0] = result_seedEdges;
	result[1] = result_influEdges;
	return result;
}

set<int> *selection(int K, int edgesNum, int M,
		map<int, set<int>> edgesAndInfluEdges, string input_partition_file,
		string parType) {
	clock_t time0, time1, time2, time3;
	double duration1 = 0.0;
	time0 = clock();
//  Step1: get the partition info by kMeans
//	int hour = 9;
//	int M = 100;
//	string city = "chengdu";
//	string baseFile = "dataset/";
//	string pType = "max";
//	string input_kMeansPartition = baseFile + city + "/partition/kMeans_"
//			+ to_string(M) + "_" + to_string(hour) + "_" + pType + ".txt";

//	Method 1: follow the Ping's work
//	list<set<int>> mClusters = getThetaPartition(theta, lambda, input_file,
//			edge_influEdges, edge_density_map, new_edge_density_map);
//	M = mClusters.size();

//	Method 2: partition the graph by the traffic density/volume
//	cluster_info clusterInfo = getDensityPartition(M, edgesNum,
//			input_densityPartition);
//	map<int, set<int>> mClusters = clusterInfo.mClusters;
//	map<int, int> edgeAndPartition = clusterInfo.edgeAndPartition;
//	cout << "hhhhh: " << mClusters.size() << " " << edgeAndPartition.size()
//			<< endl;

//	Method 3: partition the graph by the k-means method
//	cluster_info clusterInfo = getKMeansPartition(M, input_partition_file,
//			edgesAndInfluEdges);
	cluster_info clusterInfo;
	if (parType == "kmeans") {
		clusterInfo = getKMeansPartition(M, input_partition_file,
				edgesAndInfluEdges);
	} else if (parType == "grid") {
		clusterInfo = getKMeansPartition(M, input_partition_file,
				edgesAndInfluEdges);
	} else {
		clusterInfo = getPartitionMinCut(M, input_partition_file);
	}
	map<int, set<int>> mClusters = clusterInfo.mClusters;
	map<int, int> edgeAndPartition = clusterInfo.edgeAndPartition;
	cout << "finish loading the clusters info: " << mClusters.size() << endl;

//	Step2: start the selection process
	set<int> result_seedEdges;
	set<int> result_influEdges;
	double result_influScore = 0.0;
	map<int, double> clusterAndEstScore;
	map<int, list<int>> sorted_lists;
	map<int, int> edgeAndCluster;

	bool edges_flag[edgesNum];
	bool influEdges_flag[edgesNum];

	list<int> max_list;

	for (int i = 0; i < edgesNum; i++) {
		edges_flag[i] = false;
		influEdges_flag[i] = false;
		edgeAndEstScore_d[i] = edgesAndInfluEdges[i].size();
	}

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
			if (edgeAndEstScore_d[edgeIdInCluster] > estScoreInCluster) {
				estScoreInCluster = edgeAndEstScore_d[edgeIdInCluster];
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
			set<int> influEdges = edgesAndInfluEdges[edgeIdInCluster];
			set<int>::iterator iter_set_influEdge;
			double edgeActScore_temp = 0.0;
			for (iter_set_influEdge = influEdges.begin();
					iter_set_influEdge != influEdges.end();
					iter_set_influEdge++) {
				if (influEdges_flag[*iter_set_influEdge] == false) {
//					edgeActScore_temp +=
//							new_edge_density_map[*iter_set_influEdge]
//									- edge_density_map[*iter_set_influEdge];
					edgeActScore_temp++;
				}
			}
//			edgeAndActScore[edgeIdInCluster] = edgeActScore_temp;
			edgeAndEstScore_d[edgeIdInCluster] = edgeActScore_temp;
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
					if (edgeAndEstScore_d[*iter_list_s] > edgeActScore_max) {
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
			if (edgeActScore_max < edgeAndEstScore_d[est_max_edgeId]) {
				set<int> influEdges = edgesAndInfluEdges[est_max_edgeId];
				set<int>::iterator iter_set_influEdge;
				double edgeActScore_temp = 0.0;
				for (iter_set_influEdge = influEdges.begin();
						iter_set_influEdge != influEdges.end();
						iter_set_influEdge++) {
					if (influEdges_flag[*iter_set_influEdge] == false) {
						edgeActScore_temp++;
//						edgeActScore_temp +=
//								new_edge_density_map[*iter_set_influEdge]
//										- edge_density_map[*iter_set_influEdge];
					}
				}
//				update the upper bound of other edges with large scores (lines 16-18)
				edgeAndEstScore_d[est_max_edgeId] = edgeActScore_temp;
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
			clusterAndEstScore[m] = edgeAndEstScore_d[*sorted_lists[m].begin()];
		}

//		push into seedEdges set
		result_seedEdges.insert(edgeId_max);
		cout << edgeId_max << ";";
		edges_flag[edgeId_max] = true;
//		cout << result_seedEdges.size() << ": " << clusterId_max << ": "
//				<< edgeId_max << ": " << result_influEdges.size() << ": "
//				<< edge_influEdges[edgeId_max].size() << endl;
		mClusters[clusterId_max].erase(edgeId_max);
		sorted_lists[clusterId_max].remove(edgeId_max);
		result_influScore += edgeActScore_max;
		set<int> influEdges = edgesAndInfluEdges[edgeId_max];
		set<int>::iterator iter_set_influEdge;
		for (iter_set_influEdge = influEdges.begin();
				iter_set_influEdge != influEdges.end(); iter_set_influEdge++) {
			if (influEdges_flag[*iter_set_influEdge] == false) {
				result_influEdges.insert(*iter_set_influEdge);
				influEdges_flag[*iter_set_influEdge] = true;
				cout << *iter_set_influEdge << ",";
			}
		}
		cout << endl;

	} // end while loop to get K seed edges

	double score = 0.0;
	set<int> result_influEdges1;
	vector<int>::iterator iter_result;

	time1 = clock();
	double duration = (double) (time1 - time0) / CLOCKS_PER_SEC;
	string runtime_s = "Selection: " + to_string(duration) + " "
			+ to_string(duration1) + " seconds" + "\n"
			+ to_string(result_influEdges.size()) + " "
			+ to_string(result_influEdges1.size()) + " "
			+ to_string(result_influScore) + " " + to_string(score);

	cout << runtime_s << endl;
	set<int> *result = new set<int> [3];
	result[0] = result_seedEdges;
	result[1] = result_influEdges;
	return result;
}

set<int> *partition(int K, int edgesNum, int M,
		map<int, set<int>> edgesAndInfluEdges, string input_partition_file,
		string parType) {
//	Method 3: partition the graph by the k-means method
	cluster_info clusterInfo;
	if (parType == "kmeans") {
		clusterInfo = getKMeansPartition(M, input_partition_file,
				edgesAndInfluEdges);
	} else if (parType == "grid") {
		clusterInfo = getKMeansPartition(M, input_partition_file,
				edgesAndInfluEdges);
	} else {
		clusterInfo = getPartitionGrid(M, input_partition_file);
	}
	map<int, set<int>> mClusters = clusterInfo.mClusters; // clusterId(0): edges in the cluster
//	map<int, int> edgeAndPartition = clusterInfo.edgeAndPartition;

	clock_t time0, time1, time2, time3, time4;
	time0 = clock();
	set<int> result_seedEdges;
	set<int> result_influEdges;
	map<int, int> edgeAndInfluNum;

//	S \xi
	map<int, map<int, set<int>>> seedEdgesInOneBatch;
//	\Phi(S) \mathbb{I}
	map<int, map<int, set<int>>> seedEdgesInSomeBatch;

	map<int, map<int, set<int>>> influEdgesInOneBatch;
	map<int, map<int, set<int>>> influEdgesInSomeBatch;

	for (int i = 0; i <= M; i++) {
		for (int k = 0; k <= K; k++) {
			set<int> edges;
			seedEdgesInOneBatch[i][k] = edges;
			seedEdgesInSomeBatch[i][k] = edges;
			influEdgesInOneBatch[i][k] = edges;
			influEdgesInSomeBatch[i][k] = edges;
		}
	}

	for (int i = 1; i <= M; i++) {
		for (int k = 1; k <= K; k++) {
			set<int> s_hat = seedEdgesInSomeBatch[i - 1][K - k];
			set<int> i_hat = influEdgesInSomeBatch[i - 1][K - k];

			set<int> *resultByBatchSel = greedySelForPartition(k, mClusters[i],
					edgesNum, edgesAndInfluEdges, s_hat, i_hat);

			seedEdgesInOneBatch[i][k] = resultByBatchSel[0];
			influEdgesInOneBatch[i][k] = resultByBatchSel[1];
			int max_influNum = -1;
			int max_j = -1;
			for (int j = 0; j <= k; j++) {
//				int temp_influNum = influEdgesInSomeBatch[i - 1][k - j].size()
//						+ influEdgesInOneBatch[i][j].size();
				set<int> union_influEdges;
				set<int> iEdges1 = influEdgesInSomeBatch[i - 1][k - j];
				set<int> iEdges2 = influEdgesInOneBatch[i][j];
				set_union(iEdges1.begin(), iEdges1.end(), iEdges2.begin(),
						iEdges2.end(),
						inserter(union_influEdges, union_influEdges.begin()));
				int temp_influNum = union_influEdges.size();

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

			cout << "max_j: " << max_j << ";i: " << i << ";k: " << k
					<< "; seedNum: " << seedEdgesInOneBatch[i][k].size()
					<< "; influNum: " << influEdgesInOneBatch[i][k].size()
					<< "; seedNum: " << seedEdgesInSomeBatch[i][k].size()
					<< "; influNum: " << influEdgesInSomeBatch[i][k].size()
					<< endl;
		}
	}
	result_seedEdges = seedEdgesInSomeBatch[M][K];
	result_influEdges = influEdgesInSomeBatch[M][K];

	time1 = clock();
	set<int> result_influEdges1;
	set<int>::iterator iter_result;
	for (iter_result = result_seedEdges.begin();
			iter_result != result_seedEdges.end(); iter_result++) {
		int edgeId = *iter_result;
		cout << edgeId << ";";
		set<int> influEdges = edgesAndInfluEdges[edgeId];
		set<int>::iterator iter_influEdges;
		for (iter_influEdges = influEdges.begin();
				iter_influEdges != influEdges.end(); iter_influEdges++) {
			int influEdgeId = *iter_influEdges;
			if (result_influEdges1.find(influEdgeId)
					== result_influEdges1.end()) {
				cout << influEdgeId << ",";
				result_influEdges1.insert(influEdgeId);
			}
		}
		cout << endl;
	}
	cout << "actual score result is: " << result_influEdges1.size() << endl;

	double duration = (double) (time1 - time0) / CLOCKS_PER_SEC;
	string runtime_p = "Partition: " + to_string(duration) + " seconds" + "\n"
			+ to_string(result_influEdges.size()) + " "
			+ to_string(result_influEdges1.size()) + " "
			+ to_string(result_seedEdges.size());
	cout << runtime_p << endl;

	set<int> *result = new set<int> [3];
	result[0] = result_seedEdges;
	result[1] = result_influEdges;
	return result;
}

set<int> *greedySelForPartition(int k, set<int> edges, int edgesNum,
		map<int, set<int>> edge_influEdges, set<int> s_hat, set<int> i_hat) {
	set<int> result_seedEdges;
	set<int> result_influEdges;
	set<int> *result = new set<int> [3];

	bool edges_flag[edgesNum];
	bool influEdges_flag[edgesNum];

	for (int i = 0; i < edgesNum; i++) {
		edges_flag[i] = false;
		influEdges_flag[i] = false;
	}
	set<int>::iterator iter_s_hat;
	set<int>::iterator iter_i_hat;
	int edgesInSeeds = 0;
	int edgesInflu = 0;
	for (iter_s_hat = s_hat.begin(); iter_s_hat != s_hat.end(); iter_s_hat++) {
		if (edges.find(*iter_s_hat) != edges.end()) {
			edgesInSeeds++;
		}
		edges_flag[*iter_s_hat] = true;
	}
	for (iter_i_hat = i_hat.begin(); iter_i_hat != i_hat.end(); iter_i_hat++) {
		if (edges.find(*iter_i_hat) != edges.end()) {
			edgesInflu++;
		}
		influEdges_flag[*iter_i_hat] = true;
	}
	if (edges.size() - edgesInSeeds < k) {
		result[0] = result_seedEdges;
		result[1] = result_influEdges;
	} else {
		set<int>::iterator iter_list;

		set<int> result_influEdgesInSomeBatches = i_hat;
		int result_influScore = 0;
		while (result_seedEdges.size() < k) {
			int maxEdge;
			double maxInfluScore = -1;
			for (iter_list = edges.begin(); iter_list != edges.end();
					iter_list++) {
				int tempEdge = *iter_list;
				if (edges_flag[tempEdge] == false) { // check whether the current tempEdge has been included in the previous batches
					double tempInfluScore = 0;
					set<int>::iterator iter_set;
					for (iter_set = edge_influEdges[tempEdge].begin();
							iter_set != edge_influEdges[tempEdge].end();
							iter_set++) {
						if (influEdges_flag[*iter_set] == false) {
							tempInfluScore++;
						}
					}
					if (tempInfluScore > maxInfluScore) {
						maxEdge = tempEdge;
						maxInfluScore = tempInfluScore;
					}
				}
			}
//			cout << "****" << maxEdge << " " << maxInfluScore << endl;
//			can be improved
			result_seedEdges.insert(maxEdge);
			edges_flag[maxEdge] = true;
			edges.erase(maxEdge);
			set<int>::iterator iter;
			for (iter = edge_influEdges[maxEdge].begin();
					iter != edge_influEdges[maxEdge].end(); iter++) {
				if (influEdges_flag[*iter] == false) {
					result_influEdges.insert(*iter);
//					result_influEdgesInSomeBatches.insert(*iter);
					influEdges_flag[*iter] = true;
				}
			}
			result_influScore += maxInfluScore;
		}
		result[0] = result_seedEdges;
		result[1] = result_influEdges; // the influenced edges in the i-th batches
//		result[2] = result_influEdgesInSomeBatches; // the influenced edges in the first i batches
		cout << result_influScore << " ; actual: " << result_influEdges.size()
				<< endl;
	}
	return result;
}

set<int> *topK_minCut(int K, int edgesNum, string city,
		map<int, set<int>> edgesAndInfluEdges, string input_partition_file,
		map<int, Edge> edges_map, map<int, double> edge_density_map) {
	set<int> result_seedEdges;
	set<int> result_influEdges;
	set<int> *result = new set<int> [3];
	clock_t time0, time1, time2, time3, time4;
	cout << "The number of edges is: " << edgesNum << endl;

//	get the edge cuts after partitions
	set<int> edgesCuts = getMinCutsEdges(edgesNum, edges_map,
			input_partition_file);
	cout << "The number of cuts is: " << edgesCuts.size() << endl;

	vector<pair<int, double>> vtMap;
	set<int>::iterator it;
	for (it = edgesCuts.begin(); it != edgesCuts.end(); it++) {
		int edgeId = *it;
		vtMap.push_back(make_pair(edgeId, edge_density_map[edgeId]));
	}

	sort(vtMap.begin(), vtMap.end(),
			[](const pair<int, double> &x, const pair<int, double> &y) -> int {
				return x.second > y.second;
			});

	int result_influScore = 0;
	bool influEdges_flag[edgesNum];

	for (int i = 0; i < edgesNum; i++) {
		influEdges_flag[i] = false;
	}
	set<int>::iterator iter_set;
	int seedEdges = 0;
	double resultInfluScore = 0;
	for (auto it = vtMap.begin(); it != vtMap.end(); it++) {
		int i = it->first; // i is the edgeId.
		cout << i << ";";
		result_seedEdges.insert(i);
		double tempInfluScore = 0;
		if (edgesAndInfluEdges[i].size() != 0) {
			set<int>::iterator iter_set;
			for (iter_set = edgesAndInfluEdges[i].begin();
					iter_set != edgesAndInfluEdges[i].end(); iter_set++) {
				if (influEdges_flag[*iter_set] == false) {
					tempInfluScore++;
					result_influEdges.insert(*iter_set);
					cout << *iter_set << ",";
					influEdges_flag[*iter_set] = true;
				}
			}
		}
		cout << endl;
		resultInfluScore += tempInfluScore;

		seedEdges++;
		if (seedEdges == K) {
			break;
		}
	}
	result[0] = result_seedEdges;
	result[1] = result_influEdges;
	cout << "The number of influenced edges is: " << resultInfluScore << endl;
	return result;
}

set<int> getMinCutsEdges(int edgesNum, map<int, Edge> edges,
		string input_node_partition) {
	set<int> notInPartitionEdges;
	map<int, int> nodeAndPartition = getNodePartition(input_node_partition);
	map<int, set<int>> partitionAndEdges;
	for (int i = 0; i < edgesNum; i++) {
		Edge edge = edges[i];
		int edgeId = edge.edgeId;
		int souceNode = edge.sourceNode;
		int desNode = edge.desNode;
		int sourcePid = nodeAndPartition[souceNode];
		int desPid = nodeAndPartition[desNode];
		if (sourcePid != desPid) {
			notInPartitionEdges.insert(edgeId);
		}
	}
	return notInPartitionEdges;
}
