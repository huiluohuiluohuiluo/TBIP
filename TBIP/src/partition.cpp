#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iterator>
#include <map>
#include <set>
#include <bits/stdc++.h>

#include "time.h"
#include "edge.h"
#include "fileProcess.h"
//#include "cluster.h"

using namespace std;

/**
 * Given a set, return all the subclusters
 */
list<set<int>> getSubClusters(set<int> edges) {
	list<set<int>> subclusters;
	set<int>::iterator iter_set;
	int length = edges.size();
	int masks[length];
	for (int i = 0; i < length; i++) {
		masks[i] = (1 << i);
	}
	for (int i = 0; i < (1 << length); i++) {
		set<int> subcluster;
		int j = 0;
		for (iter_set = edges.begin(); iter_set != edges.end(); iter_set++) {
			if ((masks[j] & i) != 0) {
				subcluster.insert(*iter_set);
			}
			j++;
		}
		if (subcluster.size() > 0) {
			subclusters.push_back(subcluster);
		}
	}
	return subclusters;
}

/**
 * given a set of edges, compute their influence scores
 */
double getEdgesInfluScore(double lambda, set<int> edges,
		map<int, set<int>> edge_influEdges, map<int, double> edge_density_map,
		map<int, double> new_edge_density_map) {

	set<int>::iterator iter_set1;
	set<int> influEdges;
	for (iter_set1 = edges.begin(); iter_set1 != edges.end(); iter_set1++) {
		set<int>::iterator iter_s;
		for (iter_s = edge_influEdges[*iter_set1].begin();
				iter_s != edge_influEdges[*iter_set1].end(); iter_s++) {
			influEdges.insert(*iter_s);
		}
	}

	set<int>::iterator iter_set2;
	double score = 0.0;
	for (iter_set2 = influEdges.begin(); iter_set2 != influEdges.end();
			iter_set2++) {
		if (new_edge_density_map[*iter_set2] > lambda
				&& edge_density_map[*iter_set2] >= 0
				&& new_edge_density_map[*iter_set2]
						> edge_density_map[*iter_set2]) {
			score += new_edge_density_map[*iter_set2]
					- edge_density_map[*iter_set2];
		}
	}
	return score;
}

/**
 * given an edge, compute their influence scores
 */
double getEdgeInfluScore(double lambda, set<int> influEdges,
		map<int, double> edge_density_map,
		map<int, double> new_edge_density_map) {
	double score = 0.0;
	set<int>::iterator iter_set;

	for (iter_set = influEdges.begin(); iter_set != influEdges.end();
			iter_set++) {
		if (new_edge_density_map[*iter_set] > lambda
				&& edge_density_map[*iter_set] >= 0
				&& new_edge_density_map[*iter_set]
						> edge_density_map[*iter_set]) {
			score += new_edge_density_map[*iter_set]
					- edge_density_map[*iter_set];
		}
	}
	return score;
}

set<int> getEdgesInfluEdges(set<int> edges,
		map<int, set<int>> edge_influEdges) {
	set<int> influEdges;
	set<int>::iterator iter_set1;
	for (iter_set1 = edges.begin(); iter_set1 != edges.end(); iter_set1++) {
		set<int> e_influEdges = edge_influEdges[*iter_set1];
		set<int>::iterator iter_set2;
		for (iter_set2 = e_influEdges.begin(); iter_set2 != e_influEdges.end();
				iter_set2++) {
			influEdges.insert(*iter_set2);
		}
	}
	return influEdges;
}

double getJaccardDistance(set<int> cEdges1, set<int> cEdges2) {
	double ratio = 0.0;
	set<int> intersectionEdges;
	set<int> unionEdges;

	set_intersection(cEdges1.begin(), cEdges1.end(), cEdges2.begin(),
			cEdges2.end(),
			inserter(intersectionEdges, intersectionEdges.begin()));

	set_union(cEdges1.begin(), cEdges1.end(), cEdges2.begin(), cEdges2.end(),
			inserter(unionEdges, unionEdges.begin()));

	if (unionEdges.size() > 0) {
		ratio = (double) intersectionEdges.size() / (double) unionEdges.size();
	}
	return ratio;
}

double getOverlapRatio(double lambda, double theta, set<int> cEdges1,
		set<int> cEdges2, map<int, set<int>> edge_influEdges,
		map<int, double> edge_density_map,
		map<int, double> new_edge_density_map) {
	double ratio = 0.0;
	bool overlap_flag = false;

	set<int> influEdges1 = getEdgesInfluEdges(cEdges1, edge_influEdges);
	set<int> influEdges2 = getEdgesInfluEdges(cEdges2, edge_influEdges);
	set<int> intersectionEdges;
	set_intersection(influEdges1.begin(), influEdges1.end(),
			influEdges2.begin(), influEdges2.end(),
			inserter(intersectionEdges, intersectionEdges.begin()));
	if (intersectionEdges.size() > 0) {
		overlap_flag = true;
	}
	if (overlap_flag == true) {
		double influScore2 = getEdgesInfluScore(lambda, cEdges2,
				edge_influEdges, edge_density_map, new_edge_density_map);
		// given a list cluster, return a list of subclusters

		list<set<int>> subClusters_list = getSubClusters(cEdges1);
		list<set<int>>::iterator iter_list;
		int i = 0;
		for (iter_list = subClusters_list.begin();
				iter_list != subClusters_list.end(); iter_list++) {
			set<int> subsetEdges = *iter_list;
			set<int> unionEdges;
			set_union(subsetEdges.begin(), subsetEdges.end(), cEdges2.begin(),
					cEdges2.end(), inserter(unionEdges, unionEdges.begin()));

			double influScore1 = getEdgesInfluScore(lambda, subsetEdges,
					edge_influEdges, edge_density_map, new_edge_density_map);
			double influScore12 = getEdgesInfluScore(lambda, unionEdges,
					edge_influEdges, edge_density_map, new_edge_density_map);

			double temp_ratio = 0.0;
			if (influScore1 != 0) {
				temp_ratio = (influScore1 + influScore2 - influScore12)
						/ influScore1;
				i++;
			}
			if (temp_ratio > ratio) {
				ratio = temp_ratio;
			}
		}
	}
	return ratio;
}

list<set<int>> getDensityPartition(int M, vector<Edge> edges,
		string input_node_partition) {
	list<set<int>> mClusters;
	map<int, int> nodeAndPartition = getNodePartition(input_node_partition);
	map<int, set<int>> partitionAndEdges;
	vector<Edge>::iterator iter_vec;
	for (iter_vec = edges.begin(); iter_vec != edges.end(); iter_vec++) {
		Edge edge = *iter_vec;
		int edgeId = edge.edgeId;
		int souceNode = edge.sourceNode;
		int partitionId = nodeAndPartition[souceNode];
		set<int> value;
		if (partitionAndEdges.find(partitionId) != partitionAndEdges.end()) {
			value = partitionAndEdges[partitionId];
		}
		value.insert(edgeId);
		partitionAndEdges[partitionId] = value;
	}
	for (int i = 0; i < M; i++) {
		mClusters.push_back(partitionAndEdges[i]);
	}
	return mClusters;
}

list<set<int>> getKMeansPartition_list(string input_file,
		map<int, set<int>> old_edge_influEdges) {
	int edgesNum = old_edge_influEdges.size();
	ifstream file(input_file);
	list<set<int>> mClusters;
	map<int, set<int>> edge_influEdges;
	set<int> edges;
	string line;
	while (getline(file, line)) {
		vector<string> line_vec;
		split(line, line_vec, ';');
		set<int> cEdges;
		int cEdgesNum = line_vec.size();
		for (int i = 0; i < cEdgesNum; i++) {
			int edgeId = stoi(line_vec[i]);
			cEdges.insert(edgeId);
		}
		mClusters.push_back(cEdges);
	}
	return mClusters;
}

map<int, int> getEdgePartition_list(string input_file,
		map<int, set<int>> old_edge_influEdges) {
	int edgesNum = old_edge_influEdges.size();
	ifstream file(input_file);
	map<int, int> edgeAndPartition;
	map<int, set<int>> edge_influEdges;
	set<int> edges;
	string line;
	int partitionId = 0;
	while (getline(file, line)) {
		vector<string> line_vec;
		split(line, line_vec, ';');
		int cEdgesNum = line_vec.size();
		for (int i = 0; i < cEdgesNum; i++) {
			int edgeId = stoi(line_vec[i]);
			edgeAndPartition[edgeId] = partitionId;
		}
		partitionId++;
	}
	return edgeAndPartition;
}

map<int, int> getEdgePartition(int M, vector<Edge> edges,
		string input_node_partition) {
	map<int, int> nodeAndPartition = getNodePartition(input_node_partition);
	map<int, int> edgeAndPartition;
	vector<Edge>::iterator iter_vec;
	for (iter_vec = edges.begin(); iter_vec != edges.end(); iter_vec++) {
		Edge edge = *iter_vec;
		int edgeId = edge.edgeId;
		int souceNode = edge.sourceNode;
		int partitionId = nodeAndPartition[souceNode];
		edgeAndPartition[edgeId] = partitionId;
	}
	return edgeAndPartition;
}

list<int> getPartitionEdges(int M, vector<Edge> edges,
		string input_node_partition) {
	list<int> notInPartitionEdges;
	map<int, int> nodeAndPartition = getNodePartition(input_node_partition);
	map<int, set<int>> partitionAndEdges;
	vector<Edge>::iterator iter_vec;
	for (iter_vec = edges.begin(); iter_vec != edges.end(); iter_vec++) {
		Edge edge = *iter_vec;
		int edgeId = edge.edgeId;
		int souceNode = edge.sourceNode;
		int desNode = edge.desNode;
		int sourcePid = nodeAndPartition[souceNode];
		int desPid = nodeAndPartition[desNode];
		if (sourcePid != desPid) {
			notInPartitionEdges.push_back(edgeId);
		}
	}
	return notInPartitionEdges;
}

list<set<int>> getThetaPartition(double theta, double lambda, string input_file,
		map<int, set<int>> old_edge_influEdges,
		map<int, double> edge_density_map,
		map<int, double> new_edge_density_map) {
	int edgesNum = old_edge_influEdges.size();
	ifstream file(input_file);
	map<int, set<int>> edge_influEdges;
	list<set<int>> mClusters;
	if (!file) {
		string t_file = "dataset/chengdu/partition/1.txt";
		ifstream temp_file(t_file);
		string line;
		set<int> traverseEdges;
		while (getline(temp_file, line)) {
			vector<string> line_vec;
			split(line, line_vec, ';');
			set<int> cEdges;
			int cEdgesNum = line_vec.size();
			for (int i = 0; i < cEdgesNum; i++) {
				int edgeId = stoi(line_vec[i]);
				cEdges.insert(edgeId);
				traverseEdges.insert(edgeId);
			}
//			mClusters.push_back(cEdges);
		}

		ofstream output(input_file, ios::app);
		cout << "Create partition file ..." << endl;
		for (int i = 0; i < edgesNum; i++) {
//		for (int i = 0; i < 10; i++) {
			if (old_edge_influEdges[i].size() > 0) { // filter some edges that can not influence other edges
				set<int> cEdges;
				cEdges.insert(i);
				set<int> old_influEdges = old_edge_influEdges[i];
				set<int> influEdges;
				set<int>::iterator iter_influEdges;
				for (iter_influEdges = old_influEdges.begin();
						iter_influEdges != old_influEdges.end();
						iter_influEdges++) {
					if (new_edge_density_map[*iter_influEdges] > lambda
							&& edge_density_map[*iter_influEdges] >= 0
							&& new_edge_density_map[*iter_influEdges]
									> edge_density_map[*iter_influEdges]) {
						influEdges.insert(*iter_influEdges);
					}
				}
				if (influEdges.size() > 0) {
					edge_influEdges[i] = influEdges;
					set<int>::iterator iter_trav = find(traverseEdges.begin(),
							traverseEdges.end(), i); // iterator is used to save the idx

					if (iter_trav == traverseEdges.end()) {
						mClusters.push_back(cEdges);
					}
				}
			}
		}
		cout << "mClusters size is: " << mClusters.size() << endl;
		cout << "edge_influEdges size is: " << edge_influEdges.size() << endl;
		list<set<int>>::iterator iter_list;
		bool flag = false;
		while (true) {
			list<set<int>>::iterator iter_list1;
			list<set<int>>::iterator iter_list2;
			int i = 0;
			for (iter_list1 = mClusters.begin(); iter_list1 != mClusters.end();
					iter_list1++) {
				set<int> cEdges1 = *iter_list1;
				i++;
				int j = 0;
				for (iter_list2 = mClusters.begin();
						iter_list2 != mClusters.end(); iter_list2++) {
					set<int> cEdges2 = *iter_list2;
					j++;
					if (cEdges1 != cEdges2) {
						set<int> unionEdges;
						set_union(cEdges1.begin(), cEdges1.end(),
								cEdges2.begin(), cEdges2.end(),
								inserter(unionEdges, unionEdges.begin()));
						double ratio = getOverlapRatio(lambda, theta, cEdges1,
								cEdges2, edge_influEdges, edge_density_map,
								new_edge_density_map);
						// if the overlap between two clusters is larger than theta, then these two clusters should be merged
//						cout << "ratio is: " << ratio << " " << i << " " << j
//								<< endl;
						if (ratio >= theta) {
//						if (ratio >= 0) {
							// merge two clusters: cluster 1 and cluster 2
							set<int> edges;
							set<int>::iterator iter_set1;
							set<int>::iterator iter_set2;
							for (iter_set1 = cEdges1.begin();
									iter_set1 != cEdges1.end(); iter_set1++) {
								edges.insert(*iter_set1);
							}

							for (iter_set2 = cEdges2.begin();
									iter_set2 != cEdges2.end(); iter_set2++) {
								edges.insert(*iter_set2);
							}
//							cout << "Before ---- mClusters.size() is: "
//									<< mClusters.size() << endl;
							mClusters.remove(cEdges1);
							mClusters.remove(cEdges2);
							mClusters.push_back(edges);
							cout << "After ---- mClusters.size() is: "
									<< mClusters.size() << endl;
							flag = true;
							iter_list1 = mClusters.end();
							iter_list2 = mClusters.end();
							break;
						}
					}
				} // end second loop
				if (flag == false) {
//					cout << "Temp Before ---- mClusters.size() is: "
//							<< mClusters.size() << endl;
					mClusters.remove(cEdges1);
//					temp_mClusters.push_back(cEdges1);
					set<int>::iterator iter_c_edges;
					string line;
					for (iter_c_edges = cEdges1.begin();
							iter_c_edges != cEdges1.end(); iter_c_edges++) {
						line += to_string(*iter_c_edges) + ";";
					}
					line.pop_back();
//					output << line << endl;
//					cout << "Temp After ---- mClusters.size() is: "
//							<< mClusters.size() << endl;
					flag = true;
				}
				break;
			} // end first loop
			if (flag == true) {
				flag = false;
			}
			if (mClusters.size() == 1) {
				break;
			}
		}

//		write into files
//		cout << "Before mClusters.size(): " << mClusters.size() << endl;
//		ofstream output(input_file, ios::app);
//		list<set<int>>::iterator iter_list_t;
//		for (iter_list_t = temp_mClusters.begin();
//				iter_list_t != temp_mClusters.end(); iter_list_t++) {
//			mClusters.push_back(*iter_list_t);
//		}
//		cout << "temp_mClusters.size(): " << temp_mClusters.size() << endl;
//		cout << "After mClusters.size(): " << mClusters.size() << endl;
//		list<set<int>>::iterator iter_list_c;
//		for (iter_list_c = mClusters.begin(); iter_list_c != mClusters.end();
//				iter_list_c++) {
//			string line;
//			set<int> c_edges = *iter_list_c;
//			set<int>::iterator iter_c_edges;
//			for (iter_c_edges = c_edges.begin(); iter_c_edges != c_edges.end();
//					iter_c_edges++) {
//				line += to_string(*iter_c_edges) + ";";
//			}
//			line.pop_back();
//			output << line << endl;
//		}
	} else {
		string line;
		while (getline(file, line)) {
			vector<string> line_vec;
			split(line, line_vec, ';');
			set<int> cEdges;
			int cEdgesNum = line_vec.size();
			for (int i = 0; i < cEdgesNum; i++) {
				int edgeId = stoi(line_vec[i]);
				cEdges.insert(edgeId);
			}
			mClusters.push_back(cEdges);
		}
	}
	return mClusters;
}

map<int, set<int>> getKMeansPartition(int M, double lambda, string input_file,
		map<int, set<int>> old_edge_influEdges,
		map<int, double> edge_density_map,
		map<int, double> new_edge_density_map) {
	int edgesNum = old_edge_influEdges.size();
	ifstream file(input_file);
	map<int, set<int>> mClusters;
	map<int, set<int>> edge_influEdges;
	set<int> edges;
	if (!file) {
		for (int i = 0; i < edgesNum; i++) {
			if (old_edge_influEdges[i].size() > 0) { // filter some edges that can not influence other edges
				set<int> old_influEdges = old_edge_influEdges[i];
				set<int> influEdges;
				set<int>::iterator iter_influEdges;
				for (iter_influEdges = old_influEdges.begin();
						iter_influEdges != old_influEdges.end();
						iter_influEdges++) {
					if (new_edge_density_map[*iter_influEdges] > lambda
							&& edge_density_map[*iter_influEdges] >= 0
							&& new_edge_density_map[*iter_influEdges]
									> edge_density_map[*iter_influEdges]) {
						influEdges.insert(*iter_influEdges);
					}
				}
				if (influEdges.size() > 0) {
					edges.insert(i);
					edge_influEdges[i] = influEdges;
				}
			}
		}
		cout << "actual edge_influEdges size is: " << edge_influEdges.size()
				<< endl;

		bool flag = true;
		set<int> edgesCenter;
		set<int> temp_edgesCenter;
		int iterator = 0;
		while (flag) {
			if (mClusters.size() == 0) {
				list<set<int>>::iterator iter_list;
				int i = 0;
				set<int>::iterator iter_set_edges;
				for (iter_set_edges = edges.begin();
						iter_set_edges != edges.end(); iter_set_edges++) {
					int edgeId = *iter_set_edges;
					if (i < M) {
						set<int> edgesInCluster;
						edgesInCluster.insert(edgeId);
						mClusters[edgeId] = edgesInCluster;
						edgesCenter.insert(edgeId);
						temp_edgesCenter.insert(edgeId);
					} else {
						set<int>::iterator iter_set_edgesCenter;
						int edgeIdCenter_select;
//						double max_jaccardDis = -1;
						double min_jaccardDis = 100.0;
						for (iter_set_edgesCenter = edgesCenter.begin();
								iter_set_edgesCenter != edgesCenter.end();
								iter_set_edgesCenter++) {
							int edgeIdCenter = *iter_set_edgesCenter;
							double temp_jaccardDis = getJaccardDistance(
									edge_influEdges[edgeId],
									edge_influEdges[edgeIdCenter]);
//							if (temp_jaccardDis > max_jaccardDis) {
//								max_jaccardDis = temp_jaccardDis;
//								edgeIdCenter_select = edgeIdCenter;
//							}
							if (temp_jaccardDis < min_jaccardDis) {
								min_jaccardDis = temp_jaccardDis;
								edgeIdCenter_select = edgeIdCenter;
							}
						}
						mClusters[edgeIdCenter_select].insert(edgeId);
					}
					i++;
				}
			} else {
//				update all the center edgeIds
				set<int>::iterator iter_set_edgesCenter;
				bool flag_replace = false;
				for (iter_set_edgesCenter = edgesCenter.begin();
						iter_set_edgesCenter != edgesCenter.end();
						iter_set_edgesCenter++) {
					int edgeIdCenter = *iter_set_edgesCenter;
					set<int> edgesInCluster = mClusters[edgeIdCenter];
					set<int>::iterator iter_set_edgesInCluster;
					int new_edgeIdCenter;
//					double max_influScore = 0.0;
					double min_influScore = 1000000.0;
					for (iter_set_edgesInCluster = edgesInCluster.begin();
							iter_set_edgesInCluster != edgesInCluster.end();
							iter_set_edgesInCluster++) {
						int edgeIdInCluster = *iter_set_edgesInCluster;
						double temp_influScore = getEdgeInfluScore(lambda,
								edge_influEdges[edgeIdInCluster],
								edge_density_map, new_edge_density_map);
//						if (temp_influScore > max_influScore) {
//							max_influScore = temp_influScore;
//							new_edgeIdCenter = edgeIdInCluster;
//						}
						if (temp_influScore < min_influScore) {
							min_influScore = temp_influScore;
							new_edgeIdCenter = edgeIdInCluster;
						}
					}
//					update the center edgeId
					if (new_edgeIdCenter != edgeIdCenter) {
						temp_edgesCenter.erase(edgeIdCenter);
						temp_edgesCenter.insert(new_edgeIdCenter);
						flag_replace = true;
					}
				}
//				update all the clusters based on the new found center edges
				edgesCenter.clear();
				edgesCenter = temp_edgesCenter;
				if (flag_replace == true) {
					mClusters.clear();
					set<int>::iterator iter_set_edgeIdCenter;
					for (iter_set_edgeIdCenter = edgesCenter.begin();
							iter_set_edgeIdCenter != edgesCenter.end();
							iter_set_edgeIdCenter++) {
						int edgeIdCenter = *iter_set_edgeIdCenter;
						set<int> edgesInCluster;
						edgesInCluster.insert(edgeIdCenter);
						mClusters[edgeIdCenter] = edgesInCluster;
					}
					set<int>::iterator iter_set_edge;
					for (iter_set_edge = edges.begin();
							iter_set_edge != edges.end(); iter_set_edge++) {
						int edgeIdNotCenter = *iter_set_edge;
						set<int>::iterator iter_set_find = find(
								edgesCenter.begin(), edgesCenter.end(),
								edgeIdNotCenter);
						if (iter_set_find == edgesCenter.end()) {
//							double max_jaccardDis = -1;
							double min_jaccardDis = 100.0;
							int edgeIdCenter_select;
							set<int>::iterator iter_set_edgesCenter2;
							for (iter_set_edgesCenter2 = edgesCenter.begin();
									iter_set_edgesCenter2 != edgesCenter.end();
									iter_set_edgesCenter2++) {
								int edgeIdCenter = *iter_set_edgesCenter2;
								double temp_jaccardDis = getJaccardDistance(
										edge_influEdges[edgeIdNotCenter],
										edge_influEdges[edgeIdCenter]);
//								if (temp_jaccardDis > max_jaccardDis) {
//									max_jaccardDis = temp_jaccardDis;
//									edgeIdCenter_select = edgeIdCenter;
//								}
								if (temp_jaccardDis < min_jaccardDis) {
									min_jaccardDis = temp_jaccardDis;
									edgeIdCenter_select = edgeIdCenter;
								}
							}
							mClusters[edgeIdCenter_select].insert(
									edgeIdNotCenter);
						}
					}
				} else {
					cout << "iterator: " << iterator << " ; mClusters.size(): "
							<< mClusters.size() << "; edgesCenter.size(): "
							<< edgesCenter.size()
							<< "; temp_edgesCenter.size(): "
							<< temp_edgesCenter.size() << endl;
					break;
				}
			}
			cout << "iterator: " << iterator << " ; mClusters.size(): "
					<< mClusters.size() << "; edgesCenter.size(): "
					<< edgesCenter.size() << "; temp_edgesCenter.size(): "
					<< temp_edgesCenter.size() << endl;
			iterator++;
		}
//		write into files
		cout << "After mClusters.size(): " << mClusters.size() << endl;
		ofstream output(input_file, ios::app);
		set<int>::iterator iter_set_center;
		int total_num = 0;
		for (iter_set_center = edgesCenter.begin();
				iter_set_center != edgesCenter.end(); iter_set_center++) {
			string line;
			int center = *iter_set_center;
			set<int>::iterator iter_set_edgeInCluster;
			for (iter_set_edgeInCluster = mClusters[center].begin();
					iter_set_edgeInCluster != mClusters[center].end();
					iter_set_edgeInCluster++) {
				total_num++;
				line += to_string(*iter_set_edgeInCluster) + ";";
			}
			line.pop_back();
			output << line << endl;
		}
		cout << "The total_num is: " << total_num << endl;
	} else {
		string line;
		int line_num = 0;
		while (getline(file, line)) {
			vector<string> line_vec;
			split(line, line_vec, ';');
			set<int> cEdges;
			int cEdgesNum = line_vec.size();
			for (int i = 0; i < cEdgesNum; i++) {
				int edgeId = stoi(line_vec[i]);
				cEdges.insert(edgeId);
			}
//			int edgeIdCenter = stoi(line_vec[0]);
			mClusters[line_num] = cEdges;
			line_num++;
		}
	}
	return mClusters;
}
